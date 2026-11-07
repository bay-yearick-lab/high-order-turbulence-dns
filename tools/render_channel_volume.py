#!/usr/bin/env python
"""
render_channel_volume.py -- volume-rendered animation of the compressible
turbulent channel DNS (M = 1.5, Re = 3000, 55 x 80 x 55 B-spline collocation grid,
periodic in x and z, isothermal no-slip walls at y = 0 and y = 2).

Instead of isosurfaces (which look lumpy on a coarse grid), each snapshot is
  1. upsampled 3x with cubic interpolation (periodic in x/z, cubic spline onto a
     uniform wall-normal grid in y), lightly Gaussian-smoothed,
  2. differentiated to the fluctuating vorticity magnitude |omega'| (vorticity of
     u' = u - <u>_xz(y), v, w -- so the laminar mean shear does not light up the walls),
  3. volume-rendered with a transfer function that is transparent for weak vorticity
     and glows for strong vortical structures; a faint u' plane sits on the lower wall.

The animation sweeps laminar start -> transition burst -> sustained turbulence with a
slow camera orbit. No text overlays.

Outputs (OUT_DIR): card_channel_volume.gif (400x400), channel_volume.gif (900x500)

Usage: python render_channel_volume.py [--data DIR] [--out DIR] [--still IDX]
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pyvista as pv
from scipy import ndimage
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from render_channel import list_iters, load  # noqa: E402

DATA_DIR = os.environ.get("CHANNEL_DATA_DIR", "/tmp/claude-1000/phdviz/turb/turbulence-data")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
FFMPEG = os.environ.get("FFMPEG", "ffmpeg")
BG_HEX = "#0B1626"
UP = 3              # upsampling factor
SIGMA = 1.1         # Gaussian smoothing on the fine grid (cells)
NFRAMES = 72        # 7.2 s at 10 fps
FPS = 10
ORBIT = 22.0        # total camera azimuth sweep (deg)
pv.OFF_SCREEN = True


# ----------------------------------------------------------------------------
def fine_fields(d):
    x, y, z = d["x"], d["y"], d["z"]
    nx, ny, nz = len(x), len(y), len(z)
    dx, dz = x[1] - x[0], z[1] - z[0]
    ny_f = UP * (ny - 1) + 1
    yf = np.linspace(y[0], y[-1], ny_f)
    dxf, dzf, dyf = dx / UP, dz / UP, yf[1] - yf[0]

    def fine(a):
        a = ndimage.zoom(a, (UP, 1, UP), order=3, mode="grid-wrap", grid_mode=True)
        a = CubicSpline(y, a, axis=1)(yf)
        return ndimage.gaussian_filter(a, SIGMA, mode=("wrap", "nearest", "wrap"))

    ubar = d["u"].mean(axis=(0, 2), keepdims=True)
    F = [fine(d["u"] - ubar), fine(d["v"]), fine(d["w"])]
    Ubar = CubicSpline(y, ubar[0, :, 0])(yf)[None, :, None]

    def ddx(f):
        return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * dxf)

    def ddz(f):
        return (np.roll(f, -1, 2) - np.roll(f, 1, 2)) / (2 * dzf)

    def ddy(f):
        return np.gradient(f, dyf, axis=1)

    uf, vf, wf = F
    wx = ddy(wf) - ddz(vf)
    wy = ddz(uf) - ddx(wf)
    wz = ddx(vf) - ddy(uf)
    om = np.sqrt(wx ** 2 + wy ** 2 + wz ** 2)
    # Q-criterion of the full velocity (mean shear dU/dy added to du/dy)
    dUdy = np.gradient(Ubar, dyf, axis=1)
    G = [[ddx(uf), ddy(uf) + dUdy, ddz(uf)], [ddx(vf), ddy(vf), ddz(vf)], [ddx(wf), ddy(wf), ddz(wf)]]
    Q = np.zeros_like(om)
    for i in range(3):
        for j in range(3):
            Q -= 0.5 * G[i][j] * G[j][i]       # Q = -1/2 du_i/dx_j du_j/dx_i (incompressible form)
    return {"om": om, "Q": Q, "uf": uf}, (dxf, dyf, dzf)


def to_image(a, spacing):
    """(nx, ny, nz) solver-ordered array -> ImageData with world X=x, Y=z, Z=y."""
    nx, ny, nz = a.shape
    a = np.concatenate([a, a[:1]], axis=0)
    a = np.concatenate([a, a[:, :, :1]], axis=2)       # periodic closure
    g = pv.ImageData(dimensions=(nx + 1, nz + 1, ny), spacing=(spacing[0], spacing[2], spacing[1]))
    g.point_data["s"] = a.transpose(1, 2, 0).ravel()  # x fastest, then z, then y
    return g


# ----------------------------------------------------------------------------
def opacity_tf(n=64, lo=0.16, gamma=1.7, top=0.9):
    import vtk
    f = vtk.vtkPiecewiseFunction()
    for t in np.linspace(0, 1, n):
        f.AddPoint(float(t), float(top * np.clip((t - lo) / (1 - lo), 0, 1) ** gamma))
    return f


def add_volume(p, vol, unit=0.05, **tf):
    a = p.add_volume(vol, scalars="s", cmap="inferno", clim=(0.0, 1.0), opacity="linear",
                     shade=False, mapper="gpu", blending="composite", show_scalar_bar=False)
    prop = a.GetProperty()
    prop.SetScalarOpacity(opacity_tf(**tf))
    prop.SetScalarOpacityUnitDistance(unit)
    prop.SetInterpolationTypeToLinear()
    m = a.GetMapper()
    m.SetAutoAdjustSampleDistances(False)
    m.SetSampleDistance(0.004)
    if hasattr(m, "SetUseJittering"):
        m.SetUseJittering(True)
    return a


def plane_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("floor", ["#0d1a30", "#1c3a63", "#3f78ad", "#a9cdee"])


def scalar(F, ref):
    """Volume scalar in [0, 1]: positive Q normalised by the turbulent reference rms(Q+)(y)
    (fixed in time, so the laminar start stays dark), soft-compressed for continuous cores."""
    q = np.clip(F["Q"], 0, None) / (ref["k"] * ref["qrms"])
    return np.clip(q, 0, 1) ** ref["pow"]


def render(F, sp, size, azim, ref, plane=True):
    W, H = size
    p = pv.Plotter(off_screen=True, window_size=(W, H), lighting="none")
    p.set_background(BG_HEX)
    vol = to_image(scalar(F, ref), sp)
    Lx, Lz, Ly = [b for b in vol.bounds[1::2]]
    add_volume(p, vol, unit=ref["unit"], lo=ref["lo"], gamma=ref["gamma"], top=ref["top"])
    if plane:   # lower-wall u' plane (y ~ 0.05): near-wall streaks as a faint floor
        j = int(round(0.05 / sp[1]))
        pl = to_image(F["uf"][:, j:j + 1, :], sp).extract_surface(algorithm="dataset_surface")
        pl.points[:, 2] = 0.0
        p.add_mesh(pl, scalars="s", cmap=plane_cmap(), clim=(-0.25, 0.25), opacity=ref["plane_op"],
                   show_scalar_bar=False, lighting=False)
    box = pv.Box(bounds=(0, Lx, 0, Lz, 0, Ly)).extract_feature_edges(5)
    p.add_mesh(box, color="#5f7396", opacity=0.35, line_width=1.0, lighting=False)
    c = np.array([Lx / 2, Lz / 2, Ly / 2])
    a = np.radians(azim)
    el = np.radians(ref["elev"])
    dist = ref["dist_wide"] if W > H else ref["dist_sq"]
    dirn = np.array([np.cos(el) * np.cos(a), np.cos(el) * np.sin(a), np.sin(el)])
    p.camera.position = tuple(c + dist * dirn)
    p.camera.focal_point = tuple(c)
    p.camera.up = (0, 0, 1)
    p.camera.view_angle = 30
    img = p.screenshot(return_img=True)
    p.close()
    return img


STYLE = dict(k=7.0, pow=1.0, unit=0.18, lo=0.08, gamma=1.2, top=0.9, plane_op=0.45,
             elev=26.0, dist_wide=12.3, dist_sq=15.5, azim0=-128.0)


def reference(data_dir, its, kind):
    F, _ = fine_fields(load(data_dir, its[int(0.8 * len(its))], kind))
    q = np.clip(F["Q"], 0, None)
    qrms = np.sqrt((q ** 2).mean(axis=(0, 2), keepdims=True))
    qrms = ndimage.gaussian_filter1d(qrms, 4, axis=1, mode="nearest")
    return dict(STYLE, qrms=qrms)


def frame_plan(its):
    # snapshots: spend a few frames on the laminar start, more on the burst, then turbulence
    n = len(its)
    s = np.linspace(0, 1, NFRAMES)
    idx = np.round(s * (n - 1)).astype(int)
    azim = ORBIT * (0.5 - 0.5 * np.cos(np.pi * s))
    return idx, azim


def encode(frames_dir, out, size, fps, max_mb):
    W, H = size
    vf = f"scale={W}:{H}:flags=lanczos,split[a][b];[a]palettegen=max_colors=256:stats_mode=full[p];[b][p]paletteuse=dither=none"
    subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", os.path.join(frames_dir, "f%03d.png"), "-vf", vf, "-loop", "0", out],
                   check=True)
    mb = os.path.getsize(out) / 1e6
    print(f"{out}: {mb:.2f} MB")
    if mb > max_mb:
        print("WARNING: over size budget")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DATA_DIR)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--still", type=int, nargs="*", help="render stills for these snapshot indices")
    ap.add_argument("--set", nargs="*", help="override STYLE entries, e.g. lo=0.2")
    ap.add_argument("--only", choices=["card", "banner"], help="render just one of the two GIFs")
    ap.add_argument("--ss", type=int, default=2, help="supersampling factor")
    args = ap.parse_args()
    its, kind = list_iters(args.data)
    ref = reference(args.data, its, kind)
    for kv in args.set or []:
        k_, v_ = kv.split("=")
        ref[k_] = float(v_)
    sizes = {"card": (400, 400), "banner": (900, 500)}
    if args.only:
        sizes = {args.only: sizes[args.only]}
    from PIL import Image
    if args.still is not None:
        for i in args.still:
            F, sp = fine_fields(load(args.data, its[i], kind))
            for k, (W, H) in sizes.items():
                img = render(F, sp, (W * args.ss, H * args.ss), ref["azim0"] + ORBIT / 2, ref)
                Image.fromarray(img).resize((W, H), Image.LANCZOS).save(
                    os.path.join(args.out, f"_still_{k}_{i:03d}.png"))
        return
    idx, azim = frame_plan(its)
    tmp = {k: tempfile.mkdtemp(prefix=f"chvol_{k}_") for k in sizes}
    for f, (i, az) in enumerate(zip(idx, azim)):
        F, sp = fine_fields(load(args.data, its[i], kind))
        for k, (W, H) in sizes.items():
            img = render(F, sp, (W * args.ss, H * args.ss), ref["azim0"] + az, ref)
            Image.fromarray(img).resize((W, H), Image.LANCZOS).save(os.path.join(tmp[k], f"f{f:03d}.png"))
        print(f"frame {f + 1}/{len(idx)}", flush=True)
    if "card" in sizes:
        encode(tmp["card"], os.path.join(args.out, "card_channel_volume.gif"), sizes["card"], FPS, 5.0)
    if "banner" in sizes:
        encode(tmp["banner"], os.path.join(args.out, "channel_volume.gif"), sizes["banner"], FPS, 8.0)
    for t in tmp.values():
        shutil.rmtree(t)


if __name__ == "__main__":
    main()
