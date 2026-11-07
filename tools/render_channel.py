#!/usr/bin/env python
"""
render_channel.py -- visuals from a 3D compressible turbulent channel-flow DNS
(B-spline collocation solver; case M = 1.5, Re = 3000, periodic in x and z,
isothermal no-slip walls at y = 0 and y = 2).

Input (DATA_DIR) -- one of the solver's per-variable snapshot layouts:
  * <DATA_DIR>/{u,v,w}-parquet/{u,v,w}_<iter>.parquet   (columns x, y, z, value)
  * <DATA_DIR>/{u,v,w}_<iter>.out                       (ASCII: x y z value)
Points are listed x-fastest, then y, then z (Fortran order) on a tensor grid:
uniform in x (streamwise) and z (spanwise), wall-clustered in y (wall-normal).

Derived fields (computed here, not stored by the solver):
  Q = 1/2 (|Omega|^2 - |S|^2)   (2nd-order periodic differences in x, z;
                                 non-uniform 2nd-order differences in y)
  u' = u - <u>_{xz}(y)          (streamwise fluctuation about the plane mean)

Outputs: channel_hero.png, channel_panels.png, channel.gif

Usage: python render_channel.py [--data DIR] [--out DIR] [--only hero,panels,gif]
"""
import argparse
import glob
import os
import re

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage

# ----------------------------------------------------------------------------
DATA_DIR = os.environ.get("CHANNEL_DATA_DIR", "data/turbulence-data")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

CASE = "Compressible turbulent channel  ·  M = 1.5,  Re = 3000"
BG = (9, 14, 30)
BG_HEX = "#090e1e"
BOX_COLOR = "#6a7aa0"
TXT = (225, 230, 240)
SUB = (150, 160, 185)
CMAP_ISO = "RdBu_r"      # used when Q-isosurfaces are coloured by u'
ISO_COLORS = ("#3d5f9e", "#9fb6d9", "#f4ead2")   # Q-isosurfaces: wall distance 0 -> 1
CMAP_PLANE = "RdBu_r"    # wall-parallel planes coloured by u'
UP_CLIM = 0.25           # symmetric colour limit for u'
Q_LEVEL = 1.5            # isosurface level, in units of the snapshot's Q_rms
pv.OFF_SCREEN = True


# ----------------------------------------------------------------------------
# data
def list_iters(data_dir):
    pq = glob.glob(os.path.join(data_dir, "u-parquet", "u_*.parquet"))
    out = glob.glob(os.path.join(data_dir, "u_*.out"))
    files = pq or out
    its = sorted(int(re.search(r"_(\d+)\.", os.path.basename(f)).group(1)) for f in files)
    if not its:
        raise SystemExit(f"no u_<iter>.parquet / u_<iter>.out snapshots in {data_dir}")
    return its, ("parquet" if pq else "out")


def _read_var(data_dir, var, it, kind):
    if kind == "parquet":
        import pyarrow.parquet as pq
        t = pq.read_table(os.path.join(data_dir, f"{var}-parquet", f"{var}_{it:010d}.parquet"))
        a = np.column_stack([np.asarray(t.column(i), dtype=np.float64) for i in range(t.num_columns)])
    else:
        a = np.loadtxt(os.path.join(data_dir, f"{var}_{it:010d}.out"))
    a = a[:, -4:]                                    # x, y, z, value
    x, y, z = (np.unique(a[:, i]) for i in range(3))
    return x, y, z, a[:, 3].reshape(len(x), len(y), len(z), order="F")


def load(data_dir, it, kind, variables="uvw"):
    d = {}
    for v in variables:
        x, y, z, d[v] = _read_var(data_dir, v, it, kind)
    d.update(x=x, y=y, z=z)
    return d


def derive(d):
    """Q-criterion and streamwise fluctuation on the solver grid."""
    x, y, z = d["x"], d["y"], d["z"]
    dx, dz = x[1] - x[0], z[1] - z[0]

    def grads(f):
        return ((np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * dx),
                np.gradient(f, y, axis=1),
                (np.roll(f, -1, 2) - np.roll(f, 1, 2)) / (2 * dz))

    G = np.array([grads(d[c]) for c in "uvw"])       # G[i, j] = d u_i / d x_j
    S = 0.5 * (G + G.transpose(1, 0, 2, 3, 4))
    O = 0.5 * (G - G.transpose(1, 0, 2, 3, 4))
    d["Q"] = 0.5 * ((O ** 2).sum((0, 1)) - (S ** 2).sum((0, 1)))
    d["up"] = d["u"] - d["u"].mean(axis=(0, 2), keepdims=True)
    return d


def to_grid(d, up=2, sigma=0.6):
    """Periodically closed, lightly smoothed, upsampled RectilinearGrid.
    World axes: X = x (streamwise), Y = z (spanwise), Z = y (wall-normal)."""
    x, y, z = d["x"], d["y"], d["z"]
    Lx = x[1] - x[0] + x[-1] - x[0]
    Lz = z[1] - z[0] + z[-1] - z[0]
    xs = np.linspace(0, Lx, up * len(x) + 1)
    zs = np.linspace(0, Lz, up * len(z) + 1)
    # y: upsample the index coordinate so the wall stretching is preserved
    ys = np.interp(np.linspace(0, len(y) - 1, up * (len(y) - 1) + 1), np.arange(len(y)), y)
    g = pv.RectilinearGrid(xs, zs, ys)
    for name in ("Q", "up", "u"):
        a = d[name]
        if sigma > 0:
            a = ndimage.gaussian_filter(a, sigma, mode=("wrap", "nearest", "wrap"))
        a = np.pad(a, ((0, 1), (0, 0), (0, 1)), mode="wrap")
        fx = (len(xs) - 1) / (a.shape[0] - 1)
        fy = (len(ys) - 1) / (a.shape[1] - 1)
        fz = (len(zs) - 1) / (a.shape[2] - 1)
        a = ndimage.zoom(a, (fx, fy, fz), order=3, mode="nearest", grid_mode=False)
        a = a[:len(xs), :len(ys), :len(zs)]
        g.point_data[name] = a.transpose(1, 2, 0).ravel()
    return g


# ----------------------------------------------------------------------------
# rendering
def lights(p):
    p.add_light(pv.Light(position=(1.0, -0.6, 1.2), focal_point=(0, 0, 0),
                         intensity=0.95, light_type="camera light"))
    p.add_light(pv.Light(position=(-1.0, 0.3, 0.4), focal_point=(0, 0, 0),
                         intensity=0.35, light_type="camera light"))
    p.add_light(pv.Light(position=(0.0, 1.0, -0.6), focal_point=(0, 0, 0),
                         intensity=0.25, light_type="camera light"))


def iso(g, level, ymax=None):
    s = g.contour([level], scalars="Q")
    if ymax is not None and s.n_points:
        s = s.clip(normal="z", origin=(0, 0, ymax), invert=True)
    if s.n_points:
        s = s.compute_normals(consistent_normals=True, split_vertices=False)
    return s


def wall_plane(g, yval):
    """Wall-parallel slice at wall distance yval, carrying u'."""
    return g.slice(normal="z", origin=(0, 0, yval))


def render(g, size, cam, q_level, ymax=None, plane_y=None, plane_shift=0.0,
           box=True, ssao=True, parallel=False, iso_by="height"):
    p = pv.Plotter(off_screen=True, window_size=size, lighting="none")
    p.set_background(BG_HEX)
    lights(p)
    s = iso(g, q_level, ymax)
    if s.n_points:
        if iso_by == "up":
            kw = dict(scalars="up", cmap=CMAP_ISO, clim=(-UP_CLIM, UP_CLIM))
        else:
            s["wall distance"] = np.minimum(s.points[:, 2], 2.0 - s.points[:, 2])
            kw = dict(scalars="wall distance", cmap=iso_cmap(), clim=(0.0, 0.8))
        p.add_mesh(s, smooth_shading=True, ambient=0.22, diffuse=0.8, specular=0.45,
                   specular_power=40, show_scalar_bar=False, **kw)
    if plane_y is not None:
        pl = wall_plane(g, plane_y)
        pl.points[:, 2] -= plane_shift
        p.add_mesh(pl, scalars="up", cmap=CMAP_PLANE, clim=(-UP_CLIM, UP_CLIM),
                   ambient=0.55, diffuse=0.5, specular=0.0, show_scalar_bar=False)
    if box:
        b = g.outline() if ymax is None else pv.Box(bounds=(*g.bounds[:4], g.bounds[4], ymax)).outline()
        p.add_mesh(b, color=BOX_COLOR, line_width=1.2, opacity=0.8)
    p.camera_position = cam["pos"]
    p.camera.view_angle = cam.get("angle", 30.0)
    if parallel:
        p.enable_parallel_projection()
        p.camera.parallel_scale = cam["scale"]
    if ssao:
        p.enable_ssao(radius=0.12, bias=0.01, kernel_size=128, blur=True)
    p.enable_anti_aliasing("ssaa")
    img = p.screenshot(return_img=True)
    p.close()
    return Image.fromarray(img)


def iso_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("iso", ISO_COLORS)


def camera(center, dirn, dist, angle=30.0, up=(0, 0, 1)):
    d = np.asarray(dirn, float)
    d /= np.linalg.norm(d)
    c = np.asarray(center, float)
    return {"pos": [tuple(c + d * dist), tuple(c), up], "angle": angle}


# ----------------------------------------------------------------------------
# typography
def font(size, weight="regular"):
    from matplotlib import font_manager
    names = {"regular": "Lato-Regular", "light": "Lato-Light",
             "medium": "Lato-Medium", "bold": "Lato-Semibold"}[weight]
    for c in (f"/usr/share/fonts/truetype/lato/{names}.ttf",
              "/usr/share/fonts/truetype/lato/Lato-Regular.ttf",
              font_manager.findfont("DejaVu Sans")):
        if os.path.exists(c):
            return ImageFont.truetype(c, size)
    return ImageFont.load_default()


def math_font(size):
    from matplotlib import font_manager
    return ImageFont.truetype(font_manager.findfont("DejaVu Sans"), size)


def colorbar(width, label, fnt, clim=UP_CLIM, cmap=CMAP_ISO, bar_h=10, ticks=None):
    import matplotlib
    cm = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
    H = bar_h + 2 * fnt.size + 14
    im = Image.new("RGBA", (width, H), (0, 0, 0, 0))
    grad = (cm(np.linspace(0, 1, width))[:, :3] * 255).astype(np.uint8)
    im.paste(Image.fromarray(np.repeat(grad[None], bar_h, axis=0)), (0, 0))
    dr = ImageDraw.Draw(im)
    ty = bar_h + 4
    lo, mid, hi = ticks or (f"−{clim:g}", "0", f"+{clim:g}")
    dr.text((0, ty), lo, font=fnt, fill=SUB)
    if mid:
        dr.text((width // 2, ty), mid, font=fnt, fill=SUB, anchor="ma")
    dr.text((width, ty), hi, font=fnt, fill=SUB, anchor="ra")
    dr.text((width // 2, ty + fnt.size + 4), label, font=fnt, fill=SUB, anchor="ma")
    return im


def vignette(img, strength=0.3):
    W, H = img.size
    yy, xx = np.mgrid[0:H, 0:W]
    r = np.sqrt(((xx - W / 2) / (W / 2)) ** 2 + ((yy - H / 2) / (H / 2)) ** 2) / np.sqrt(2)
    f = 1 - strength * np.clip(r, 0, 1) ** 2.2
    a = np.asarray(img).astype(np.float32)
    bg = np.array(BG, np.float32)
    return Image.fromarray((bg + (a - bg) * f[..., None]).clip(0, 255).astype(np.uint8))


def grid_label(d):
    return f"{len(d['x'])} × {len(d['y'])} × {len(d['z'])} grid"


def legends(img, x0, y0, fnt, gap=40, width=200):
    """u' bar (plane) and wall-distance bar (structures), side by side."""
    a = colorbar(width, "near-wall plane: u′", fnt)
    b = colorbar(width, "vortex cores: wall distance y", fnt, cmap=iso_cmap(),
                 ticks=("0", "", "0.8"))
    img.paste(a, (x0, y0), a)
    img.paste(b, (x0 + width + gap, y0), b)


# ----------------------------------------------------------------------------
# deliverables
HERO_DIR = (0.55, -1.0, 0.55)
PLANE_Y = 0.03           # wall distance of the near-wall plane (drawn on the wall)


def hero_camera(g, dist=9.7):
    b = g.bounds
    c = ((b[0] + b[1]) / 2, (b[2] + b[3]) / 2 - 0.2, 0.24)
    return camera(c, HERO_DIR, dist, angle=30)


def make_hero(d, out, W=1600, H=900, ss=2):
    g = to_grid(d, up=3)
    qrms = np.sqrt(np.mean(d["Q"] ** 2))
    img = render(g, (W * ss, H * ss), hero_camera(g), Q_LEVEL * qrms, ymax=1.0,
                 plane_y=PLANE_Y, plane_shift=PLANE_Y)
    img = img.resize((W, H), Image.LANCZOS) if ss > 1 else img
    img = vignette(img, 0.3)
    dr = ImageDraw.Draw(img)
    dr.text((44, 36), CASE, font=font(32, "medium"), fill=TXT)
    dr.text((46, 82), f"DNS, B-spline collocation  ·  {grid_label(d)}  ·  lower half-channel",
            font=font(19, "regular"), fill=SUB)
    dr.text((46, 108), f"Q-criterion isosurfaces (Q = {Q_LEVEL:g} × rms(Q)) over near-wall streaks",
            font=math_font(16), fill=SUB)
    legends(img, W - 484, 46, math_font(14))
    img.save(out, optimize=True)
    print("wrote", out)


def _periodic_upsample(a, f, axes):
    """Cubic upsampling of a 2D slice, periodic along the listed axes."""
    pad = [(3, 3) if ax in axes else (0, 0) for ax in range(2)]
    b = np.pad(a, pad, mode="wrap")
    b = ndimage.zoom(b, f, order=3, mode="nearest")
    sl = tuple(slice(3 * f, -3 * f) if ax in axes else slice(None) for ax in range(2))
    return b[sl]


def make_panels(d, out, W=2400, H=800):
    """(a) Q-structures, top view  (b) near-wall streaks  (c) u' in a y-z cut."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    g = to_grid(d, up=3)
    qrms = np.sqrt(np.mean(d["Q"] ** 2))
    x, y, z, upf = d["x"], d["y"], d["z"], d["up"]
    Lx = x[-1] + (x[1] - x[0])
    Lz = z[-1] + (z[1] - z[0])

    top, bot = 100, 110
    PW = W // 3
    PH = H - top - bot
    canvas = Image.new("RGB", (W, H), BG)

    # common plotting box (pixels inside each panel) for the x-z panels
    mx = 70
    bw = PW - 2 * mx                      # box width
    bh = int(round(bw * Lz / Lx))         # physical aspect of x-z plane
    by = (PH - bh) // 2

    # (a) top view of Q-structures, rendered at exactly the x-z box aspect
    ss = 2
    b = g.bounds
    c = ((b[0] + b[1]) / 2, (b[2] + b[3]) / 2, 0.5)
    cam = camera(c, (0, 0, 1), 10, up=(0, 1, 0))
    cam["scale"] = Lz / 2
    a_img = render(g, (bw * ss, bh * ss), cam, Q_LEVEL * qrms, ymax=1.0, parallel=True,
                   box=False, ssao=True, plane_y=PLANE_Y, plane_shift=PLANE_Y)
    a_img = a_img.resize((bw, bh), Image.LANCZOS)

    fam = [f.name for f in font_manager.fontManager.ttflist if f.name == "Lato"]
    plt.rcParams.update({"font.family": fam[0] if fam else "DejaVu Sans", "font.size": 13,
                         "text.color": "#e1e6f0", "axes.labelcolor": "#969fb9",
                         "xtick.color": "#969fb9", "ytick.color": "#969fb9",
                         "axes.edgecolor": "#3a4566", "axes.linewidth": 0.8})
    bgc = np.array(BG) / 255

    def panel(img_arr, extent, xlabel, ylabel, box, interp="bicubic", raster=False):
        fig = plt.figure(figsize=(PW / 100, PH / 100), dpi=100, facecolor=bgc)
        l, btm, w_, h_ = box
        ax = fig.add_axes([l / PW, btm / PH, w_ / PW, h_ / PH], facecolor=bgc)
        if raster:
            ax.imshow(img_arr, extent=extent, origin="upper", aspect="auto", interpolation="lanczos")
        else:
            ax.imshow(img_arr, extent=extent, origin="lower", aspect="auto", cmap=CMAP_PLANE,
                      vmin=-UP_CLIM, vmax=UP_CLIM, interpolation=interp)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(length=3)
        fig.canvas.draw()
        im = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3])
        plt.close(fig)
        return im

    xz_box = (mx, PH - by - bh, bw, bh)
    canvas.paste(panel(np.asarray(a_img), (0, Lx, 0, Lz), "x  (streamwise)", "z  (spanwise)",
                       xz_box, raster=True), (0, top))

    # (b) wall-parallel plane closest to y = 0.05
    j = int(np.argmin(np.abs(y - 0.05)))
    sl = _periodic_upsample(upf[:, j, :].T, 6, axes=(0, 1))
    canvas.paste(panel(sl, (0, Lx, 0, Lz), "x  (streamwise)", "z  (spanwise)", xz_box), (PW, top))

    # (c) cross-stream (z-y) plane at mid-length; resample y onto a uniform axis
    i = len(x) // 2
    sl = upf[i]                                   # (ny, nz)
    yu = np.linspace(0, 2, 4 * len(y))
    col = np.array([np.interp(yu, y, sl[:, k]) for k in range(sl.shape[1])]).T
    col = _periodic_upsample(col, 4, axes=(1,))[:, :]
    ch = int(round(bw * 2.0 / Lz))
    canvas.paste(panel(col, (0, Lz, 0, 2), "z  (spanwise)", "y  (wall-normal)",
                       (mx, (PH - ch) // 2, bw, ch)), (2 * PW, top))

    dr = ImageDraw.Draw(canvas)
    dr.text((44, 30), CASE, font=font(30, "medium"), fill=TXT)
    dr.text((W - 44, 40), f"DNS  ·  {grid_label(d)}  ·  colour: streamwise velocity fluctuation u′ "
            "(plane-averaged mean removed)", font=font(18), fill=SUB, anchor="ra")
    caps = [f"(a)  vortex cores, Q = {Q_LEVEL:g} × rms(Q), lower half (top view)",
            f"(b)  near-wall streaks, wall-parallel plane y = {y[j]:.3f}",
            f"(c)  cross-stream plane x = {x[i]:.2f}"]
    for k, cap in enumerate(caps):
        dr.text((k * PW + PW // 2, H - bot + 18), cap, font=math_font(17), fill=TXT, anchor="ma")
        if k:
            dr.line([(k * PW, top + 20), (k * PW, H - bot - 10)], fill=(28, 36, 58), width=1)
    legends(canvas, W // 2 - 220, H - 56, math_font(13))
    canvas.save(out, optimize=True)
    print("wrote", out)


def make_gif(data_dir, its, kind, out, W=720, H=410, fps=12, max_mb=8.0, colors=160):
    """Fixed absolute isosurface level (Q_LEVEL x rms(Q) of the last frame) so the
    growth and decay of vortical activity through transition is visible."""
    frames = []
    fnt, fm = font(18, "medium"), math_font(13)
    d_last = derive(load(data_dir, its[-1], kind))
    level = Q_LEVEL * np.sqrt(np.mean(d_last["Q"] ** 2))
    for n, it in enumerate(its):
        d = derive(load(data_dir, it, kind))
        g = to_grid(d, up=2)
        img = render(g, (W * 2, H * 2), hero_camera(g), level, ymax=1.0,
                     plane_y=PLANE_Y, plane_shift=PLANE_Y, ssao=False)
        img = vignette(img.resize((W, H), Image.LANCZOS), 0.3)
        dr = ImageDraw.Draw(img)
        dr.text((20, 16), CASE, font=fnt, fill=TXT)
        dr.text((20, 42), f"Q = {level:.2f} isosurfaces over near-wall u′  ·  {grid_label(d)}",
                font=fm, fill=SUB)
        # iteration progress bar
        x0, x1, yb = 20, W - 20, H - 18
        frac = (it - its[0]) / max(its[-1] - its[0], 1)
        dr.line([(x0, yb), (x1, yb)], fill=(40, 50, 78), width=3)
        dr.line([(x0, yb), (x0 + frac * (x1 - x0), yb)], fill=(214, 96, 77), width=3)
        dr.text((x1, yb - 8), f"iteration {it:,}", font=fm, fill=SUB, anchor="rb")
        frames.append(img)
        print(f"  frame {n + 1}/{len(its)}  it={it}", flush=True)
    strip = Image.new("RGB", (W, H * 5))
    for k, idx in enumerate(np.linspace(0, len(frames) - 1, 5).astype(int)):
        strip.paste(frames[idx], (0, k * H))
    pal = strip.quantize(colors=colors, method=Image.Quantize.MEDIANCUT)
    q = [f.quantize(palette=pal, dither=Image.Dither.FLOYDSTEINBERG) for f in frames]
    dur = [int(1000 / fps)] * len(q)
    dur[-1] = 1200
    q[0].save(out, save_all=True, append_images=q[1:], duration=dur, loop=0, optimize=True)
    mb = os.path.getsize(out) / 1e6
    print("wrote", out, f"{mb:.2f} MB")
    if mb > max_mb:
        print(f"WARNING: GIF exceeds {max_mb} MB; reduce frames or size")


# ----------------------------------------------------------------------------
# hero v2: Q computed on a cubic-upsampled velocity field
V2_UP = 3                 # upsampling factor before differentiation
V2_Q = 0.55               # isosurface level, units of rms(Q) on the fine grid
V2_SIGMA = 1.2            # Gaussian smoothing (fine-grid cells) before differentiation
V2_CMAP = "RdBu_r"
V2_UCLIM = (0.35, 1.15)   # colour range for streamwise velocity u on the cores


def fine_fields(d, up=V2_UP, sigma=V2_SIGMA):
    from scipy.interpolate import CubicSpline
    x, y, z = d["x"], d["y"], d["z"]
    nx, ny, nz = len(x), len(y), len(z)
    Lx = nx * (x[1] - x[0])
    Lz = nz * (z[1] - z[0])
    xf = np.arange(nx * up + 1) * Lx / (nx * up)
    zf = np.arange(nz * up + 1) * Lz / (nz * up)
    sf = np.linspace(0, ny - 1, up * (ny - 1) + 1)
    yf = CubicSpline(np.arange(ny), y)(sf)

    def fine(a):
        # periodic cubic spline in x and z (closed with the first plane), spline in index-y
        ax = np.concatenate([a, a[:1]], axis=0)
        a1 = CubicSpline(np.r_[x, Lx], ax, axis=0, bc_type="periodic")(xf)
        az = np.concatenate([a1, a1[:, :, :1]], axis=2)
        a2 = CubicSpline(np.r_[z, Lz], az, axis=2, bc_type="periodic")(zf)
        a3 = CubicSpline(np.arange(ny), a2, axis=1)(sf)
        if sigma > 0:
            a3 = ndimage.gaussian_filter(a3, sigma, mode=("wrap", "nearest", "wrap"))
        return a3

    F = {c: fine(d[c]) for c in "uvw"}
    dxf, dzf = xf[1] - xf[0], zf[1] - zf[0]

    def grads(f):
        g0 = np.gradient(f, dxf, axis=0)
        g2 = np.gradient(f, dzf, axis=2)
        return g0, np.gradient(f, yf, axis=1), g2

    G = np.array([grads(F[c]) for c in "uvw"])
    S = 0.5 * (G + G.transpose(1, 0, 2, 3, 4))
    O = 0.5 * (G - G.transpose(1, 0, 2, 3, 4))
    Q = 0.5 * ((O ** 2).sum((0, 1)) - (S ** 2).sum((0, 1)))
    del G, S, O
    g = pv.RectilinearGrid(xf, zf, yf)
    g.point_data["Q"] = Q.transpose(1, 2, 0).ravel()
    g.point_data["u"] = F["u"].transpose(1, 2, 0).ravel()
    return g, float(np.sqrt(np.mean(Q[1:-1, 2:-2, 1:-1] ** 2)))


def render_v2(g, size, cam, level, ymin=None, ymax=None, box=True):
    p = pv.Plotter(off_screen=True, window_size=size, lighting="none")
    p.set_background(BG_HEX)
    lights(p)
    s = g.contour([level], scalars="Q")
    if ymax is not None:
        s = s.clip(normal="z", origin=(0, 0, ymax), invert=True)
    if ymin is not None:
        s = s.clip(normal="z", origin=(0, 0, ymin), invert=False)
    s = s.connectivity("all")
    rid = np.asarray(s.cell_data["RegionId"])
    keep = np.where(np.bincount(rid) >= 120)[0]          # drop speckles (< 120 triangles)
    s = s.extract_cells(np.isin(rid, keep)).extract_surface(algorithm="dataset_surface")
    s = s.compute_normals(consistent_normals=True, split_vertices=False)
    p.add_mesh(s, scalars="u", cmap=V2_CMAP, clim=V2_UCLIM, smooth_shading=True,
               ambient=0.2, diffuse=0.78, specular=0.55, specular_power=45,
               show_scalar_bar=False)
    if box:
        b = g.bounds
        p.add_mesh(pv.Box(bounds=b).outline(), color=BOX_COLOR, line_width=1.0, opacity=0.6)
    p.camera_position = cam["pos"]
    p.camera.view_angle = cam.get("angle", 30.0)
    p.enable_ssao(radius=0.15, bias=0.01, kernel_size=128, blur=True)
    p.enable_anti_aliasing("ssaa")
    img = Image.fromarray(p.screenshot(return_img=True))
    p.close()
    return img


V2_CAM = dict(center=(3.6, 2.05, 0.95), dirn=(-1.0, -0.42, 0.30), dist=5.6, angle=40)


def make_hero_v2(d, out, out_wide, W=1600, H=900, ss=2, cam=None, level=V2_Q):
    g, qrms = fine_fields(d)
    c = cam or V2_CAM
    img = render_v2(g, (W * ss, H * ss), camera(c["center"], c["dirn"], c["dist"],
                                                angle=c["angle"]), level * qrms)
    img = vignette(img.resize((W, H), Image.LANCZOS), 0.32)
    dr = ImageDraw.Draw(img)
    dr.text((44, 36), "Compressible turbulent channel, M = 1.5", font=font(30, "medium"), fill=TXT)
    img.save(out, optimize=True)
    print("wrote", out)
    if out_wide:
        Hw = 560
        top = int(0.5 * (H - Hw)) + 30
        wide = img.crop((0, top, W, top + Hw))
        wide = vignette(wide, 0.0)
        # redraw title in the crop
        base = img.crop((0, top, W, top + Hw))
        wide = Image.fromarray(np.asarray(base))
        if top > 36:
            ImageDraw.Draw(wide).text((44, 30), "Compressible turbulent channel, M = 1.5",
                                      font=font(28, "medium"), fill=TXT)
        wide.save(out_wide, optimize=True)
        print("wrote", out_wide)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA_DIR)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--only", default="hero,hero2,panels,gif")
    ap.add_argument("--iter", type=int, default=None, help="snapshot for still images (default: last)")
    ap.add_argument("--gif-frames", type=int, default=116)
    args = ap.parse_args()
    its, kind = list_iters(args.data)
    os.makedirs(args.out, exist_ok=True)
    todo = set(args.only.split(","))
    it = args.iter if args.iter is not None else its[-1]
    if todo & {"hero", "hero2", "panels"}:
        d = derive(load(args.data, it, kind))
        print(f"snapshot it={it}: {grid_label(d)}, Q_rms={np.sqrt(np.mean(d['Q']**2)):.3g}, "
              f"u' rms={d['up'].std():.3g}")
        if "hero" in todo:
            make_hero(d, os.path.join(args.out, "channel_hero.png"))
        if "hero2" in todo:
            make_hero_v2(d, os.path.join(args.out, "channel_hero_v2.png"),
                         os.path.join(args.out, "channel_hero_v2_wide.png"))
        if "panels" in todo:
            make_panels(d, os.path.join(args.out, "channel_panels.png"))
    if "gif" in todo and len(its) > 1:
        sel = [its[k] for k in np.linspace(0, len(its) - 1, min(args.gif_frames, len(its))).astype(int)]
        make_gif(args.data, sel, kind, os.path.join(args.out, "channel.gif"))


if __name__ == "__main__":
    main()
