"""Regenerate the README figures in assets/.

    python tools/make_readme_assets.py

Inputs (all in the repository):
  assets/channel/channel_render.png   channel DNS render (tools/render_channel.py)
  data/channel_spectra/*.csv           1-D streamwise energy spectra, turbulent channel
  python/bspline_filter                2-D Burgers errors are recomputed here
Numbers for the robustness, cost and MMS figures are the published tables of the
paper and dissertation, listed below.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
STYLE = ROOT / "tools" / "lab.mplstyle"
sys.path.insert(0, str(ROOT / "python"))

INK, TEAL, AMBER, SLATE, CORAL, GRAY = "#1B2A41", "#0F9D8A", "#F2A541", "#3D5A80", "#E07A5F", "#8D99AE"

# Paper Table 11 / dissertation Table 14.5: time to floating-point exception, Re = 3000, t in [0, 100].
BLOWUP = [("Optimized 13pt-II", 4.13), ("Optimized 13pt-I", 4.55), ("Optimized 11pt-I", 6.21),
          ("Compact, 8th-order boundary", 8.11), ("Optimized 11pt-II", 8.46),
          ("Compact, 6th-order boundary", 10.28), ("Compact, 4th-order boundary", 11.43),
          ("Optimized 9pt-I", 13.34), ("Optimized 9pt-II", 14.37)]
# Paper Table 12 / dissertation Table 14.6: time per Runge-Kutta step, 60 x 80 x 60 grid.
COST = [("Unfiltered", 1.731), ("Constant\n$\\lambda$", 1.735), ("Optimized\n$\\lambda$", 1.881),
        ("Optimized\n$\\lambda$ + $\\vec{\\lambda}$", 2.110)]
# Dissertation Table 14.3: MMS L-infinity errors, degree-7 B-splines, 3-D compressible Navier-Stokes.
MMS_N = np.array([40, 60, 80, 100, 120, 140, 160, 180])
MMS = {
    "$\\rho$": [3.58e-1, 1.45e-2, 1.88e-3, 4.15e-4, 1.27e-4, 4.85e-5, 2.15e-5, 1.05e-5],
    "$T$": [3.21e0, 1.73e-1, 2.20e-2, 4.76e-3, 1.41e-3, 5.37e-4, 2.38e-4, 1.16e-4],
    "$u_1$, $u_3$": [1.14e0, 4.56e-2, 5.49e-3, 1.16e-3, 3.44e-4, 1.32e-4, 5.92e-5, 2.91e-5],
    "$u_2$": [1.81e0, 1.03e-1, 1.34e-2, 2.93e-3, 8.74e-4, 3.39e-4, 1.51e-4, 7.42e-5],
}


def _font(size: int, bold: bool = False):
    names = ["Lato-Bold.ttf", "DejaVuSans-Bold.ttf"] if bold else ["Lato-Regular.ttf", "DejaVuSans.ttf"]
    for d in ("/usr/share/fonts/truetype/lato", "/usr/share/fonts/truetype/dejavu"):
        for n in names:
            p = Path(d) / n
            if p.exists():
                return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def hero() -> None:
    """1600 x 560 banner: channel DNS render with title and subtitle."""
    W, H, bg = 1600, 560, (9, 14, 30)
    src = Image.open(ASSETS / "channel" / "channel_render.png").convert("RGB")
    scene = src.crop((140, 145, 1600, 845))
    sw = 1000
    scene = scene.resize((sw, int(scene.height * sw / scene.width)), Image.LANCZOS)
    canvas = Image.new("RGB", (W, H), bg)
    x0, y0 = W - sw - 10, (H - scene.height) // 2 + 10
    canvas.paste(scene, (x0, y0))
    d = ImageDraw.Draw(canvas)
    d.text((56, 62), "high-order-turbulence-dns", font=_font(54, bold=True), fill=(236, 240, 246))
    for i, line in enumerate(["High-order B-spline DNS of compressible",
                              "turbulence with provably stable,",
                              "boundary-consistent filtering and",
                              "GCV-selected filter strength"]):
        d.text((58, 140 + 38 * i), line, font=_font(28), fill=(200, 208, 222))
    for i, line in enumerate(["Channel flow, M = 1.5, Re = 3000",
                              "55 x 80 x 55 grid, iteration 138,000"]):
        d.text((58, H - 100 + 30 * i), line, font=_font(21), fill=(150, 162, 186))
    d.line([(58, 312), (330, 312)], fill=(15, 157, 138), width=4)
    canvas.save(ASSETS / "hero.png", optimize=True)


def card_robustness() -> None:
    """Time to divergence at Re = 3000 versus the B-spline filter."""
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    names = [n for n, _ in BLOWUP] + ["Compact, 2nd-order centered", "B-spline filter (any strength)"]
    times = [t for _, t in BLOWUP] + [100.0, 100.0]
    colors = [AMBER] * len(BLOWUP) + [GRAY, TEAL]
    y = np.arange(len(names))
    ax.barh(y, times, color=colors, height=0.68)
    for yi, t, c in zip(y, times, colors):
        label = f"diverged at t = {t:.2f}" if t < 100 else "stable to t = 100"
        ax.text(t + 1.2, yi, label, va="center", fontsize=10, color=INK)
    ax.set_yticks(y, names, fontsize=10)
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("simulated time before floating-point exception")
    ax.grid(axis="y", visible=False)
    ax.set_title("Channel flow at Re = 3000: B-spline filter stays stable,\n"
                 "9 of 10 finite-difference filters diverge", loc="left")
    fig.text(0.01, -0.035, "Degree-7 B-splines, 60 x 80 x 60 grid, M = 1.5. Gray: stable but 2nd-order "
             "accurate at the wall.", fontsize=8.5, color=GRAY)
    fig.savefig(ASSETS / "card_robustness.png")
    plt.close(fig)


def card_accuracy_cost() -> None:
    """2-D Burgers error (recomputed) and cost per Runge-Kutta step (published timing)."""
    from bspline_filter import burgers2d

    cases = ["unfiltered", "optimized-scalar", "optimized-vector"]
    errs = [burgers2d.solve(c)["linf_u"] for c in cases]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.4), gridspec_kw={"wspace": 0.38})
    a.bar([0, 1, 2], errs, color=[AMBER, SLATE, TEAL], width=0.62)
    a.set_yscale("log")
    a.set_ylim(5e-3, 6e-1)
    a.set_xticks([0, 1, 2], ["Unfiltered", "Optimized\n$\\lambda$", "Optimized\n$\\vec{\\lambda}$"])
    for i, e in enumerate(errs):
        a.text(i, e * 1.25, f"{e:.2e}", ha="center", fontsize=10.5, color=INK)
    a.set_ylabel("$L_\\infty$ error at t = 1")
    a.set_title("(a) 2-D Burgers, 20 x 20, degree 7", fontsize=12, loc="left")
    a.grid(axis="x", visible=False)
    base = COST[0][1]
    rel = [(t / base - 1.0) * 100 for _, t in COST]
    b.bar(range(4), [t for _, t in COST], color=[AMBER, SLATE, SLATE, TEAL], width=0.62)
    for i, (r, (_, t)) in enumerate(zip(rel, COST)):
        b.text(i, t + 0.04, "baseline" if i == 0 else f"+{r:.1f}%", ha="center", fontsize=10.5,
               color=INK)
    b.set_xticks(range(4), [n for n, _ in COST], fontsize=10)
    b.set_ylim(0, 2.5)
    b.set_ylabel("seconds per Runge-Kutta step")
    b.set_title("(b) channel DNS, 60 x 80 x 60", fontsize=12, loc="left")
    b.grid(axis="x", visible=False)
    fig.suptitle(f"{errs[0] / errs[2]:.0f}x lower error; filter adds 0.2% to 22% per step",
                 x=0.06, ha="left", fontsize=15, fontweight="bold", color=INK, y=1.02)
    fig.savefig(ASSETS / "card_accuracy_cost.png")
    plt.close(fig)


def card_mms() -> None:
    """Sixth-order convergence of the filtered 3-D compressible solver (MMS)."""
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    markers = ["o", "s", "^", "D"]
    for (name, e), m, c in zip(MMS.items(), markers, [TEAL, SLATE, CORAL, AMBER]):
        ax.loglog(MMS_N, e, marker=m, color=c, lw=2, ms=6, label=name)
    ref = 1.6e-5 * (MMS_N / 180.0) ** -6
    ax.loglog(MMS_N[2:], ref[2:], "--", color=INK, lw=1.4, label="slope 6")
    ax.set_xlabel("points per direction (N x N x N)")
    ax.set_ylabel("$L_\\infty$ error")
    ax.set_xticks(MMS_N[::1], [str(n) for n in MMS_N])
    ax.minorticks_off()
    ax.legend(ncol=2, loc="lower left")
    order = np.log(np.array(MMS["$\\rho$"][-2]) / MMS["$\\rho$"][-1]) / np.log(180 / 160)
    ax.set_title(f"Sixth-order convergence of the filtered 3-D compressible\n"
                 f"Navier-Stokes solver (observed {order:.2f} for density)", loc="left")
    fig.text(0.01, -0.035, "Manufactured solutions, degree-7 B-splines, 6th-derivative penalty, "
             "GCV-optimized vector strength (wall-normal) and scalar strengths (periodic).", fontsize=8.5, color=GRAY)
    fig.savefig(ASSETS / "card_mms.png")
    plt.close(fig)


def card_spectra() -> None:
    """Streamwise energy spectra of the under-resolved channel versus reference DNS."""
    d = ROOT / "data" / "channel_spectra"
    scale = 3.26239  # normalization used in the published figure
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    spec = [("unfiltered", "Unfiltered", AMBER, "--"),
            ("constant_scalar", "Constant $\\lambda^2 = 10^{-6}$", SLATE, "--"),
            ("optimized_scalar", "Optimized $\\lambda$", CORAL, "--"),
            ("optimized_vector", "Optimized $\\vec{\\lambda}$", TEAL, "-")]
    k, e = np.loadtxt(d / "streamwise_reference_dns.csv", unpack=True)
    ax.loglog(k, scale * e, color=INK, lw=2.2, label="Reference DNS (120 x 180 x 120)")
    for f, lab, c, ls in spec:
        k, e = np.loadtxt(d / f"streamwise_{f}.csv", unpack=True)
        ax.loglog(k, scale * e, ls, color=c, lw=2.4 if f == "optimized_vector" else 1.8, label=lab)
    ax.set_xlabel("streamwise wavenumber $k_1$")
    ax.set_ylabel("energy spectrum")
    ax.set_ylim(5e-10, 1e-4)
    ax.legend(loc="lower left", fontsize=10)
    ax.set_title("Optimized vector filter recovers the DNS spectrum\non a 60 x 80 x 60 grid", loc="left")
    fig.text(0.01, -0.035, "Turbulent channel, M = 1.5, Re = 3000, degree-7 B-splines, |x2| = 0.94.",
             fontsize=8.5, color=GRAY)
    fig.savefig(ASSETS / "card_spectra.png")
    plt.close(fig)


def main() -> None:
    """Write every README asset."""
    plt.style.use(str(STYLE))
    hero()
    card_robustness()
    card_mms()
    card_spectra()
    card_accuracy_cost()


if __name__ == "__main__":
    main()
