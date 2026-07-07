"""Consistent figure policy for every notebook (matplotlib).

Colormap policy:
  * signed fields (vorticity, PV anomaly, height anomaly) -> diverging, symmetric
  * non-negative fields (speed, energy, temperature)      -> perceptually uniform
Always create fig,ax explicitly (never bare plt.*) so exports are deterministic.
"""
import numpy as np
import matplotlib.pyplot as plt

SEQ = "viridis"      # non-negative / magnitude fields
DIV = "RdBu_r"       # signed fields, use symmetric limits


def new_fig(w=5.0, h=4.0):
    fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
    return fig, ax


def field(ax, f, grid=None, signed=False, title="", cbar=True, **kw):
    """imshow a 2D field with the house colormap policy."""
    if signed:
        m = np.max(np.abs(f)) or 1.0
        kw.setdefault("cmap", DIV); kw.setdefault("vmin", -m); kw.setdefault("vmax", m)
    else:
        kw.setdefault("cmap", SEQ)
    ext = None if grid is None else [0, grid.L, 0, grid.L]
    im = ax.imshow(f.T, origin="lower", extent=ext, aspect="equal", **kw)
    ax.set_title(title)
    if cbar:
        ax.figure.colorbar(im, ax=ax, shrink=0.85)
    return im
