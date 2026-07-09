"""Potential-vorticity field builders for inversion/advection demos (pure NumPy).

Pairs with `gfdlib.spectral.Grid` and `gfdlib.timestep.ifrk4_step` exactly as
in the barotropic-turbulence chapter -- these functions only construct
initial q(x,y) fields; inversion and evolution reuse existing primitives.
"""
import numpy as np


def gaussian_blob(grid, x0, y0, amp, sigma):
    """A single Gaussian PV anomaly at (x0,y0), periodic-safe.

    Uses the minimum-image (nearest periodic copy) distance so a blob placed
    near a domain edge doesn't get spuriously clipped by wraparound.
    """
    dx = (grid.x - x0 + grid.L / 2) % grid.L - grid.L / 2
    dy = (grid.y - y0 + grid.L / 2) % grid.L - grid.L / 2
    return amp * np.exp(-(dx ** 2 + dy ** 2) / (2 * sigma ** 2))


def staircase_pv(grid, n_steps, amp, sharpness=15.0):
    """An idealized PV staircase: n_steps alternating high/low bands in y,
    each separated by a sharp riser -- a smoothed square wave,
    tanh(sharpness*sin(2*pi*n_steps*y/L)), which is exact in the
    sharpness->infinity limit. Larger `sharpness` gives sharper risers.
    """
    return amp * np.tanh(sharpness * np.sin(2 * np.pi * n_steps * grid.y / grid.L))
