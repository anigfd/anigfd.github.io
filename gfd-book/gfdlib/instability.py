"""Barotropic instability: linear growth rates and nonlinear roll-up ICs
(pure NumPy).

The linear part is a genuinely new numerical primitive for this book: a
matrix eigenvalue solve (not spectral time-stepping). The nonlinear part
reuses gfdlib.spectral.Grid and gfdlib.timestep.ifrk4_step exactly as in
ch18/ch07 -- only the initial condition (a periodic double shear layer) is
new.
"""
import numpy as np


def _laplacian_1d(n, dy):
    """Second-derivative matrix, Dirichlet BCs, interior points only."""
    return (np.diag(-2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1)) / dy ** 2


def growth_rate(y, U, k, beta=0.0):
    """Fastest-growing normal-mode growth rate k*c_i for the barotropic
    Rayleigh-Kuo equation

        (U-c)(phi''-k^2 phi) + (beta-U'')phi = 0

    on a channel with rigid walls (phi=0) just outside y[0],y[-1]. U'' is
    computed with the SAME discrete second-derivative operator used for
    phi, for exact numerical consistency between the base-state curvature
    and the perturbation operator.

    Rewritten as a generalized eigenvalue problem A phi = c B phi with
    A = U*B + diag(beta-U''), B = D^2-k^2*I; converted to the standard form
    B^-1 A (B is invertible for k>0) and solved with plain
    numpy.linalg.eigvals -- no SciPy, no generalized eigensolver needed.
    """
    n = len(y)
    dy = y[1] - y[0]
    D2 = _laplacian_1d(n, dy)
    Upp = D2 @ U
    B = D2 - k ** 2 * np.eye(n)
    A = U[:, None] * B + np.diag(beta - Upp)
    c = np.linalg.eigvals(np.linalg.solve(B, A))
    return k * max(float(np.max(c.imag)), 0.0)


def growth_rate_curve(y, U, k_values, beta=0.0):
    """growth_rate evaluated at every k in k_values."""
    return np.array([growth_rate(y, U, k, beta) for k in k_values])


def double_shear_layer(grid, delta, v_pert=0.05, n_pert=1):
    """Vorticity field for a periodic double shear layer: two opposite-signed
    tanh jets sharing a doubly-periodic domain (so the profile closes up),
    the standard barotropic-instability roll-up test case. A small
    sinusoidal v-perturbation at zonal wavenumber n_pert seeds the
    instability so it rolls up within a practical run time.
    """
    y1, y2 = 0.25 * grid.L, 0.75 * grid.L
    zeta = -(1.0 / delta) * (
        1.0 / np.cosh((grid.y - y1) / delta) ** 2
        - 1.0 / np.cosh((grid.y - y2) / delta) ** 2
    )
    zeta += v_pert * (2 * np.pi * n_pert / grid.L) * np.cos(2 * np.pi * n_pert * grid.x / grid.L)
    return zeta
