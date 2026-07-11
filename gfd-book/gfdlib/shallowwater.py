"""Linearized rotating shallow water on a beta-plane (pure NumPy).

Nondimensionalized by the deformation radius L_R = c/f0 (length), f0^-1
(time), and the mean depth H (height, so eta_hat = eta/H and c_hat = 1).
Symbols and sign conventions match NOTATION.md.
"""
import numpy as np


def coriolis(grid, beta, y0=None):
    """f_hat(y) = 1 + beta*(y - y0); y0 centers the beta term on the domain."""
    if y0 is None:
        y0 = 0.5 * grid.L
    return 1.0 + beta * (grid.y - y0)


def rhs(state, grid, f):
    """RHS of d(state)/dt for state=(eta,u,v) stacked as one (3,n,n) array.

        eta_t + u_x + v_y = 0
        u_t - f v = -eta_x
        v_t + f u = -eta_y

    f can be a scalar (f-plane) or grid.y-shaped array (beta-plane, see
    `coriolis`). Linear and exact in spectral space except for the f*u, f*v
    terms, which are evaluated in physical space since f=f(y).
    """
    eta, u, v = state
    eta_hat, u_hat, v_hat = grid.fft(eta), grid.fft(u), grid.fft(v)
    eta_x = grid.ifft(grid.ddx(eta_hat))
    eta_y = grid.ifft(grid.ddy(eta_hat))
    div = grid.ifft(grid.ddx(u_hat) + grid.ddy(v_hat))
    return np.stack([-div, f * v - eta_x, -f * u - eta_y])


def potential_vorticity(state, grid):
    """Linear PV anomaly q = zeta - eta (zeta = v_x - u_y).

    Obeys dq/dt = -beta*v (see NOTATION.md): exactly conserved pointwise on
    the f-plane (beta=0). Geostrophic adjustment redistributes eta and u,v
    but cannot change q — the balanced end state is fixed by the initial q.
    """
    eta, u, v = state
    u_hat, v_hat = grid.fft(u), grid.fft(v)
    zeta = grid.ifft(grid.ddx(v_hat) - grid.ddy(u_hat))
    return zeta - eta


def invert_pv(q, grid):
    """Geostrophic height field carrying PV anomaly q: solve

        (nabla^2 - 1) eta_bal = q   <=>   eta_bal_hat = -q_hat/(k^2+1)

    i.e. the balance in which zeta_g = nabla^2 eta_bal (geostrophic vorticity)
    and q = zeta_g - eta_bal simultaneously hold.
    """
    q_hat = grid.fft(q)
    return grid.ifft(-q_hat / (grid.k2 + 1.0))


def divergence(state, grid):
    """delta = u_x + v_y, computed spectrally (used by ch. 24).

    The single sharpest slow/fast discriminator in rotating shallow water:
    the geostrophic (slow) mode has delta = 0 identically, while
    inertia-gravity waves are MADE of divergence -- so the RMS of delta is
    a direct meter of how far off the slow manifold a state is. A
    geostrophically balanced state (u = -eta_y, v = eta_x) returns zero to
    spectral accuracy; see tests.
    """
    eta, u, v = state
    return grid.ifft(grid.ddx(grid.fft(u)) + grid.ddy(grid.fft(v)))
