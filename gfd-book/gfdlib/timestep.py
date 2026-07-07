"""Time integrators (pure NumPy). Operate on flat state arrays or spectral fields."""
import numpy as np


def rk4(rhs, y, dt, t=0.0):
    """One classical RK4 step of dy/dt = rhs(t, y)."""
    k1 = rhs(t, y)
    k2 = rhs(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = rhs(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = rhs(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def etdrk_diffusion_step(F, rhs_nl, dt, k2, nu, t=0.0):
    """IMEX: exact linear diffusion factor + explicit RK2 for the nonlinear part.

    Solves dF/dt = -nu*k2*F + rhs_nl(t,F) in spectral space. The diffusion is
    handled with an integrating factor (unconditionally stable for that term),
    which lets the notebook use a comfortable dt at 128^2-class resolution.
    """
    E = np.exp(-nu * k2 * dt)
    N1 = rhs_nl(t, F)
    F1 = E * (F + dt * N1)
    N2 = rhs_nl(t + dt, F1)
    return E * F + dt * 0.5 * (E * N1 + N2)
