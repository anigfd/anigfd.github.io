"""Time integrators (pure NumPy). Operate on flat state arrays or spectral fields."""
import numpy as np


def rk4(rhs, y, dt, t=0.0):
    """One classical RK4 step of dy/dt = rhs(t, y)."""
    k1 = rhs(t, y)
    k2 = rhs(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = rhs(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = rhs(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def ifrk4_step(F, rhs_nl, dt, L, t=0.0):
    """Integrating-factor RK4: exact diagonal linear propagator + RK4 nonlinear.

    Solves dF/dt = L*F + rhs_nl(t,F) in spectral space, where L is any diagonal
    (possibly complex) linear operator — e.g. L = -nu*k2**n_nu + 1j*beta*kx/k2
    for hyperviscous dissipation plus exact Rossby-wave propagation. The linear
    term is unconditionally stable (integrating factor); the nonlinear term gets
    classical RK4, whose stability region covers the imaginary axis, so pure
    advection at nu=0 is stable at CFL-limited dt (an RK2 IMEX is not: it slowly
    amplifies oscillatory modes and blows up over thousands of steps).
    """
    E = np.exp(L * (0.5 * dt))
    E2 = E * E
    k1 = rhs_nl(t, F)
    k2 = rhs_nl(t + 0.5 * dt, E * (F + 0.5 * dt * k1))
    k3 = rhs_nl(t + 0.5 * dt, E * F + 0.5 * dt * k2)
    k4 = rhs_nl(t + dt, E2 * F + dt * E * k3)
    return E2 * F + (dt / 6.0) * (E2 * k1 + 2.0 * E * (k2 + k3) + k4)
