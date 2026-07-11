"""Inertial oscillations and the Taylor-Proudman limit (pure NumPy).

Unforced, frictionless horizontal momentum on an f-plane, du/dt=f*v,
dv/dt=-f*u, has the EXACT solution (u+iv)(t) = (u0+iv0)*exp(-i*f*t): a
clockwise (f>0) rotation at constant speed. Integrating once more gives a
circular trajectory of radius |u0+iv0|/f and period 2*pi/f -- the textbook
"inertial oscillation" seen throughout upper-ocean velocity records.

Taylor-Proudman is not a separate solver here: it is exactly the b=0
special case of gfdlib.balance's thermal-wind relation (f*u_z=-b_y) -- ch03
reuses gfdlib.balance directly for that part rather than re-deriving it.
"""
import numpy as np


def inertial_rhs(state, f):
    """RHS of du/dt=f*v, dv/dt=-f*u for state=(u,v) -- pass to
    gfdlib.timestep.rk4 for a numerical (non-analytic) integration check.
    """
    u, v = state
    return np.array([f * v, -f * u])


def inertial_trajectory(t, u0, v0, f, x0=0.0, y0=0.0):
    """Exact analytic solution. Writing w=u+iv, dw/dt=-i*f*w gives
    w(t)=w0*exp(-i*f*t); integrating once more (w0 constant) gives the
    position. Returns (x(t), y(t), u(t), v(t)), each broadcasting over t.
    """
    w0 = u0 + 1j * v0
    w = w0 * np.exp(-1j * f * t)
    z0 = x0 + 1j * y0
    z = z0 + (w0 / (-1j * f)) * (np.exp(-1j * f * t) - 1.0)
    return z.real, z.imag, w.real, w.imag
