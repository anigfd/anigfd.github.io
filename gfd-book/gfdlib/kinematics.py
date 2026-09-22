"""Kinematics of a prescribed 2D flow: material derivative, trajectories vs.
streamlines, and the Okubo-Weiss strain/vorticity decomposition (pure NumPy).

The flow is a classic textbook combination -- solid-body strain plus a
smooth (Lamb-Oseen-style) vortex, both individually divergence-free -- used
throughout the vortex-dynamics literature to illustrate when a vortex core
survives an ambient strain field and when it is torn into filaments. The
strain is STEADY, so a particle's trajectory differs from any instantaneous
streamline only through the vortex's own rotation.
"""
import numpy as np


def vortex_velocity(x, y, Gamma, sigma):
    """Tangential velocity of a Lamb-Oseen-like smooth vortex at the origin:
    v_theta(r) = (Gamma/(2*pi*r)) * (1 - exp(-r^2/(2*sigma^2))), regular at
    r=0 (Taylor-expands to solid-body rotation Gamma*r/(4*pi*sigma^2), so
    the core spins as a rigid disk and the far field decays to the ordinary
    Gamma/(2*pi*r) irrotational vortex). Returns Cartesian (u,v).
    """
    r2 = x ** 2 + y ** 2
    r = np.sqrt(np.where(r2 == 0, 1.0, r2))
    v_theta = (Gamma / (2 * np.pi * r)) * (1.0 - np.exp(-r2 / (2 * sigma ** 2)))
    v_theta = np.where(r2 == 0, 0.0, v_theta)
    return -v_theta * y / r, v_theta * x / r


def velocity_field(x, y, alpha, Gamma, sigma):
    """Steady strain (rate alpha, axes along x,y) plus a Lamb-Oseen vortex
    at the origin -- both individually divergence-free, so the sum is too:
    u = -alpha*x + u_vortex, v = alpha*y + v_vortex.
    """
    u_vort, v_vort = vortex_velocity(x, y, Gamma, sigma)
    return -alpha * x + u_vort, alpha * y + v_vort


def okubo_weiss(x, y, alpha, Gamma, sigma, h=1e-4):
    """W = Sn^2 + Ss^2 - zeta^2 (normal strain, shear strain, vorticity),
    from central differences of velocity_field -- safer than hand-deriving
    the vortex's messy analytic derivatives (the same central-difference
    pattern gfdlib.rossby uses for its ray equations, for the same reason:
    correctness then rests only on velocity_field itself).

    W<0: vorticity dominates (elliptic core, a coherent vortex survives);
    W>0: strain dominates (hyperbolic region, material lines stretch and
    fold into filaments).
    """
    u_xp, v_xp = velocity_field(x + h, y, alpha, Gamma, sigma)
    u_xm, v_xm = velocity_field(x - h, y, alpha, Gamma, sigma)
    u_yp, v_yp = velocity_field(x, y + h, alpha, Gamma, sigma)
    u_ym, v_ym = velocity_field(x, y - h, alpha, Gamma, sigma)
    u_x = (u_xp - u_xm) / (2 * h)
    v_x = (v_xp - v_xm) / (2 * h)
    u_y = (u_yp - u_ym) / (2 * h)
    v_y = (v_yp - v_ym) / (2 * h)
    Sn = u_x - v_y
    Ss = v_x + u_y
    zeta = v_x - u_y
    return Sn ** 2 + Ss ** 2 - zeta ** 2
