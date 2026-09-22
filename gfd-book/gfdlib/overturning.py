"""Buoyancy-driven overturning circulation: Munk's (1966) "abyssal
recipes" vertical balance, and a minimal 2D differential-heating-driven
overturning cell (Hadley-cell / meridional-overturning-circulation
mechanism), pure NumPy.

Both are derived here, not looked up. The 2D model reuses
gfdlib.convection.ChannelGrid's already-validated Boussinesq
vorticity-streamfunction machinery (poisson_solve, jacobian, laplacian)
unchanged, replacing Rayleigh-Benard's uniform bottom heating with a
differential, surface-concentrated heating pattern -- warm at one end of
the (periodic) x-domain, cool at the other -- which drives a SINGLE
overturning cell instead of convective rolls (the classic "horizontal
convection" mechanism, Rossby 1965): rising motion over the heated column,
sinking over the cooled one, closing via return flow at the surface and
at depth.
"""
import numpy as np


def abyssal_profile(z, w, kappa, H, T_bottom=0.0, T_top=1.0):
    """Munk's (1966) exact steady 1D balance w*dT/dz=kappa*d^2T/dz^2 on
    [0,H], T(0)=T_bottom, T(H)=T_top: an upwelling w (positive=upward)
    advecting the background stratification against diffusive mixing by
    kappa. General solution T(z)=A+B*exp(w*z/kappa); A,B fixed by the 2
    BCs via a small numpy.linalg.solve (not a memorized closed form).
    Large w/kappa sweeps the BOTTOM value through most of the interior,
    compressing the transition to the top value into a thin layer near
    z=H -- the deep ocean's observed near-uniform abyssal temperature,
    with the thermocline confined near the surface.
    """
    lam = w / kappa
    M = np.array([[1.0, 1.0], [1.0, np.exp(lam * H)]])
    A, B = np.linalg.solve(M, np.array([T_bottom, T_top]))
    z = np.asarray(z)
    return A + B * np.exp(lam * z)


def surface_heating(grid, Q0, z_power=3):
    """Idealized differential surface buoyancy forcing on a
    gfdlib.convection.ChannelGrid: Q(x,z)=Q0*cos(2*pi*x/Lx)*z^z_power --
    warm ("equator"/low-latitude) at x=0, cool ("pole"/high-latitude) at
    x=Lx/2, concentrated near the surface z=1 (z_power>1 suppresses deep
    heating) rather than uniform with depth, analogous to surface
    radiative/buoyancy forcing rather than a deep heat source.
    """
    xprof = np.cos(2 * np.pi * grid.x / grid.Lx)[:, None]
    zprof = grid.z[None, :] ** z_power
    return Q0 * xprof * zprof


def rhs_overturning(state, grid, Ra, Pr, Q):
    """RHS of d(state)/dt for state=(zeta,theta) stacked as one (2,nx,nz)
    array, on a gfdlib.convection.ChannelGrid:

        zeta_t  = -J(psi,zeta)  + Pr*lap(zeta) + Pr*Ra*theta_x
        theta_t = -J(psi,theta) +    lap(theta) + Q(x,z)
        nabla^2 psi = zeta

    Identical structure to convection.rhs_boussinesq, but with the fixed
    background-stratification advection term "+psi_x" replaced by an
    explicit, spatially-varying heating field Q (theta is no longer a
    perturbation about a background vertical gradient -- the circulation
    is driven entirely by Q, self-consistently coupled to the flow through
    the buoyancy term Pr*Ra*theta_x, exactly as in Rossby's (1965)
    horizontal-convection problem).
    """
    zeta, theta = state
    psi = grid.poisson_solve(zeta)
    dzeta = -grid.jacobian(psi, zeta) + Pr * grid.laplacian(zeta) + Pr * Ra * grid.ddx(theta)
    dtheta = -grid.jacobian(psi, theta) + grid.laplacian(theta) + Q
    return np.stack([dzeta, dtheta])
