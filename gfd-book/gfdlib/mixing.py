"""Eddy transport & mixing: a passive tracer stirred by 2D turbulence, and
the effective-diffusivity diagnostic it defines (pure NumPy).

A tracer's total concentration is split into an imposed, unbounded MEAN
gradient plus a periodic perturbation: c_total = Gamma*y + c'(x,y,t). c'
obeys the same advection-diffusion physics as any scalar, with the mean
gradient appearing as a source term -Gamma*v (v the meridional velocity
from the SAME vorticity field that is independently evolving under its own
barotropic dynamics -- rhs_coupled stacks both fields and reuses
gfdlib.spectral.Grid exactly as ch18's vorticity solver does):

    c'_t + J(psi,c') = -Gamma*v + kappa*lap(c')     (kappa folded into the
                                                       caller's linear op L,
                                                       exactly as ch18 does
                                                       for hyperviscosity)

The eddy (turbulent) diffusivity is defined by the down-gradient flux-
gradient relation K_eff = -mean(v'c')/Gamma. Because the c' equation is
LINEAR in c', K_eff must be exactly independent of Gamma for a fixed flow
realization -- verified numerically (not assumed) to be Gamma-independent
to machine precision, and to vanish exactly when v=0 identically (no
possible eddy flux without a velocity field). See gfdlib tests.
"""
import numpy as np


def rhs_coupled(state_hat, grid, Gamma):
    """Nonlinear tendency for state_hat=(zeta_hat, c_hat): vorticity
    self-advects under its own streamfunction (ch18's -J(psi,zeta), no
    tracer feedback -- the tracer is passive), and the tracer perturbation
    is advected by that same streamfunction plus the mean-gradient source
    -Gamma*v. Dissipation (hyperviscosity on zeta, diffusivity on c) is left
    to the caller's diagonal linear operator L, exactly as in ch18/ch16.
    """
    zeta_hat, c_hat = state_hat
    psi_hat = grid.invert_laplacian(zeta_hat)
    v_hat = grid.ddx(psi_hat)
    dzeta_hat = -grid.jacobian(psi_hat, zeta_hat)
    dc_hat = -grid.jacobian(psi_hat, c_hat) - Gamma * v_hat
    return np.stack([dzeta_hat, dc_hat])


def effective_diffusivity(v, c_prime, Gamma):
    """K_eff = -mean(v*c')/Gamma: the eddy diffusivity implied by the
    domain-averaged down-gradient tracer flux, for physical-space fields v,
    c_prime on the same doubly-periodic grid (so the mean is a true,
    homogeneous spatial average)."""
    return -np.mean(v * c_prime) / Gamma
