"""Symmetric and inertial instability of a geostrophically-balanced front:
the Ertel-PV sign criterion, derived here directly (pure NumPy).

Consider a zonal thermal-wind front U=U(y,z), f ∂U/∂z = -∂b/∂y, with
buoyancy b=b(y,z), N^2=∂b/∂z. Frictionless, adiabatic motion in the (y,z)
plane conserves absolute momentum M=U-fy and buoyancy b along parcel
trajectories, so the front is stable to slantwise (symmetric) displacements
iff the Ertel PV q=(absolute vorticity)*grad(b) has the same sign as f:

    q = (f - dU/dy) N^2 - f (dU/dz)^2
      = (f - dU/dy) N^2 - f (db/dy)^2/f^2        [thermal wind: dU/dz=-(db/dy)/f]

Writing Ro=(dU/dy)/f (cross-front relative-vorticity Rossby number, note
sign) and Ri=N^2/(dU/dz)^2 (gradient Richardson number based on the
thermal-wind shear), q/(f N^2) = (1-Ro) - 1/Ri. This module uses the sign
convention Ro=-(dU/dy)/f (so f(1+Ro) is the absolute vertical vorticity,
matching the usual frontal-dynamics literature, e.g. Thomas/Tandon/
Mahadevan 2013); instability is q*f<0, i.e. (1+Ro) < 1/Ri.

Limits, each independently checkable:
  * N^2 -> large (stratification dominates): (1+Ro)<0 -- pure INERTIAL
    instability (absolute vertical vorticity changes sign), independent of
    Ri.
  * Ro=0 (no horizontal shear): Ri<1 -- the classic symmetric-instability
    threshold for a front with zero relative vorticity.
  * Ri<0 (N^2<0): unconditionally unstable -- ordinary GRAVITATIONAL
    (convective) instability, regardless of Ro.
"""
import numpy as np

STABLE, SYMMETRIC, INERTIAL, GRAVITATIONAL = 0, 1, 2, 3


def ertel_pv_ratio(Ro, Ri):
    """q/(f*N^2) = (1+Ro) - 1/Ri for Ri>0 (N^2>0); undefined sign meaning
    for Ri<=0, handled separately by classify (gravitational instability)."""
    return (1.0 + Ro) - 1.0 / Ri


def critical_ri(Ro):
    """The symmetric-instability marginal curve Ri_c=1/(1+Ro), valid for
    Ro>-1 (for Ro<=-1 the front is inertially unstable at every Ri>0)."""
    return 1.0 / (1.0 + Ro)


def classify(Ro, Ri):
    """Stability regime at a point (Ro,Ri) as one of STABLE, SYMMETRIC,
    INERTIAL, GRAVITATIONAL. Vectorized: Ro, Ri may be arrays (e.g. a
    meshgrid), and the return value is an integer array of the same shape,
    suitable for a discrete-colormap map plot."""
    Ro, Ri = np.broadcast_arrays(np.asarray(Ro, dtype=float), np.asarray(Ri, dtype=float))
    out = np.full(Ro.shape, STABLE, dtype=int)
    out = np.where(Ri <= 0, GRAVITATIONAL, out)
    out = np.where((Ri > 0) & (1.0 + Ro <= 0), INERTIAL, out)
    unstable_balanced = (Ri > 0) & (1.0 + Ro > 0) & ((1.0 + Ro) < 1.0 / Ri)
    out = np.where(unstable_balanced, SYMMETRIC, out)
    return out
