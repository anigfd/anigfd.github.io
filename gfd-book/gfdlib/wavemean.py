"""Wave-mean-flow interaction, first pass: Stokes drift and Rossby-wave
critical layers (pure NumPy).

Two independent pieces, both first-order-in-amplitude effects a wave has
on the mean state it rides on:

  * Stokes drift -- even though the EULERIAN mean of a wave's velocity is
    zero by definition, the LAGRANGIAN mean (the actual net drift of a
    fluid parcel) is not, at second order in wave amplitude. Derived here
    from scratch for a single 1D Fourier wave, not looked up.
  * Rossby-wave critical layers -- a wave riding on a sheared mean flow
    U(y) develops a singular structure (local wavenumber and, by wave-
    action conservation, amplitude, both diverging) as it approaches the
    latitude where the Doppler-shifted phase speed matches U(y). This is
    the linear, inviscid ray-theory signature of the process that -- once
    any dissipation or nonlinearity is added -- deposits the wave's
    momentum into the mean flow there (the physical content of "wave-mean
    interaction", explored further in ch. 23's wave-activity framework).
"""
import numpy as np


def wave_velocity(x, t, A, k, omega):
    """u'(x,t) = A*cos(k*x - omega*t): a single Fourier wave, Eulerian-mean
    velocity identically zero at every fixed x."""
    return A * np.cos(k * x - omega * t)


def stokes_drift(A, k, omega):
    """Lagrangian-mean drift, O(A^2), for wave_velocity above.

    Derivation: a parcel's displacement from its mean position x0, to
    leading order, is xi(t)=integral of wave_velocity(x0,t)dt =
    (A/omega)*sin(omega*t-k*x0). The Lagrangian mean velocity is the time
    average of wave_velocity(x0+xi,t) Taylor-expanded to O(xi):

        u_L = <wave_velocity(x0,t)> + <xi * d(wave_velocity)/dx|_{x0,t}>

    The first term is the Eulerian mean, zero by construction. The second
    (the Stokes drift) works out to A^2*k/(2*omega) after using
    <sin^2>=1/2 -- positive, i.e. in the direction of phase propagation,
    consistent with the classical surface-gravity-wave result.
    """
    return A ** 2 * k / (2.0 * omega)


def rossby_shear_dispersion(y, k, l, U_func, beta):
    """omega(y,k,l) = U(y)*k - beta*k/(k^2+l^2): the Doppler-shifted
    barotropic Rossby dispersion relation on a background zonal shear
    U(y) (reduces to gfdlib.qg's dispersion at U=0)."""
    return U_func(y) * k - beta * k / (k ** 2 + l ** 2)


def ray_rhs_shear(state, k, beta, dUdy_func):
    """Hamilton's ray equations for state=(y,l), at fixed (conserved) k,
    from rossby_shear_dispersion:

        dy/dt = d(omega)/dl = 2*beta*k*l/(k^2+l^2)^2
        dl/dt = -d(omega)/dy = -k*dU/dy

    (the beta term carries no EXPLICIT y-dependence at fixed k,l, since l
    is the independent ray variable -- only the mean-flow term contributes
    to -d(omega)/dy). As the ray approaches a critical layer where
    U(y)=omega/k, the beta term must shrink to keep omega fixed, forcing
    k^2+l^2 -> infinity: l diverges and dy/dt -> 0, the ray asymptotically
    approaching but never (in finite time, in this inviscid linear theory)
    reaching the critical layer -- verified directly by integration, not
    just asserted (see tests and the ch13 notebook).
    """
    y, l = state
    K2 = k ** 2 + l ** 2
    dydt = 2.0 * beta * k * l / K2 ** 2
    dldt = -k * dUdy_func(y)
    return np.array([dydt, dldt])
