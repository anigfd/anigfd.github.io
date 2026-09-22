"""Rossby-wave WKB ray tracing on a sphere (pure NumPy).

Stationary barotropic Rossby-wave rays on a background solid-body zonal
flow U(phi) = Omega_s*a*cos(phi) -- a uniform additional rotation
superimposed on the planetary rotation Omega. This is the classic
Hoskins & Karoly (1981) setup: adding solid-body rotation is equivalent to
simply rescaling the effective planetary rotation rate to Omega+Omega_s
(both U and the meridional gradient of absolute vorticity scale the same
way), and the resulting ray paths are exact great circles -- verified
numerically (see tests), not assumed.

Ray equations use the zonal wavenumber index n and meridional index m,
canonically conjugate to longitude lambda and latitude phi (n is exactly
conserved by the zonal symmetry of U and beta_eff, so it is a fixed
parameter here rather than a fourth evolved state). Physical wavenumbers
are kx=n/(a*cos(phi)), ky=m/a.
"""
import numpy as np


def stationary_wavenumber2(Omega_s, Omega=1.0, a=1.0):
    """K_s^2 = beta_eff/U for a STATIONARY (omega=0) ray. Independent of
    latitude for solid-body U -- that latitude-independence is exactly what
    makes the ray paths great circles."""
    return 2.0 * (Omega + Omega_s) / (Omega_s * a ** 2)


def dispersion_omega(phi, n, m, Omega_s, Omega=1.0, a=1.0):
    """Local barotropic Rossby-wave dispersion relation on the background
    flow U(phi)=Omega_s*a*cos(phi):

        omega = U*kx - beta_eff*kx/(kx^2+ky^2)

    with beta_eff(phi) = 2*(Omega+Omega_s)*cos(phi)/a the meridional
    gradient of absolute vorticity of the basic state.
    """
    kx = n / (a * np.cos(phi))
    ky = m / a
    U = Omega_s * a * np.cos(phi)
    beta_eff = 2.0 * (Omega + Omega_s) * np.cos(phi) / a
    K2 = kx ** 2 + ky ** 2
    return U * kx - beta_eff * kx / K2


def ray_rhs(state, n, Omega_s, Omega=1.0, a=1.0, h=1e-6):
    """d(lambda,phi,m)/dt via Hamilton's ray equations, using central
    differences of dispersion_omega (safer than a hand-derived symbolic
    ray equation -- correctness rests only on dispersion_omega itself,
    which is independently checked against the great-circle result).
    """
    lam, phi, m = state
    domega_dm = (dispersion_omega(phi, n, m + h, Omega_s, Omega, a)
                 - dispersion_omega(phi, n, m - h, Omega_s, Omega, a)) / (2 * h)
    domega_dphi = (dispersion_omega(phi + h, n, m, Omega_s, Omega, a)
                   - dispersion_omega(phi - h, n, m, Omega_s, Omega, a)) / (2 * h)
    domega_dn = (dispersion_omega(phi, n + h, m, Omega_s, Omega, a)
                 - dispersion_omega(phi, n - h, m, Omega_s, Omega, a)) / (2 * h)
    return np.array([domega_dn, domega_dm, -domega_dphi])


def launch_state(lam0, phi0, alpha0, Omega_s, Omega=1.0, a=1.0):
    """Initial (lambda, phi, m) and the conserved n for a stationary ray
    launched from (lam0,phi0) at angle alpha0 (radians, 0=east, pi/2=north)."""
    Ks = np.sqrt(stationary_wavenumber2(Omega_s, Omega, a))
    kx0, ky0 = Ks * np.cos(alpha0), Ks * np.sin(alpha0)
    n = kx0 * a * np.cos(phi0)
    m0 = ky0 * a
    return np.array([lam0, phi0, m0]), n


def great_circle_deviation(lam, phi):
    """Max |r.n_hat| for the ray's unit-sphere position r(t), where n_hat is
    the plane normal fit from the first and midpoint samples -- ~0 for an
    exact great circle, and a direct, coordinate-free correctness check.
    """
    r = np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)], axis=-1)
    n_hat = np.cross(r[0], r[len(r) // 4])
    n_hat = n_hat / np.linalg.norm(n_hat)
    return np.abs(r @ n_hat)
