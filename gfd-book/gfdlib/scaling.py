"""Dimensionless numbers and reduced-model regime classification (pure
NumPy). See NOTATION.md's "Dimensionless numbers" table for definitions;
this module turns those definitions into small, testable functions plus a
regime classifier for the Ro-Bu plane used in chs. 2 and 4.
"""
import numpy as np


def rossby_number(U, f, L):
    """Ro = U/(f*L): advective time scale / rotation time scale."""
    return U / (f * L)


def ekman_number(nu, f, H):
    """Ek = nu/(f*H^2): viscous diffusion time scale / rotation time scale."""
    return nu / (f * H ** 2)


def froude_number(U, N, H):
    """Fr = U/(N*H): advective time scale / buoyancy-oscillation time scale."""
    return U / (N * H)


def reynolds_number(U, L, nu):
    """Re = U*L/nu: inertial force / viscous force."""
    return U * L / nu


def deformation_radius(N, H, f):
    """L_R = N*H/f."""
    return N * H / f


def burger_number(N, H, f, L):
    """Bu = (N*H/(f*L))^2 = (L_R/L)^2: (deformation radius / length scale)^2."""
    return (deformation_radius(N, H, f) / L) ** 2


def coriolis_parameter(latitude_deg, Omega=7.292e-5):
    """f = 2*Omega*sin(latitude); Omega defaults to Earth's rotation rate
    (rad/s). f=0 exactly at the equator -- callers must guard against that
    before dividing by f (Ro, Ek, Bu all blow up there, correctly: the
    equatorial beta-plane is a genuinely different regime, outside this
    module's scope).
    """
    return 2.0 * Omega * np.sin(np.radians(latitude_deg))


def classify_regime(Ro, Bu):
    """Coarse reduced-model regime label from (Ro, Bu), for the ch04
    regime map:

      0 "unbalanced"      Ro > 1        -- inertia-gravity waves dominate;
                                            no steady balance to exploit.
      1 "QG, baroclinic"  Ro<1, Bu>1    -- L<L_R: stratification controls
                                            the response (chs. 8, 16).
      2 "QG, barotropic"  Ro<1, Bu<1    -- L>L_R: rotation controls the
                                            response, flow tends columnar
                                            (Taylor-Proudman, ch. 3; ch. 18).

    A pedagogical simplification of the classical Ro-Bu regime diagram
    (Vallis, AOFD, Fig. 5.1-style): real transitions are smooth, not the
    sharp lines drawn here, and Bu=1 (L~L_R) is where baroclinic
    instability (ch. 16) is most efficient rather than a hard boundary.
    """
    Ro = np.asarray(Ro, dtype=float)
    Bu = np.asarray(Bu, dtype=float)
    return np.where(Ro > 1.0, 0, np.where(Bu > 1.0, 1, 2))
