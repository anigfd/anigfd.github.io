"""Geostrophic, hydrostatic, and thermal-wind balance (pure NumPy).

A meridional (y,z) cross-section through an idealized buoyancy front. Every
quantity here is a closed-form expression -- no discretization, no time
stepping -- so gradients and integrals are evaluated exactly rather than
approximated. See NOTATION.md for the balance sign conventions.
"""
import numpy as np


def buoyancy_field(y, z, db, Ly, Htrop, N2=0.0):
    """b(y,z) = N2*z + db*tanh(y/Ly)*cos(pi*z/(2*Htrop)).

    The background term N2*z carries no meridional gradient (and therefore
    drives no thermal-wind shear); it only tilts the isopycnals for realism.
    The frontal term changes sign at z=Htrop -- an idealized tropopause.
    """
    return N2 * z + db * np.tanh(y / Ly) * np.cos(np.pi * z / (2 * Htrop))


def dbdy_frontal(y, z, db, Ly, Htrop):
    """Exact d(b)/dy = (db/Ly) * sech^2(y/Ly) * cos(pi*z/(2*Htrop))."""
    return (db / Ly) / np.cosh(y / Ly) ** 2 * np.cos(np.pi * z / (2 * Htrop))


def thermal_wind_u(y, z, db, Ly, Htrop, f, u0=0.0):
    """u_g(y,z) = u0 - (1/f) * integral_0^z dbdy_frontal(y,z') dz'.

    Closed form: u0 - (db/(f*Ly)) * sech^2(y/Ly) * (2*Htrop/pi) * sin(pi*z/(2*Htrop)).
    The integrand changes sign at z=Htrop, so u_g has an extremum there at
    every y -- the geostrophic jet core sits exactly where the meridional
    buoyancy gradient reverses.
    """
    sech2 = 1.0 / np.cosh(y / Ly) ** 2
    integral_g = (2 * Htrop / np.pi) * np.sin(np.pi * z / (2 * Htrop))
    return u0 - (db / (f * Ly)) * sech2 * integral_g
