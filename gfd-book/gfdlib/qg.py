"""Single-layer QG on a beta-plane (pure NumPy).

Deliberately thin: the QGPV Helmholtz operator (nabla^2 - 1/L_R^2) is
mathematically identical to the linear shallow-water PV inversion in
gfdlib.shallowwater once nondimensionalized by L_R (L_R_hat=1) -- so
gfdlib.shallowwater.invert_pv IS the QG psi-from-q inversion, no new
primitive needed. gfdlib.pv.gaussian_blob places the initial vortex, and
gfdlib.timestep.ifrk4_step integrates it exactly as in ch18, with a
deformation-radius-modified linear operator (see dispersion_omega below).
"""
import numpy as np


def dispersion_omega(kx, ky, beta):
    """QG Rossby-wave dispersion relation (L_R_hat=1):

        omega = -beta*kx/(kx^2+ky^2+1)

    The "+1" (absent from the barotropic relation in ch18) is the
    deformation-radius screening term. Unlike the barotropic case, where
    |omega| grows without bound as k->0, this is bounded: |omega|<=beta/2,
    attained at |k|=1 (i.e. a wavelength comparable to L_R).
    """
    return -beta * kx / (kx ** 2 + ky ** 2 + 1.0)
