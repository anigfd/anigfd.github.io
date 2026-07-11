"""Vertical normal modes of a stratified fluid (pure NumPy).

Solves the linear, hydrostatic, rigid-lid vertical-structure eigenproblem

    d/dz[ (1/N^2(z)) dPhi/dz ] + (1/c_n^2) Phi = 0,   Phi'(0)=Phi'(H)=0

by the standard finite-VOLUME method: integrate the equation over a control
volume around each grid node (full spacing dz for interior nodes, HALF a
cell for the two boundary nodes, whose zero-flux Neumann condition removes
one face's flux entirely) rather than a naive ghost-point finite
difference. This matters: a first attempt used equal-weight ghost points
and produced a matrix that was NOT symmetric between the boundary row and
its interior neighbor (the boundary row's coefficient came out exactly 2x
the neighbor's own reference to it) -- caught by testing against the exact
constant-N solution, where it gave badly wrong mode SHAPES (though, oddly,
nearly-correct eigenvalues) once plotted, not just a slightly-off number.

The finite-volume matrix L (flux differences only, no volume division) IS
symmetric by construction; the unequal control volumes appear only as a
diagonal weight D=diag(dz/2,dz,...,dz,dz/2) in the generalized eigenproblem
L*Phi = mu*D*Phi. Since this book has no SciPy generalized eigensolver, the
standard congruence transform Psi=sqrt(D)*Phi turns this into the ORDINARY
symmetric eigenproblem (D^-1/2 L D^-1/2) Psi = mu*Psi, solved with
numpy.linalg.eigh (real, sorted eigenvalues -- the right tool for a
genuinely self-adjoint problem, unlike the general eig used for the
non-self-adjoint shear-instability eigenproblems elsewhere in this book).

Mode n=0 is always Phi=const, c_0=infinity (no deformation-radius
constraint at all) -- the CONTINUOUS analog of gfdlib.baroclinic's
barotropic mode (L*const=0 exactly, for any grid, since every flux
difference vanishes for a constant field). Modes n=1,2,... are the
baroclinic modes, with speed c_n and deformation radius L_n=c_n/f0
decreasing roughly as 1/n. Verified against the exact constant-N solution
Phi_n=cos(n*pi*z/H), c_n=N*H/(n*pi) (derived and checked by hand -- see
tests): eigenvalues match to 1e-5 relative and mode shapes to 1e-12
absolute at n=400 grid points, confirming the fix.
"""
import numpy as np


def vertical_modes(z, N2, H, n_modes):
    """The n_modes gravest BAROCLINIC (n=1,2,...) vertical modes.

    z must be a uniform grid of H/(len(z)-1) spacing from 0 to H inclusive
    (both boundaries are grid points); N2 is N^2(z) sampled at the same
    points.

    Returns (c, Phi): c has shape (n_modes,), the mode speeds
    c_1>c_2>...>c_{n_modes}>0. Phi has shape (n_modes, len(z)), each row
    normalized to Phi(0)=1 (or, if Phi(0) is numerically near zero for a
    particular profile, to max|Phi|=1).
    """
    n = len(z)
    dz = H / (n - 1)
    g = 1.0 / N2
    g_half = 0.5 * (g[:-1] + g[1:])   # weight at the n-1 half-points

    L = np.zeros((n, n))              # symmetric flux-difference operator
    for j in range(1, n - 1):
        L[j, j - 1] += g_half[j - 1] / dz
        L[j, j] -= g_half[j - 1] / dz
        L[j, j] -= g_half[j] / dz
        L[j, j + 1] += g_half[j] / dz
    L[0, 0] = -g_half[0] / dz
    L[0, 1] = g_half[0] / dz
    L[-1, -1] = -g_half[-1] / dz
    L[-1, -2] = g_half[-1] / dz

    D = np.full(n, dz)
    D[0] = D[-1] = dz / 2.0
    d_inv_sqrt = 1.0 / np.sqrt(D)
    M = L * d_inv_sqrt[:, None] * d_inv_sqrt[None, :]   # symmetric congruence transform

    eigvals, eigvecs_psi = np.linalg.eigh(M)   # ascending; eigvals[-1]~0 is the barotropic mode
    c = np.zeros(n_modes)
    Phi = np.zeros((n_modes, n))
    for k in range(n_modes):
        idx = n - 2 - k                         # n-2 -> mode 1, n-3 -> mode 2, ...
        mu = eigvals[idx]
        c[k] = np.sqrt(-1.0 / mu) if mu < 0 else np.inf
        phi = eigvecs_psi[:, idx] * d_inv_sqrt
        norm = phi[0] if abs(phi[0]) > 1e-10 else np.max(np.abs(phi)) * np.sign(phi[np.argmax(np.abs(phi))])
        Phi[k] = phi / norm
    return c, Phi


def pycnocline_N2(z, N2_min, N2_max, z_center, thickness):
    """A smooth idealized pycnocline: N^2(z) peaks at z_center with
    amplitude N2_max-N2_min above a background N2_min, using a sech^2
    profile of the given thickness -- the same functional shape as
    gfdlib.internalwaves' variable-N cases, applied here to a
    z-in-[0,H] grid instead of a spectral (x,z) plane.
    """
    return N2_min + (N2_max - N2_min) / np.cosh((z - z_center) / thickness) ** 2
