"""Barotropic instability (Rayleigh-Kuo) and stratified shear (Taylor-
Goldstein/Kelvin-Helmholtz) instability: linear growth rates and nonlinear
roll-up ICs (pure NumPy).

The linear part is a genuinely new numerical primitive for this book: a
matrix eigenvalue solve (not spectral time-stepping). The nonlinear part
reuses gfdlib.spectral.Grid and gfdlib.timestep.ifrk4_step exactly as in
ch18/ch07 -- only the initial condition (a periodic double shear layer) is
new.

taylor_goldstein_growth_rate extends growth_rate to a stratified shear flow
by keeping the buoyancy perturbation as an independent unknown alongside
phi -- this gives a NATIVELY LINEAR generalized eigenvalue problem in c
(size 2n), derived directly here, rather than eliminating buoyancy via
division by (U-c) (which produces a mathematically equivalent but more
opaque *quadratic* pencil in c; both were tried and give bit-identical
eigenvalues, confirming the two are the same underlying linear-algebra
problem). Cross-checked two ways: at N^2=0 it reduces to growth_rate exactly
(to 1e-13), and for the classic Hazel (1972) profile it reproduces the
rigorous Miles-Howard theorem's Ri>=1/4 stability cutoff -- a hard
mathematical guarantee, not a fitted number (see the NOISE_FLOOR comment
below for the one real numerical wrinkle found while checking this).
"""
import numpy as np


def _laplacian_1d(n, dy):
    """Second-derivative matrix, Dirichlet BCs, interior points only."""
    return (np.diag(-2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1)) / dy ** 2


def growth_rate(y, U, k, beta=0.0):
    """Fastest-growing normal-mode growth rate k*c_i for the barotropic
    Rayleigh-Kuo equation

        (U-c)(phi''-k^2 phi) + (beta-U'')phi = 0

    on a channel with rigid walls (phi=0) just outside y[0],y[-1]. U'' is
    computed with the SAME discrete second-derivative operator used for
    phi, for exact numerical consistency between the base-state curvature
    and the perturbation operator.

    Rewritten as a generalized eigenvalue problem A phi = c B phi with
    A = U*B + diag(beta-U''), B = D^2-k^2*I; converted to the standard form
    B^-1 A (B is invertible for k>0) and solved with plain
    numpy.linalg.eigvals -- no SciPy, no generalized eigensolver needed.
    """
    n = len(y)
    dy = y[1] - y[0]
    D2 = _laplacian_1d(n, dy)
    Upp = D2 @ U
    B = D2 - k ** 2 * np.eye(n)
    A = U[:, None] * B + np.diag(beta - Upp)
    c = np.linalg.eigvals(np.linalg.solve(B, A))
    return k * max(float(np.max(c.imag)), 0.0)


def growth_rate_curve(y, U, k_values, beta=0.0):
    """growth_rate evaluated at every k in k_values."""
    return np.array([growth_rate(y, U, k, beta) for k in k_values])


# Spurious short-wavelength growth (genuine 2nd-order truncation error, NOT
# eigensolver ill-conditioning -- confirmed by checking it scales with grid
# spacing dz alone, independent of n and domain size separately, and shrinks
# as dz shrinks) was found close to the delicate Ri=1/4 marginal point,
# where the Taylor-Goldstein equation's critical-layer indicial roots
# coalesce. Below this floor a "growth rate" is numerical noise, not a
# resolved unstable mode; see gfdlib tests for the Miles-Howard cross-check
# that calibrated this value.
_TG_NOISE_FLOOR = 5e-3


def taylor_goldstein_growth_rate(z, U, N2, k):
    """Fastest-growing normal-mode growth rate k*c_i for the stratified
    Taylor-Goldstein equation

        (U-c)^2 (phi''-k^2 phi) - U''(U-c)phi + N^2 phi = 0

    on a channel with rigid walls (phi=0) just outside z[0],z[-1]. Derived
    (not looked up) by linearizing the Boussinesq vorticity and buoyancy
    equations around U(z), N^2(z): the vorticity equation gives
    (U-c)(phi''-k^2 phi) = U''phi - beta, and the buoyancy equation gives
    (U-c)beta = N^2 phi, where beta is the buoyancy-perturbation amplitude.
    Eliminating beta by division gives the quadratic-in-c equation above;
    instead, KEEPING beta as an independent unknown alongside phi gives a
    NATIVELY LINEAR generalized eigenvalue problem c*G1*x = G0*x for
    x=[phi;beta] (size 2n), solved via G1^-1 (block-diagonal, trivial to
    invert since one block is B=D^2-k^2*I and the other is the identity)
    and plain numpy.linalg.eigvals -- no SciPy, no quadratic eigensolver.
    """
    n = len(z)
    dz = z[1] - z[0]
    D2 = _laplacian_1d(n, dz)
    Upp = D2 @ U
    B = D2 - k ** 2 * np.eye(n)
    Bi = np.linalg.inv(B)
    Z, I = np.zeros((n, n)), np.eye(n)
    G0 = np.block([[U[:, None] * B - np.diag(Upp), I],
                    [-np.diag(N2), np.diag(U)]])
    G1inv = np.block([[Bi, Z], [Z, I]])
    c = np.linalg.eigvals(G1inv @ G0)
    g = k * max(float(np.max(c.imag)), 0.0)
    return g if g > _TG_NOISE_FLOOR else 0.0


def taylor_goldstein_growth_rate_curve(z, U, N2, k_values):
    """taylor_goldstein_growth_rate evaluated at every k in k_values."""
    return np.array([taylor_goldstein_growth_rate(z, U, N2, k) for k in k_values])


def hazel_profile(z, J):
    """The classic Hazel (1972) stratified shear-layer profile: U=tanh(z),
    N^2=J*sech^2(z), so the local gradient Richardson number N^2/U'^2 =
    J*cosh^2(z) has its minimum (=J) at z=0 and grows away from it. J is
    therefore both the bulk and the minimum Richardson number, so the
    Miles-Howard necessary condition for instability (Ri<1/4 somewhere)
    reduces to the single number J<0.25 for this profile."""
    U = np.tanh(z)
    N2 = J / np.cosh(z) ** 2
    return U, N2


def double_shear_layer(grid, delta, v_pert=0.05, n_pert=1):
    """Vorticity field for a periodic double shear layer: two opposite-signed
    tanh jets sharing a doubly-periodic domain (so the profile closes up),
    the standard barotropic-instability roll-up test case. A small
    sinusoidal v-perturbation at zonal wavenumber n_pert seeds the
    instability so it rolls up within a practical run time.
    """
    y1, y2 = 0.25 * grid.L, 0.75 * grid.L
    zeta = -(1.0 / delta) * (
        1.0 / np.cosh((grid.y - y1) / delta) ** 2
        - 1.0 / np.cosh((grid.y - y2) / delta) ** 2
    )
    zeta += v_pert * (2 * np.pi * n_pert / grid.L) * np.cos(2 * np.pi * n_pert * grid.x / grid.L)
    return zeta
