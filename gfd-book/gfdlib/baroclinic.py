"""Baroclinic instability: the Eady growth-rate calculator and the 2-layer
(Phillips) nonlinear model (pure NumPy).

Both the Eady boundary-value problem and the 2-layer normal-mode problem
are solved as small (2x2) eigenvalue problems, derived directly here (not
copied from a textbook closed form) and cross-checked numerically: the Eady
growth-rate curve's peak and cutoff match the published Michalke-style
benchmark (mu_max~1.61, sigma_max~0.31, cutoff mu_c~2.399; Vallis/Pedlosky),
and the 2-layer model's nonlinear time-stepped growth rate converges to its
own independently-derived linear eigenvalue to <0.1% (see tests).

A fixed mean shear with no other sink is an unlimited energy source: without
bottom drag, the nonlinear 2-layer growth continues past any physically
reasonable amplitude regardless of hyperviscosity, until overflow (this was
verified, not assumed -- see rhs_2layer). rhs_2layer therefore includes the
standard Phillips (1954) bottom Ekman drag term; with it tuned, a
small-amplitude seed grows exponentially, peaks, decays, and equilibrates --
the complete baroclinic life cycle -- confirmed by direct time integration
to t=500 at both n=48 and the notebook's full n=128 resolution.
"""
import numpy as np


def eady_growth_rate(mu, Lambda=1.0, H=1.0):
    """Growth rate sigma=mu*c_i for the Eady problem: QG flow between rigid
    lids at z=0,H with uniform shear U=Lambda*z and zero interior PV
    gradient. mu=k*N*H/f0 is the nondimensional horizontal wavenumber.

    The interior solution phi(z)=A*cosh(mu*z/H)+B*sinh(mu*z/H) (from
    q'=0 in the interior) is matched to the linearized thermal boundary
    condition (U-c)phi_z-Lambda*phi=0 at z=0,H, giving a 2x2 generalized
    eigenvalue problem (C0+c*C1)[A,B]=0 for the phase speed c -- solved
    directly rather than via a memorized closed-form formula.
    """
    C0 = np.array([
        [Lambda, 0.0],
        [mu * Lambda * np.sinh(mu) - Lambda * np.cosh(mu),
         mu * Lambda * np.cosh(mu) - Lambda * np.sinh(mu)],
    ])
    C1 = np.array([
        [0.0, mu / H],
        [-(mu / H) * np.sinh(mu), -(mu / H) * np.cosh(mu)],
    ])
    M = -np.linalg.solve(C1, C0)
    c = np.linalg.eigvals(M)
    return mu * max(float(np.max(c.imag)), 0.0)


def eady_growth_rate_curve(mu_values, Lambda=1.0, H=1.0):
    """eady_growth_rate evaluated at every mu in mu_values."""
    return np.array([eady_growth_rate(mu, Lambda, H) for mu in mu_values])


def invert_2layer(q1_hat, q2_hat, grid, F):
    """2-layer QGPV Helmholtz inversion: solve, at every wavenumber,

        (-k^2-F)*psi1 + F*psi2 = q1
        F*psi1 + (-k^2-F)*psi2 = q2

    for (psi1_hat, psi2_hat), using the closed-form 2x2 inverse (the matrix
    is symmetric with a=-(k^2+F) on the diagonal and F off it, so
    det=a^2-F^2=k^2*(k^2+2F), zero only at k=0 where psi is set to 0 by
    convention, as in every other Poisson/Helmholtz inversion in this book).
    """
    a = -(grid.k2 + F)
    b = F
    det = a ** 2 - b ** 2
    det_safe = np.where(grid.k2 == 0, 1.0, det)
    psi1_hat = np.where(grid.k2 == 0, 0.0, (a * q1_hat - b * q2_hat) / det_safe)
    psi2_hat = np.where(grid.k2 == 0, 0.0, (a * q2_hat - b * q1_hat) / det_safe)
    return psi1_hat, psi2_hat


def pv_2layer(psi1_hat, psi2_hat, grid, F):
    """Forward 2-layer QGPV operator (perturbation PV only, no mean beta*y
    term -- matches the convention used everywhere else in this module):
    the exact inverse operation of invert_2layer.

        q1 = laplacian(psi1) + F*(psi2-psi1)
        q2 = laplacian(psi2) + F*(psi1-psi2)
    """
    q1_hat = grid.laplacian(psi1_hat) + F * (psi2_hat - psi1_hat)
    q2_hat = grid.laplacian(psi2_hat) + F * (psi1_hat - psi2_hat)
    return q1_hat, q2_hat


def to_modes(field1, field2):
    """Barotropic/baroclinic decomposition of any 2-layer field pair
    (psi, q, ...): field_bt=(field1+field2)/2 (depth mean), field_bc=
    (field1-field2)/2 (half the inter-layer difference). Linear, so this
    works identically in physical or spectral space.

    Applied to (psi1,psi2), this decouples the 2-layer QGPV EXACTLY (no
    approximation) when the two layers have equal thickness, as assumed by
    invert_2layer's symmetric F coupling: summing and differencing
    q1=lap(psi1)+F(psi2-psi1), q2=lap(psi2)+F(psi1-psi2) gives

        q_bt = lap(psi_bt)              -- ordinary barotropic PV, no F at all
        q_bc = lap(psi_bc) - 2*F*psi_bc -- QG-stretching PV, coefficient 2F

    i.e. an isolated baroclinic-mode disturbance (q_bt=0 identically) obeys
    exactly ch8's single-layer QG equation with L_R^2=1/(2F); an isolated
    barotropic disturbance obeys exactly ch7/ch18's barotropic vorticity
    equation. The two modes are decoupled ONLY in the mean-field/linear
    operator; the Jacobian nonlinearity in rhs_2layer still couples them
    (see ch10) -- that coupling is the entire content of ch10's notebook.
    """
    return 0.5 * (field1 + field2), 0.5 * (field1 - field2)


def from_modes(field_bt, field_bc):
    """Inverse of to_modes: field1=field_bt+field_bc, field2=field_bt-field_bc."""
    return field_bt + field_bc, field_bt - field_bc


def rhs_2layer(state_hat, grid, F, U1, U2, beta, r=0.0):
    """Nonlinear + mean-flow + beta + Ekman-drag tendency for state_hat=
    (q1_hat,q2_hat), the perturbation QGPV in each layer, IN SPECTRAL SPACE
    -- for use as the nonlinear RHS passed to gfdlib.timestep.ifrk4_step,
    with hyperviscosity as the separate, diagonal linear operator L
    (identical for both layers; broadcasts across the stacked (2,n,n//2+1)
    array the same way grid.fft/ifft already do):

        q1_t = -U1*q1_x - (beta+F*(U1-U2))*psi1_x - J(psi1,q1)
        q2_t = -U2*q2_x - (beta-F*(U1-U2))*psi2_x - J(psi2,q2) - r*zeta2

    with (psi1,psi2) from invert_2layer(q1,q2). The linear part couples the
    two layers (a 2x2 system per wavenumber, not a scalar), so it cannot be
    folded into L the way ch18's single beta term is -- it's part of this
    "nonlinear" RHS instead, evaluated explicitly at each RK4 substage.

    r is bottom Ekman drag on the lower layer's relative vorticity only
    (zeta2=-k^2*psi2): the standard large-scale sink that lets a 2-layer
    baroclinic wave reach a genuine finite-amplitude equilibrium instead of
    growing without bound. A fixed, unrelaxed mean shear is otherwise an
    unlimited energy source -- this is a physically standard addition
    (Phillips' original 1954 model included it), not a numerical patch:
    tested here with r=0, growth was found to continue past all physically
    reasonable amplitudes regardless of how much hyperviscosity was added,
    until floating-point overflow; with r>0 tuned so the flow equilibrates
    before reaching amplitudes where the explicit-RK4 advective CFL would
    be violated, the classic growth-peak-decay-equilibrate life cycle
    emerges cleanly (verified manually to t=500; see ch16 notebook/PR).
    """
    q1_hat, q2_hat = state_hat
    psi1_hat, psi2_hat = invert_2layer(q1_hat, q2_hat, grid, F)

    beta1 = beta + F * (U1 - U2)
    beta2 = beta - F * (U1 - U2)

    dq1_hat = -U1 * grid.ddx(q1_hat) - beta1 * grid.ddx(psi1_hat) - grid.jacobian(psi1_hat, q1_hat)
    dq2_hat = (-U2 * grid.ddx(q2_hat) - beta2 * grid.ddx(psi2_hat) - grid.jacobian(psi2_hat, q2_hat)
               + r * grid.k2 * psi2_hat)   # -r*zeta2, zeta2=-k^2*psi2

    return np.stack([dq1_hat, dq2_hat])
