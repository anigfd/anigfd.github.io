"""Point vortices: the Kirchhoff Hamiltonian system, its invariants, and a
symplectic integrator to compare against the book's standard RK4 (pure
NumPy).

N point vortices with circulations Gamma_i at positions (x_i, y_i) are an
EXACT finite-dimensional reduction of 2D Euler flow (the vorticity field is
a sum of delta functions; each vortex is advected by the others' velocity
fields, never its own). Kirchhoff (1876): the system is Hamiltonian with
the unusual feature that x_i and y_i are themselves the conjugate pair --
weighted by circulation:

    Gamma_i dx_i/dt =  dH/dy_i,   Gamma_i dy_i/dt = -dH/dx_i,
    H = -(1/4pi) sum_{i<j} Gamma_i Gamma_j ln(r_ij^2).

Phase space IS physical space: a snapshot of vortex positions is a full
phase-space state, and area (weighted by Gamma) is the symplectic measure.

Noether's theorem then reads off the invariants from the symmetries of H
(it depends only on separations): translation in x -> Q = sum Gamma_i x_i,
translation in y -> P = sum Gamma_i y_i, rotation -> L = sum Gamma_i
(x_i^2+y_i^2), time -> H itself. (Conservation of Q and P can also be seen
directly: sum_i Gamma_i u_i is a double sum antisymmetric under i<->j.)

step_midpoint is the implicit midpoint rule -- symplectic for ANY
Hamiltonian, separable or not (this one is not: H mixes the conjugate
variables inseparably, which rules out the more familiar leapfrog).
Second-order accurate, solved here by fixed-point iteration. The payoff is
qualitative, not quantitative: a symplectic method's energy error stays
BOUNDED forever (backward-error analysis: it exactly integrates a nearby
Hamiltonian), while any non-symplectic method's energy drifts secularly --
verified in the tests by comparing against same-order Heun over a long
chaotic 4-vortex run.
"""
import numpy as np


def vortex_rhs(state, Gamma):
    """d(state)/dt for state = array (2, N) of positions, from the mutual
    Biot-Savart sum (each vortex advected by all OTHERS):

        u_i = -(1/2pi) sum_{j!=i} Gamma_j (y_i-y_j)/r_ij^2
        v_i = +(1/2pi) sum_{j!=i} Gamma_j (x_i-x_j)/r_ij^2
    """
    x, y = state
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    r2 = dx ** 2 + dy ** 2
    np.fill_diagonal(r2, 1.0)          # dummy; self-term removed next line
    inv_r2 = 1.0 / r2
    np.fill_diagonal(inv_r2, 0.0)
    u = -(1.0 / (2 * np.pi)) * (dy * inv_r2) @ Gamma
    v = (1.0 / (2 * np.pi)) * (dx * inv_r2) @ Gamma
    return np.stack([u, v])


def hamiltonian(state, Gamma):
    """H = -(1/4pi) sum_{i<j} Gamma_i Gamma_j ln(r_ij^2), the conserved
    interaction energy (Kirchhoff). Verified in tests to generate exactly
    vortex_rhs through Gamma_i dx_i/dt = dH/dy_i, Gamma_i dy_i/dt =
    -dH/dx_i, by conservation along integrated trajectories."""
    x, y = state
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    r2 = dx ** 2 + dy ** 2
    np.fill_diagonal(r2, 1.0)
    log_r2 = np.log(r2)                 # diagonal contributes log(1)=0
    return -(1.0 / (8 * np.pi)) * Gamma @ log_r2 @ Gamma   # 1/8: double-counts i<j


def invariants(state, Gamma):
    """(H, Q, P, L): energy, the two linear impulses Q = sum Gamma*x,
    P = sum Gamma*y, and the angular impulse L = sum Gamma*(x^2+y^2) --
    the Noether charges of translation, rotation and time symmetry."""
    x, y = state
    return (hamiltonian(state, Gamma),
            float(Gamma @ x), float(Gamma @ y),
            float(Gamma @ (x ** 2 + y ** 2)))


def step_heun(state, Gamma, dt):
    """Heun's method (explicit RK2) -- 2nd order, NOT symplectic. Included
    as the fair same-order, same-cost comparison for step_midpoint."""
    k1 = vortex_rhs(state, Gamma)
    k2 = vortex_rhs(state + dt * k1, Gamma)
    return state + 0.5 * dt * (k1 + k2)


def step_midpoint(state, Gamma, dt, n_iter=6):
    """Implicit midpoint rule, s' = s + dt*f((s+s')/2) -- 2nd order and
    symplectic for any Hamiltonian. The implicit equation is solved by
    fixed-point iteration from an explicit-Euler predictor; for the time
    steps used in this book the iteration converges to machine precision
    well within n_iter=6 sweeps (contraction factor ~ dt * |grad f|)."""
    s_new = state + dt * vortex_rhs(state, Gamma)
    for _ in range(n_iter):
        s_new = state + dt * vortex_rhs(0.5 * (state + s_new), Gamma)
    return s_new
