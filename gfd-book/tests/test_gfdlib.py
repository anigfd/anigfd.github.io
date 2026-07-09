"""Numeric correctness of gfdlib primitives. `make test` must pass before commit."""
import numpy as np
from gfdlib import spectral, timestep, diagnostics, shallowwater, internalwaves


def test_poisson_roundtrip():
    g = spectral.Grid(64)
    psi = np.sin(2 * g.x) * np.cos(3 * g.y)
    zeta = g.ifft(g.laplacian(g.fft(psi)))
    rec = g.ifft(g.invert_laplacian(g.fft(zeta)))
    rec -= rec.mean(); psi -= psi.mean()
    assert np.max(np.abs(rec - psi)) < 1e-10


def test_jacobian_self_zero():
    g = spectral.Grid(64)
    A = g.fft(np.sin(g.x) * np.sin(2 * g.y))
    assert np.max(np.abs(g.ifft(g.jacobian(A, A)))) < 1e-10


def test_rk4_decay():
    y = 1.0
    for _ in range(100):
        y = timestep.rk4(lambda t, y: -y, y, 0.01)
    assert abs(y - np.exp(-1)) < 1e-6


def test_spectrum_parseval():
    g = spectral.Grid(64)
    f = np.random.default_rng(0).standard_normal((64, 64)); f -= f.mean()
    _, E = diagnostics.isotropic_spectrum(g.fft(f), g)
    assert abs(E.sum() - (f ** 2).mean()) < 1e-6


def test_ifrk4_exact_linear_rotation():
    """With no nonlinear term, the integrating factor is the exact propagator."""
    F = np.array(1.0 + 0j)
    L = 1j * 2.0
    for _ in range(100):
        F = timestep.ifrk4_step(F, lambda t, F: 0 * F, 0.01, L)
    assert abs(F - np.exp(1j * 2.0)) < 1e-12


def test_rossby_phase_speed():
    """Linear beta term: a k=(1,0) mode moves westward at c = -beta/k^2."""
    g = spectral.Grid(64)
    beta = 5.0
    zh = g.fft(np.cos(g.x))
    L = 1j * beta * g.kx * g.k2_inv          # -beta*v with psi = invert_laplacian(zeta)
    dt, nsteps = 0.01, 100
    for _ in range(nsteps):
        zh = timestep.ifrk4_step(zh, lambda t, F: 0 * F, dt, L)
    expected = np.cos(g.x + beta * dt * nsteps)   # cos(x - c*t), c = -beta
    assert np.max(np.abs(g.ifft(zh) - expected)) < 1e-8


def test_selective_decay():
    """2D turbulence: enstrophy must decay faster than energy."""
    g = spectral.Grid(64)
    rng = np.random.default_rng(1)
    W = g.fft(rng.standard_normal((64, 64))) * np.exp(-g.k2 / (2 * 6.0 ** 2))
    W *= g.dealias
    rhs = lambda t, W: -g.jacobian(g.invert_laplacian(W), W)
    ke0, ens0 = diagnostics.energy_enstrophy(g.invert_laplacian(W), g)
    for _ in range(1000):
        W = timestep.ifrk4_step(W, rhs, 1e-3, -2e-3 * g.k2) * g.dealias
    ke1, ens1 = diagnostics.energy_enstrophy(g.invert_laplacian(W), g)
    assert (ens1 / ens0) < (ke1 / ke0)


def test_sw_geostrophic_mode_is_steady():
    """The exact geostrophic mode u=-eta_y, v=eta_x is a fixed point of rhs."""
    g = spectral.Grid(32, L=20.0)
    eta = np.sin(2 * np.pi * g.x / g.L) * np.cos(2 * np.pi * g.y / g.L)
    eta_hat = g.fft(eta)
    u = g.ifft(-g.ddy(eta_hat))
    v = g.ifft(g.ddx(eta_hat))
    f = shallowwater.coriolis(g, beta=0.0)
    d = shallowwater.rhs(np.stack([eta, u, v]), g, f)
    assert np.max(np.abs(d)) < 1e-9


def test_sw_invert_pv_roundtrip():
    """invert_pv is the exact inverse of forming q from a geostrophic eta."""
    g = spectral.Grid(32, L=20.0)
    eta_true = np.sin(2 * np.pi * g.x / g.L) * np.cos(4 * np.pi * g.y / g.L)
    zeta_g = g.ifft(g.laplacian(g.fft(eta_true)))   # geostrophic vorticity
    q = zeta_g - eta_true
    eta_rec = shallowwater.invert_pv(q, g)
    assert np.max(np.abs(eta_rec - eta_true)) < 1e-10


def test_sw_pv_conserved_on_f_plane():
    """Released from rest, q is pointwise conserved when beta=0."""
    g = spectral.Grid(32, L=20.0)
    eta0 = 0.5 * np.exp(-((g.x - g.L / 2) ** 2 + (g.y - g.L / 2) ** 2) / (2 * 2.0 ** 2))
    state = np.stack([eta0, np.zeros_like(eta0), np.zeros_like(eta0)])
    f = shallowwater.coriolis(g, beta=0.0)
    q0 = shallowwater.potential_vorticity(state, g)
    dt = 0.4 * g.dx
    for _ in range(50):
        state = timestep.rk4(lambda t, s: shallowwater.rhs(s, g, f), state, dt)
    q1 = shallowwater.potential_vorticity(state, g)
    assert np.max(np.abs(q1 - q0)) < 1e-3


def test_sw_beta_sources_pv_at_expected_rate():
    """On a beta-plane, dq/dt = -beta*v at t=0 for a field released from rest."""
    g = spectral.Grid(32, L=20.0)
    eta0 = 0.5 * np.exp(-((g.x - g.L / 2) ** 2 + (g.y - g.L / 2) ** 2) / (2 * 2.0 ** 2))
    state0 = np.stack([eta0, np.zeros_like(eta0), np.zeros_like(eta0)])
    beta = 0.1
    f = shallowwater.coriolis(g, beta=beta)
    dt = 1e-4
    q0 = shallowwater.potential_vorticity(state0, g)
    state1 = timestep.rk4(lambda t, s: shallowwater.rhs(s, g, f), state0, dt)
    q1 = shallowwater.potential_vorticity(state1, g)
    dqdt_numeric = (q1 - q0) / dt
    v0 = state0[2]
    assert np.max(np.abs(dqdt_numeric - (-beta * v0))) < 1e-2


def test_iw_dispersion_bounded_by_f_and_N():
    """omega(k,m) is a convex combination of N^2 and f^2 -- always in [f,N]."""
    N, f = 1.0, 0.2
    rng = np.random.default_rng(0)
    kx, kz = rng.uniform(0.1, 5, 200), rng.uniform(0.1, 5, 200)
    omega = internalwaves.dispersion_omega(kx, kz, N, f)
    assert np.all(omega >= f - 1e-12) and np.all(omega <= N + 1e-12)


def test_iw_beam_angle_roundtrip():
    """beam_angle inverts the omega(theta) relation it's derived from."""
    N, f = 1.0, 0.15
    theta = np.linspace(0.05, np.pi / 2 - 0.05, 20)
    omega = np.sqrt(N ** 2 * np.cos(theta) ** 2 + f ** 2 * np.sin(theta) ** 2)
    assert np.max(np.abs(internalwaves.beam_angle(omega, N, f) - theta)) < 1e-10


def test_iw_leapfrog_matches_dispersion():
    """Constant N,f: a single Fourier mode oscillates at exactly omega(k,m)."""
    g = spectral.Grid(64)   # default L=2*pi
    N, f = 1.0, 0.1
    kx, kz = 3.0, 2.0
    omega = internalwaves.dispersion_omega(kx, kz, N, f)
    N2 = np.full_like(g.x, N ** 2)
    phase = kx * g.x + kz * g.y
    dt = (2 * np.pi / omega) / 400   # 400 steps/period for leapfrog accuracy
    q, q_prev = np.cos(phase), np.cos(phase) * np.cos(omega * dt)
    nsteps = 200
    for _ in range(nsteps):
        q, q_prev = internalwaves.step_leapfrog(q, q_prev, g, N2, f, 0.0, dt), q
    expected = np.cos(phase) * np.cos(omega * nsteps * dt)
    assert np.max(np.abs(q - expected)) < 1e-3
