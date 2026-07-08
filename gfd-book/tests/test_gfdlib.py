"""Numeric correctness of gfdlib primitives. `make test` must pass before commit."""
import numpy as np
from gfdlib import spectral, timestep, diagnostics, convection


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


def test_channel_poisson_matches_laplacian():
    """ChannelGrid.poisson_solve(rhs) inverts laplacian() to within FD/spectral error."""
    g = convection.ChannelGrid(nx=16, nz=12, Lx=2 * np.sqrt(2))
    rng = np.random.default_rng(0)
    rhs = rng.standard_normal((16, 12))
    f = g.poisson_solve(rhs)
    assert np.max(np.abs(g.laplacian(f) - rhs)) < 1e-8


def test_neutral_curve_minimum_is_critical_ra():
    """neutral_ra(k) is minimized at K_C with value RA_C."""
    k = np.linspace(0.5, 5, 2000)
    assert abs(k[np.argmin(convection.neutral_ra(k))] - convection.K_C) < 1e-2
    assert abs(convection.neutral_ra(convection.K_C) - convection.RA_C) < 1e-8


def test_convection_subcritical_decays_supercritical_grows():
    """Below Ra_c, a noise perturbation decays; well above it, it grows."""
    nx, nz, Lx, Pr = 32, 16, 2 * np.sqrt(2), 1.0
    g = convection.ChannelGrid(nx, nz, Lx)
    rng = np.random.default_rng(7)
    theta0 = 1e-3 * rng.standard_normal((nx, nz))
    dt = 0.2 * min(g.dx, g.dz) ** 2

    def run(Ra, nsteps):
        state = np.stack([np.zeros((nx, nz)), theta0.copy()])
        rhs = lambda t, s: convection.rhs_boussinesq(s, g, Ra, Pr)
        amp0 = np.abs(state[1]).max()
        for _ in range(nsteps):
            state = timestep.rk4(rhs, state, dt)
        return amp0, np.abs(state[1]).max()

    amp0, amp1 = run(Ra=500.0, nsteps=800)     # below RA_C ~= 657.5
    assert amp1 < 0.5 * amp0

    amp0, amp1 = run(Ra=5000.0, nsteps=800)    # well above RA_C
    assert amp1 > 5.0 * amp0


def test_lorenz_convective_fixed_point():
    """(sqrt(b(r-1)), sqrt(b(r-1)), r-1) is a fixed point of lorenz_rhs for r>1."""
    sigma, b, r = 10.0, 8.0 / 3.0, 28.0
    c = np.sqrt(b * (r - 1))
    state = np.array([c, c, r - 1])
    assert np.max(np.abs(convection.lorenz_rhs(state, sigma, r, b))) < 1e-10
