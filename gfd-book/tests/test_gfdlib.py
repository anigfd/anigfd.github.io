"""Numeric correctness of gfdlib primitives. `make test` must pass before commit."""
import numpy as np
from gfdlib import spectral, timestep, diagnostics, shallowwater, internalwaves, convection, pv, balance, rossby, qg, instability, baroclinic


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


def test_rossby_stationary_wavenumber_gives_zero_omega():
    """At K=K_s, the local dispersion relation is exactly stationary
    (omega=0) at every latitude -- K_s's latitude-independence is what
    makes this possible for all phi simultaneously."""
    Omega_s = 0.4
    Ks2 = rossby.stationary_wavenumber2(Omega_s)
    Ks = np.sqrt(Ks2)
    for phi_deg in [10, 30, 50, 70]:
        phi = np.radians(phi_deg)
        for alpha in np.linspace(0, 2 * np.pi, 9):
            n = Ks * np.cos(alpha) * np.cos(phi)
            m = Ks * np.sin(alpha)
            omega = rossby.dispersion_omega(phi, n, m, Omega_s)
            assert abs(omega) < 1e-10


def test_rossby_ray_is_exact_great_circle():
    """Integrating the ray equations for a stationary wave on solid-body
    background rotation produces an exact great circle (Hoskins & Karoly
    1981) -- verified via a coordinate-free planarity check, not assumed."""
    Omega_s = 0.3
    lam0, phi0, alpha0 = 0.0, np.radians(30.0), np.radians(35.0)
    state, n = rossby.launch_state(lam0, phi0, alpha0, Omega_s)

    dt, nsteps = 0.01, 3000
    traj = [state.copy()]
    for _ in range(nsteps):
        state = timestep.rk4(lambda t, s: rossby.ray_rhs(s, n, Omega_s), state, dt)
        traj.append(state.copy())
    traj = np.array(traj)

    dev = rossby.great_circle_deviation(traj[:, 0], traj[:, 1])
    assert np.max(dev) < 1e-6


def test_rossby_zonal_wavenumber_index_conserved():
    """n is conserved along a ray to high precision (the numerical-diff
    ray equations should reproduce the exact zonal symmetry dn/dt=0)."""
    Omega_s = 0.3
    state, n = rossby.launch_state(0.0, np.radians(20.0), np.radians(50.0), Omega_s)
    dt = 0.01
    for _ in range(500):
        state = timestep.rk4(lambda t, s: rossby.ray_rhs(s, n, Omega_s), state, dt)
    # n itself isn't part of state (it's fixed by construction); instead
    # check the physical zonal wavenumber kx=n/(a cos phi) is consistent
    # with a still-stationary ray (omega should remain ~0 throughout).
    lam, phi, m = state
    omega_final = rossby.dispersion_omega(phi, n, m, Omega_s)
    assert abs(omega_final) < 1e-6


def test_dbdy_matches_finite_difference():
    """dbdy_frontal is the exact y-derivative of buoyancy_field's frontal part."""
    db, Ly, Htrop = 1.0, 1.0, 1.0
    y0, z0, dy = 1.3, 0.6, 1e-6
    fd = (balance.buoyancy_field(y0 + dy, z0, db, Ly, Htrop)
          - balance.buoyancy_field(y0 - dy, z0, db, Ly, Htrop)) / (2 * dy)
    exact = balance.dbdy_frontal(y0, z0, db, Ly, Htrop)
    assert abs(fd - exact) < 1e-6


def test_thermal_wind_relation_satisfied():
    """f*du_g/dz == -dbdy_frontal, the thermal-wind relation thermal_wind_u integrates."""
    db, Ly, Htrop, f = 1.0, 1.0, 1.0, 0.5
    y0, z0, dz = 1.3, 0.9, 1e-6
    dudz_fd = (balance.thermal_wind_u(y0, z0 + dz, db, Ly, Htrop, f)
               - balance.thermal_wind_u(y0, z0 - dz, db, Ly, Htrop, f)) / (2 * dz)
    lhs = f * dudz_fd
    rhs = -balance.dbdy_frontal(y0, z0, db, Ly, Htrop)
    assert abs(lhs - rhs) < 1e-5


def test_jet_core_at_tropopause():
    """u_g(y,z) has an extremum in z exactly at z=Htrop, for any fixed y."""
    db, Ly, Htrop, f = 1.0, 1.0, 1.0, 0.5
    y0 = 0.7
    z = np.linspace(0.01, 2 * Htrop - 0.01, 4001)
    u = balance.thermal_wind_u(y0, z, db, Ly, Htrop, f)
    z_at_extremum = z[np.argmax(np.abs(u))]
    assert abs(z_at_extremum - Htrop) < 1e-3


def test_gaussian_blob_periodic_at_edge():
    """A blob placed at the domain edge wraps continuously, no clipping."""
    g = spectral.Grid(64, L=10.0)
    b1 = pv.gaussian_blob(g, x0=0.0, y0=5.0, amp=1.0, sigma=1.0)
    b2 = pv.gaussian_blob(g, x0=10.0, y0=5.0, amp=1.0, sigma=1.0)
    assert np.max(np.abs(b1 - b2)) < 1e-10
    assert b1.max() > 0.99   # peak amplitude is reached somewhere (not clipped)


def test_gaussian_blob_peak_at_center():
    g = spectral.Grid(64, L=10.0)
    x0, y0 = 20 * g.dx, 30 * g.dx   # grid-aligned, so the peak is exact
    b = pv.gaussian_blob(g, x0=x0, y0=y0, amp=2.0, sigma=0.8)
    iy, ix = np.unravel_index(np.argmax(b), b.shape)
    assert abs(g.x[iy, ix] - x0) < 1e-10 and abs(g.y[iy, ix] - y0) < 1e-10
    assert abs(b.max() - 2.0) < 1e-10


def test_staircase_jets_align_with_pv_risers():
    """Inverting a PV staircase: u extrema sit at PV-gradient maxima, not
    within the well-mixed plateaus (the point of a PV staircase)."""
    g = spectral.Grid(256)
    n_steps = 3
    q = pv.staircase_pv(g, n_steps=n_steps, amp=1.0)
    psi_hat = g.invert_laplacian(g.fft(q))
    u = g.ifft(-g.ddy(psi_hat))

    q_prof, u_prof, y = q[0, :], u.mean(axis=0), g.y[0, :]
    dqdy = np.gradient(q_prof, y)

    n = len(u_prof)
    extrema = [i for i in range(n)
               if (u_prof[i] > u_prof[i - 1] and u_prof[i] > u_prof[(i + 1) % n])
               or (u_prof[i] < u_prof[i - 1] and u_prof[i] < u_prof[(i + 1) % n])]
    assert len(extrema) == 2 * n_steps   # one jet per riser, alternating sign

    mean_absdqdy = np.mean(np.abs(dqdy))
    for i in extrema:
        assert abs(dqdy[i]) > 10 * mean_absdqdy


def test_vortex_dipole_propagates():
    """An opposite-signed vortex pair self-advects (translates), unlike a
    single vortex or a like-signed pair (which merely co-rotate in place)."""
    g = spectral.Grid(64)
    xc = yc = 0.5 * g.L
    q = (pv.gaussian_blob(g, xc, yc - 0.5, 1.0, 0.35)
         + pv.gaussian_blob(g, xc, yc + 0.5, -1.0, 0.35))
    zeta_hat = g.fft(q) * g.dealias
    L_op = -1e-4 * g.k2 ** 2
    rhs = lambda t, zh: -g.jacobian(g.invert_laplacian(zh), zh)

    x0 = g.x[np.unravel_index(np.argmax(np.abs(q)), q.shape)]
    dt = 0.01
    for _ in range(400):
        zeta_hat = timestep.ifrk4_step(zeta_hat, rhs, dt, L_op) * g.dealias
    zeta_final = g.ifft(zeta_hat)
    x1 = g.x[np.unravel_index(np.argmax(np.abs(zeta_final)), zeta_final.shape)]
    assert abs(((x1 - x0 + g.L / 2) % g.L) - g.L / 2) > 0.3   # net translation


def test_qg_dispersion_bounded_by_beta_over_2():
    """QG Rossby-wave frequency is capped at beta/2, attained at |k|=1 --
    unlike the barotropic relation, which is unbounded as k->0."""
    beta = 3.0
    rng = np.random.default_rng(0)
    kx, ky = rng.uniform(-5, 5, 2000), rng.uniform(-5, 5, 2000)
    omega = qg.dispersion_omega(kx, ky, beta)
    assert np.all(np.abs(omega) <= beta / 2 + 1e-10)

    k1 = np.linspace(0.01, 5, 5000)
    omega_k1 = qg.dispersion_omega(k1, 0.0, beta)
    assert abs(k1[np.argmax(np.abs(omega_k1))] - 1.0) < 1e-2
    assert abs(np.max(np.abs(omega_k1)) - beta / 2) < 1e-3


def test_qg_helmholtz_inversion_matches_shallowwater():
    """The QGPV Helmholtz operator (L_R_hat=1) is literally
    shallowwater.invert_pv -- verify it actually inverts (nabla^2-1)psi=q."""
    g = spectral.Grid(32, L=20.0)
    psi_true = np.sin(2 * np.pi * g.x / g.L) * np.cos(4 * np.pi * g.y / g.L)
    q = g.ifft(g.laplacian(g.fft(psi_true))) - psi_true
    psi_rec = shallowwater.invert_pv(q, g)
    assert np.max(np.abs(psi_rec - psi_true)) < 1e-10


def test_qg_vortex_drifts_west_on_beta_plane():
    """An isolated QG vortex on a beta-plane sheds a Rossby wake and drifts
    westward (a beta-gyre), just as in the barotropic case but screened by
    the deformation radius via the QGPV Helmholtz inversion."""
    n, L = 96, 20.0
    g = spectral.Grid(n, L=L)
    xc = yc = 0.5 * L
    q0 = pv.gaussian_blob(g, xc, yc, 1.0, 1.0)

    beta, nu, nnu = 0.3, 1e-4, 2
    q_hat = g.fft(q0)
    L_op = -nu * g.k2 ** nnu + 1j * beta * g.kx / (g.k2 + 1.0)

    def rhs_nl(t, qh):
        psi_h = -qh / (g.k2 + 1.0)
        return -g.jacobian(psi_h, qh)

    dt = 0.01
    for _ in range(2000):
        q_hat = timestep.ifrk4_step(q_hat, rhs_nl, dt, L_op) * g.dealias

    q_final = g.ifft(q_hat)
    w = np.clip(q_final, 0, None)
    cx_final = (g.x * w).sum() / w.sum()
    drift = ((cx_final - xc + L / 2) % L) - L / 2
    assert drift < -0.5   # net westward drift


def test_eady_growth_rate_matches_published_benchmark():
    """The Eady growth-rate curve's peak location, peak value, and cutoff
    match the published benchmark (Vallis/Pedlosky): mu_max~1.61,
    sigma_max~0.31, cutoff mu_c~2.399 -- independent numbers this code was
    not tuned to reproduce."""
    mu_values = np.linspace(0.05, 3.0, 300)
    growth = baroclinic.eady_growth_rate_curve(mu_values)
    mu_max = mu_values[np.argmax(growth)]
    assert abs(mu_max - 1.61) < 0.05
    assert abs(growth.max() - 0.31) < 0.02
    unstable = mu_values[growth > 1e-6]
    assert abs(unstable.max() - 2.399) < 0.05


def test_2layer_inversion_roundtrip():
    """invert_2layer exactly inverts the forward 2-layer Helmholtz operator."""
    g = spectral.Grid(32)
    rng = np.random.default_rng(0)
    psi1 = rng.standard_normal((32, 32)); psi1 -= psi1.mean()
    psi2 = rng.standard_normal((32, 32)); psi2 -= psi2.mean()
    F = 2.0
    psi1_hat, psi2_hat = g.fft(psi1), g.fft(psi2)
    q1_hat = (-g.k2 - F) * psi1_hat + F * psi2_hat
    q2_hat = F * psi1_hat + (-g.k2 - F) * psi2_hat
    psi1_rec_hat, psi2_rec_hat = baroclinic.invert_2layer(q1_hat, q2_hat, g, F)
    assert np.max(np.abs(g.ifft(psi1_rec_hat) - psi1)) < 1e-10
    assert np.max(np.abs(g.ifft(psi2_rec_hat) - psi2)) < 1e-10


def test_2layer_nonlinear_growth_matches_linear_theory():
    """A small-amplitude single-wavenumber perturbation's time-stepped
    growth rate (rhs_2layer, spectral, via ifrk4_step with zero
    hyperviscosity) converges to an independently derived linear
    eigenvalue, once the initial condition purifies onto the dominant
    eigenmode."""
    F, U1, U2, beta = 1.3, 0.6, -0.6, 0.3

    def linear_growth_rate(k):
        beta1, beta2 = beta + F * (U1 - U2), beta - F * (U1 - U2)
        a, b = -(k ** 2 + F), F
        C0 = np.array([[U1 * a + beta1, U1 * b], [U2 * b, U2 * a + beta2]])
        C1 = np.array([[-a, -b], [-b, -a]])
        c = np.linalg.eigvals(-np.linalg.solve(C1, C0))
        return k * max(float(np.max(c.imag)), 0.0)

    k_int = 1
    sigma_pred = linear_growth_rate(float(k_int))
    assert sigma_pred > 0.1   # confirm this wavenumber is genuinely unstable

    g = spectral.Grid(48)
    eps = 1e-6
    state = np.stack([eps * np.cos(k_int * g.x), -eps * np.cos(k_int * g.x)])
    state_hat = g.fft(state)
    L_op = np.zeros_like(g.k2)   # no dissipation needed for this short, small-amplitude run
    dt = 0.005
    nsteps = int(16 / dt)
    amps, ts = [], []
    t = 0.0
    for s in range(nsteps + 1):
        if s % max(1, nsteps // 80) == 0:
            amps.append(np.abs(g.ifft(state_hat[0])).max()); ts.append(t)
        state_hat = timestep.ifrk4_step(
            state_hat, lambda tt, sh: baroclinic.rhs_2layer(sh, g, F, U1, U2, beta), dt, L_op)
        t += dt
    amps, ts = np.array(amps), np.array(ts)
    mask = (ts > 12) & (ts < 16)   # late window: mode has purified by then
    sigma_measured = np.polyfit(ts[mask], np.log(amps[mask]), 1)[0]
    assert abs(sigma_measured - sigma_pred) / sigma_pred < 0.01


def test_instability_no_inflection_point_is_stable():
    """Rayleigh's necessary condition: no sign change in U'' (beta=0) means
    every tested k must be exactly neutral (zero growth)."""
    L, n = 10.0, 200
    y = np.linspace(-L / 2, L / 2, n + 2)[1:-1]
    U = -np.cos(np.pi * y / L)   # U'' <= 0 strictly in the interior
    for k in [0.3, 0.6, 1.0, 1.5, 2.0]:
        assert instability.growth_rate(y, U, k, beta=0.0) == 0.0


def test_instability_tanh_shear_layer_peak_matches_michalke():
    """The classical U=tanh(y/delta) shear layer is unstable, with peak
    growth rate at k*delta ~= 0.4446 (Michalke 1964) -- an independent
    published number, not a value tuned to match this code."""
    delta = 1.0
    n = 400
    y = np.linspace(-15 * delta, 15 * delta, n + 2)[1:-1]
    U = np.tanh(y / delta)

    k_values = np.linspace(0.05, 1.2, 48)
    growth = instability.growth_rate_curve(y, U, k_values, beta=0.0)

    assert growth[-1] == 0.0            # short waves are stable (k*delta=1.2)
    assert growth.max() > 0.1           # genuinely unstable in between
    k_peak = k_values[np.argmax(growth)]
    assert abs(k_peak * delta - 0.4446) < 0.03


def test_double_shear_layer_vorticity_signs():
    """The two jets carry opposite-signed vorticity, peaked at their
    respective centers, consistent with two counter-propagating shears."""
    g = spectral.Grid(64)
    delta = g.L / 40
    zeta = instability.double_shear_layer(g, delta, v_pert=0.0)   # no perturbation for this check
    y1_idx = np.argmin(np.abs(g.y[0, :] - 0.25 * g.L))
    y2_idx = np.argmin(np.abs(g.y[0, :] - 0.75 * g.L))
    assert zeta[0, y1_idx] < 0
    assert zeta[0, y2_idx] > 0
