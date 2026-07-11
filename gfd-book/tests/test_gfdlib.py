"""Numeric correctness of gfdlib primitives. `make test` must pass before commit."""
import numpy as np
from gfdlib import spectral, timestep, diagnostics, shallowwater, internalwaves, convection, pv, balance, rossby, qg, instability, baroclinic, symmetric, mixing, circulation, overturning, kinematics, rotation, scaling


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


def test_abyssal_profile_bcs_and_limits():
    """Munk's abyssal-recipe profile satisfies its two BCs exactly, reduces
    to the linear pure-diffusion profile as w->0, and (for strong upward
    w) keeps the interior close to the BOTTOM value with the transition to
    the top value compressed into a thin layer near z=H."""
    w, kappa, H, T0, T1 = 1.0, 0.3, 1.0, 0.0, 1.0
    vals = overturning.abyssal_profile(np.array([0.0, H]), w, kappa, H, T0, T1)
    assert np.max(np.abs(vals - np.array([T0, T1]))) < 1e-10

    T_diffusive = overturning.abyssal_profile(np.array([0.5]), 1e-8, kappa, H, T0, T1)
    assert abs(T_diffusive[0] - 0.5) < 1e-3

    T_strong = overturning.abyssal_profile(np.array([0.3, 0.5, 0.7]), 20.0, kappa, H, T0, T1)
    assert np.all(T_strong < 0.01)


def test_abyssal_profile_satisfies_ode():
    """Independent finite-difference check: w*T'-kappa*T'' must vanish in
    the interior."""
    w, kappa, H = 1.0, 0.3, 1.0
    z = np.linspace(0, H, 2001)
    dz = z[1] - z[0]
    T = overturning.abyssal_profile(z, w, kappa, H, 0.0, 1.0)
    Tp = np.gradient(T, dz)
    Tpp = np.zeros_like(T)
    Tpp[1:-1] = (T[2:] - 2 * T[1:-1] + T[:-2]) / dz ** 2
    resid = w * Tp[1:-1] - kappa * Tpp[1:-1]
    assert np.max(np.abs(resid)) < 1e-3


def test_overturning_rises_at_heated_column_sinks_at_cooled():
    """The differential-heating-driven overturning cell must rise over the
    heated column and sink over the cooled one -- the defining physical
    signature of the Hadley-cell/MOC mechanism, not just numerical
    stability."""
    nx, nz, Lx = 32, 16, 2.0
    grid = convection.ChannelGrid(nx, nz, Lx)
    Ra, Pr, Q0 = 2e5, 1.0, 5.0
    Q = overturning.surface_heating(grid, Q0)

    rng = np.random.default_rng(0)
    zeta0 = np.zeros((nx, nz))
    theta0 = 1e-3 * rng.standard_normal((nx, nz))
    state = np.stack([zeta0, theta0])
    dt = 0.15 * min(grid.dx, grid.dz) ** 2
    rhs = lambda t, s: overturning.rhs_overturning(s, grid, Ra, Pr, Q)
    for _ in range(1500):
        state = timestep.rk4(rhs, state, dt)

    assert np.all(np.isfinite(state))
    zeta, theta = state
    psi = grid.poisson_solve(zeta)
    w = grid.ddx(psi)   # (u,w)=(-psi_z,psi_x), per convection.py's convention

    i_warm = 0                                        # grid.x[0] = 0
    i_cool = np.argmin(np.abs(grid.x - Lx / 2))        # cos(2pi x/Lx)=-1 there
    assert w[i_warm, nz // 2] > 0    # rising over the heated column
    assert w[i_cool, nz // 2] < 0    # sinking over the cooled column


def test_stommel_analytic_satisfies_bcs_and_sverdrup_limit():
    """The exact separable Stommel solution satisfies psi=0 at both x
    boundaries exactly, and recovers the classic Sverdrup interior
    solution f=1-x away from the western boundary layer as eps->0."""
    eps = 0.05
    assert abs(circulation.stommel_f(np.array([0.0, 1.0]), eps)).max() < 1e-8

    eps_small = 1e-4
    f_small = circulation.stommel_f(np.array([0.5, 0.8, 0.99]), eps_small)
    f_sverdrup = 1 - np.array([0.5, 0.8, 0.99])
    assert np.max(np.abs(f_small - f_sverdrup)) < 0.02


def test_stommel_fd_matches_analytic_with_2nd_order_convergence():
    """The general 2D finite-difference Stommel solver, cross-validated
    against the independently-derived exact analytic solution, must
    converge at 2nd order as resolution increases."""
    eps = 0.05
    errs = []
    for n in [40, 80]:
        psi_fd, x_full, y_full = circulation.solve_stommel_fd(n, n, eps)
        X, Y = np.meshgrid(x_full, y_full, indexing="xy")
        psi_exact = circulation.stommel_analytic(X, Y, eps)
        errs.append(np.max(np.abs(psi_fd - psi_exact)))
    # doubling resolution should cut error by ~4x (2nd order); allow slack
    assert errs[1] < errs[0] / 2.5
    assert errs[1] < 0.01


def test_munk_analytic_satisfies_noslip_bcs():
    """The exact separable Munk solution satisfies BOTH psi=0 and
    psi'=0 (no-slip) at both x boundaries exactly, verified via the exact
    closed-form derivative (not a finite-difference estimate, which is
    inaccurate at this sharp a boundary layer even on a fine grid)."""
    delta = 0.1
    C, roots, fp = circulation._munk_coeffs(delta)

    def fprime(x):
        return np.real(np.sum(C * roots * np.exp(roots * x)))

    for x in [0.0, 1.0]:
        assert abs(circulation.munk_f(np.array([x]), delta)[0]) < 1e-8
        assert abs(fprime(x)) < 1e-8


def test_munk_fd_matches_analytic_with_2nd_order_convergence():
    """The general 2D coupled (psi,zeta) finite-difference Munk solver
    (no-slip east/west via Thom's formula, free-slip north/south), cross-
    validated against the independently-derived exact analytic solution,
    must converge at 2nd order. (Getting the north/south BC wrong -- full
    no-slip on all four walls -- was tried and gave a resolution-
    INDEPENDENT ~30% error, immediately distinguishing a boundary-
    condition mismatch from genuine discretization error.)"""
    delta = 0.15
    errs = []
    for n in [30, 60]:
        psi_fd, x_full, y_full = circulation.solve_munk_fd(n, n, delta)
        X, Y = np.meshgrid(x_full, y_full, indexing="xy")
        psi_exact = circulation.munk_analytic(X, Y, delta)
        errs.append(np.max(np.abs(psi_fd - psi_exact)))
    assert errs[1] < errs[0] / 2.5
    assert errs[1] < 0.01


def _run_mixing(Gamma, nu=1e-6, nnu=2, kappa=1e-6, kappa_nnu=2, n=64,
                 nsteps=600, dt=0.01, k0=8, seed=1234):
    """Shared driver: decaying 2D turbulence (ch18's exact McWilliams IC)
    stirring a passive tracer perturbation against mean gradient Gamma."""
    grid = spectral.Grid(n)
    rng = np.random.default_rng(seed)
    K = np.where(grid.kmag == 0, 1e-10, grid.kmag)
    amp = K * np.sqrt(1.0 / (1.0 + (K / k0) ** 4))
    zeta_hat = grid.dealias * amp * np.exp(2j * np.pi * rng.random(K.shape))
    ke, _ = diagnostics.energy_enstrophy(grid.invert_laplacian(zeta_hat), grid)
    zeta_hat = zeta_hat * np.sqrt(0.5 / ke)

    state_hat = np.stack([zeta_hat, np.zeros_like(zeta_hat)])
    L_op = np.stack([-nu * grid.k2 ** nnu, -kappa * grid.k2 ** kappa_nnu])
    for _ in range(nsteps):
        state_hat = timestep.ifrk4_step(
            state_hat, lambda tt, sh: mixing.rhs_coupled(sh, grid, Gamma), dt, L_op
        ) * grid.dealias

    zh, ch = state_hat
    psi_hat = grid.invert_laplacian(zh)
    v = grid.ifft(grid.ddx(psi_hat))
    cprime = grid.ifft(ch)
    return mixing.effective_diffusivity(v, cprime, Gamma)


def test_mixing_keff_independent_of_gamma():
    """The tracer perturbation equation is LINEAR in c', so K_eff must be
    exactly independent of the imposed mean gradient Gamma for the same
    underlying flow -- a strong, non-phenomenological check."""
    K1 = _run_mixing(Gamma=1.0)
    K2 = _run_mixing(Gamma=2.0)
    assert abs(K1 - K2) / abs(K1) < 1e-8


def test_mixing_keff_zero_without_velocity():
    """With the velocity field identically zero, no v'c' correlation can
    develop: K_eff must be exactly zero."""
    n = 48
    grid = spectral.Grid(n)
    zeta_hat = np.zeros((n, n // 2 + 1), dtype=complex)
    state_hat = np.stack([zeta_hat, np.zeros_like(zeta_hat)])
    L_op = np.stack([np.zeros_like(grid.k2), -1e-3 * grid.k2])
    for _ in range(200):
        state_hat = timestep.ifrk4_step(
            state_hat, lambda tt, sh: mixing.rhs_coupled(sh, grid, 1.0), 0.01, L_op
        ) * grid.dealias
    zh, ch = state_hat
    psi_hat = grid.invert_laplacian(zh)
    v = grid.ifft(grid.ddx(psi_hat))
    cprime = grid.ifft(ch)
    assert abs(mixing.effective_diffusivity(v, cprime, 1.0)) < 1e-12


def test_mixing_keff_positive_for_decaying_turbulence():
    """Freely-decaying 2D turbulence stirring a mean gradient should
    produce down-gradient (K_eff>0) transport -- the generic expected
    result, checked here rather than assumed."""
    Keff = _run_mixing(Gamma=1.0, nsteps=900)
    assert Keff > 0.0


def test_symmetric_ro0_critical_ri_is_one():
    """At Ro=0 (no relative vorticity) the symmetric-instability marginal
    curve reduces to the textbook value Ri_c=1 -- an independent check on
    the sign convention, not a fitted number."""
    assert abs(symmetric.critical_ri(0.0) - 1.0) < 1e-12
    assert symmetric.ertel_pv_ratio(0.0, 0.9) < 0    # below critical -> unstable
    assert symmetric.ertel_pv_ratio(0.0, 1.1) > 0    # above critical -> stable


def test_symmetric_inertial_unconditional_for_ro_below_minus1():
    """(1+Ro)<0 means the absolute vertical vorticity has already changed
    sign: the front must be unstable at EVERY Ri>0, independent of
    stratification strength."""
    for Ri in [0.01, 1.0, 1e6]:
        assert symmetric.ertel_pv_ratio(-1.5, Ri) < 0


def test_symmetric_classify_regions():
    """classify labels each of the four regimes correctly at representative
    points, including the Ri<=0 gravitational branch and vectorized input."""
    Ro = np.array([0.0, 0.0, -1.5, 1.0])
    Ri = np.array([0.5, 2.0, 1.0, -0.1])
    labels = symmetric.classify(Ro, Ri)
    expected = [symmetric.SYMMETRIC, symmetric.STABLE,
                symmetric.INERTIAL, symmetric.GRAVITATIONAL]
    assert list(labels) == expected


def test_taylor_goldstein_matches_rayleigh_at_zero_stratification():
    """N^2=0 must reduce the Taylor-Goldstein solver exactly to the
    independent barotropic Rayleigh solver (growth_rate) -- the strongest
    available cross-check, since both are derived and implemented
    separately."""
    n = 400
    z = np.linspace(-10, 10, n + 2)[1:-1]
    U = np.tanh(z)
    k_values = np.linspace(0.05, 1.0, 30)
    g_tg = instability.taylor_goldstein_growth_rate_curve(z, U, np.zeros(n), k_values)
    g_ray = instability.growth_rate_curve(z, U, k_values, beta=0.0)
    assert np.max(np.abs(g_tg - g_ray)) < 1e-6


def test_taylor_goldstein_miles_howard_cutoff():
    """The rigorous Miles-Howard theorem (Ri>=1/4 everywhere => stable) is a
    hard mathematical guarantee, not a fitted benchmark: for the Hazel
    (1972) profile the minimum local Ri equals J, so J>=0.25 must give zero
    growth at every tested k, while J<0.25 must be genuinely unstable."""
    n = 400
    z = np.linspace(-10, 10, n + 2)[1:-1]
    k_values = np.linspace(0.05, 1.0, 30)

    for J in [0.0, 0.1, 0.2]:
        U, N2 = instability.hazel_profile(z, J)
        growth = instability.taylor_goldstein_growth_rate_curve(z, U, N2, k_values)
        assert growth.max() > 0.03, f"expected clear instability at J={J}"

    for J in [0.5, 1.0, 2.0]:
        U, N2 = instability.hazel_profile(z, J)
        growth = instability.taylor_goldstein_growth_rate_curve(z, U, N2, k_values)
        assert growth.max() == 0.0, f"Miles-Howard requires zero growth at J={J}"


def test_taylor_goldstein_j0_peak_matches_michalke():
    """At J=0 (no stratification) the Hazel profile IS the plain tanh shear
    layer, so the peak should match the same Michalke (1964) benchmark
    k*delta~=0.4446 used for the barotropic Rayleigh solver."""
    n = 400
    z = np.linspace(-10, 10, n + 2)[1:-1]
    U, N2 = instability.hazel_profile(z, 0.0)
    k_values = np.linspace(0.05, 1.0, 60)
    growth = instability.taylor_goldstein_growth_rate_curve(z, U, N2, k_values)
    k_peak = k_values[np.argmax(growth)]
    assert abs(k_peak - 0.4446) < 0.03


def test_okubo_weiss_pure_strain_exact():
    """With Gamma=0 (no vortex) the flow is exactly linear strain, so
    central differences are exact regardless of step size: W=4*alpha^2
    everywhere, strain-dominated (W>0)."""
    alpha = 0.5
    W = kinematics.okubo_weiss(3.0, -2.0, alpha=alpha, Gamma=0.0, sigma=1.0)
    assert abs(W - 4 * alpha ** 2) < 1e-8


def test_okubo_weiss_vortex_core_vorticity_dominated():
    """At the center of a pure (alpha=0) Lamb-Oseen-like vortex the flow is
    locally solid-body rotation: zero strain, vorticity=Gamma/(2*pi*sigma^2),
    so W=-zeta^2<0 (vorticity-dominated core)."""
    Gamma, sigma = 2.0, 1.0
    zeta_center = Gamma / (2 * np.pi * sigma ** 2)
    W = kinematics.okubo_weiss(0.0, 0.0, alpha=0.0, Gamma=Gamma, sigma=sigma)
    assert W < 0.0
    assert abs(W - (-zeta_center ** 2)) / zeta_center ** 2 < 0.05


def test_inertial_rhs_returns_to_start_after_one_period():
    """Numerically integrating du/dt=f*v, dv/dt=-f*u for exactly one period
    T=2*pi/f must return to the starting velocity."""
    f = 2.0
    dt = 1e-4
    nsteps = int(round((2 * np.pi / f) / dt))
    state = np.array([1.0, 0.0])
    for _ in range(nsteps):
        state = timestep.rk4(lambda t, s: rotation.inertial_rhs(s, f), state, dt)
    assert np.max(np.abs(state - [1.0, 0.0])) < 1e-4


def test_inertial_trajectory_constant_radius_and_speed():
    """The analytic solution is an exact circle: constant distance from its
    center, and constant speed (rotation conserves kinetic energy)."""
    f, u0, v0, x0, y0 = 1.3, 0.7, -0.4, 2.0, -1.0
    t = np.linspace(0, 7.3, 500)
    x, y, u, v = rotation.inertial_trajectory(t, u0, v0, f, x0, y0)
    w0 = u0 + 1j * v0
    center = (x0 + 1j * y0) + w0 / (1j * f)
    r = np.sqrt((x - center.real) ** 2 + (y - center.imag) ** 2)
    assert np.max(np.abs(r - abs(w0) / f)) < 1e-10
    assert np.max(np.abs(np.sqrt(u ** 2 + v ** 2) - abs(w0))) < 1e-10


def test_inertial_numeric_matches_analytic_trajectory():
    """Independent check: RK4-stepping the full (x,y,u,v) state must agree
    with the closed-form inertial_trajectory."""
    f, u0, v0 = 2.1, 1.0, 0.3
    dt, nsteps = 1e-3, 2000

    def rhs(t, s):
        x, y, u, v = s
        return np.array([u, v, f * v, -f * u])

    state = np.array([0.0, 0.0, u0, v0])
    for _ in range(nsteps):
        state = timestep.rk4(rhs, state, dt)
    x_a, y_a, u_a, v_a = rotation.inertial_trajectory(dt * nsteps, u0, v0, f)
    assert np.max(np.abs(state - [x_a, y_a, u_a, v_a])) < 1e-6


def test_rossby_and_burger_numbers():
    assert abs(scaling.rossby_number(U=1.0, f=1e-4, L=1e5) - 0.1) < 1e-12
    N, H, f, L = 0.01, 1000.0, 1e-4, 5e4
    Lr = scaling.deformation_radius(N, H, f)
    assert abs(Lr - N * H / f) < 1e-12
    assert abs(scaling.burger_number(N, H, f, L) - (Lr / L) ** 2) < 1e-10


def test_coriolis_parameter_known_values():
    Omega = 7.292e-5
    assert abs(scaling.coriolis_parameter(0.0)) < 1e-15
    assert abs(scaling.coriolis_parameter(30.0) - Omega) < 1e-10
    assert abs(scaling.coriolis_parameter(90.0) - 2 * Omega) < 1e-10


def test_classify_regime_known_points():
    assert scaling.classify_regime(2.0, 5.0) == 0
    assert scaling.classify_regime(0.1, 2.0) == 1
    assert scaling.classify_regime(0.1, 0.5) == 2
