# `gfdlib` — a guided tour of the shared library

`gfdlib` is the numeric core of the book: pure NumPy, Pyodide-safe (no
SciPy, no compiled extensions), one vetted implementation of each primitive.
Every notebook imports from it; none re-implements it. Every public function
has a docstring and a test in `tests/test_gfdlib.py`.

This page groups the modules by role and says which chapter uses each.
Signatures below are abbreviated — the docstrings are the reference.

## Cross-cutting infrastructure (used by nearly every chapter)

### `spectral` — pseudo-spectral grid
`Grid(n, L=2π)`: a doubly-periodic `[0,L)²` grid with cached real-FFT
wavenumber arrays, built **once per notebook** and reused every step.
Methods: `fft`/`ifft` (rfft2 wrappers; spectral fields have shape
`(n, n//2+1)`), `laplacian`, `invert_laplacian` (solve ∇²ψ = F, mean mode
zeroed), `ddx`, `ddy`, and a dealiased (2/3-rule) `jacobian(A, B)` =
AₓB_y − A_yBₓ.

### `timestep` — time integrators
- `rk4(rhs, y, dt, t)` — one classical RK4 step of dy/dt = rhs(t, y); works
  on flat state vectors or spectral fields.
- `ifrk4_step(F, rhs_nl, dt, L, t)` — integrating-factor RK4 for
  dF/dt = L·F + N(F) with diagonal linear L (e.g. hyperviscosity + exact
  Rossby propagation). The linear part is integrated exactly, so it is
  unconditionally stable; RK4's stability region covers the imaginary axis,
  so inviscid advection is stable at CFL-limited dt.

### `diagnostics` — budgets and spectra
`energy_enstrophy(psi_hat, grid)` and `isotropic_spectrum(F, grid)` — the
conservation checks and spectral diagnostics shown in the turbulence
chapters (18–19) and usable anywhere.

### `plotting` — the figure policy
`new_fig`, `field(ax, f, grid, signed=...)`, plus the shared colormaps:
diverging + symmetric limits for **signed** fields (vorticity, PV, height
anomaly), perceptually-uniform sequential for non-negative fields. Using
this module is what keeps every chapter's figures visually consistent.

## Part I — Foundations

### `kinematics` (ch. 1)
Prescribed 2D flow: `velocity_field` (steady strain + a Lamb-Oseen-like
vortex, closed form), `vortex_velocity`, and `okubo_weiss` (the
strain-vs-vorticity parameter $W$, from central differences of
`velocity_field` — same safety pattern as `rossby`'s ray equations).

### `rotation` (ch. 3)
Inertial oscillations: `inertial_rhs` (for RK4 integration) and
`inertial_trajectory` (closed-form circular solution, radius $|w_0|/f$,
period $2\pi/f$). Taylor-Proudman itself reuses `balance.thermal_wind_u`
directly (it is the $b\equiv0$ limit of thermal wind — no new solver).

### `scaling` (ch. 2, 4)
Dimensionless-number calculator: `rossby_number`, `ekman_number`,
`froude_number`, `reynolds_number`, `burger_number`, `deformation_radius`,
`coriolis_parameter`, and `classify_regime(Ro, Bu)` (the three-region
unbalanced / QG-baroclinic / QG-barotropic map ch. 4 renders).

## Part II — Balance and adjustment

### `balance` (ch. 5)
Closed-form thermal-wind machinery: `buoyancy_field` (idealized front),
`dbdy_frontal` (its exact meridional gradient), `thermal_wind_u`
(analytic integral of f u_z = −b_y). Purely diagnostic — no time stepping.

### `shallowwater` (ch. 6)
Linearized rotating shallow water on a β-plane: `coriolis(grid, beta, y0)`,
`rhs(state, grid, f)` for the (u, v, η) system, `potential_vorticity`, and
`invert_pv` — the geostrophic-adjustment and PV-conservation engine.

### `pv` (ch. 7–8)
Field builders for PV-inversion/advection demos: `gaussian_blob`,
`staircase_pv`.

## Part III — Quasi-geostrophy

### `qg` (ch. 8)
`dispersion_omega(kx, ky, beta)` — the barotropic Rossby-wave dispersion
relation ω = −βkₓ/k²; the ch. 8 notebook steps QGPV with `spectral.Grid` +
`timestep.ifrk4_step` directly.

### `rossby` (ch. 9)
WKB ray tracing on the sphere: `stationary_wavenumber2`, `dispersion_omega`,
`ray_rhs` (dλ/dt, dφ/dt, dm/dt), `launch_state`, `great_circle_deviation` —
reproduces the classic "great-circle" stationary-Rossby-wave ray paths.

## Part IV — Stratified flow and waves

### `internalwaves` (ch. 12)
`dispersion_omega(kx, kz, N, f)` (anisotropic dispersion ω² = (N²kₓ² +
f²k_z²)/k²), `beam_angle`, and `step_leapfrog` for the forced
internal-wave-beam simulation.

## Part V — Instabilities

### `convection` (ch. 14, reused in ch. 21)
Rayleigh–Bénard convection and its Lorenz-63 truncation: `neutral_ra`,
`growth_rate`, `ChannelGrid` (wall-bounded channel), `rhs_boussinesq`
(vorticity–streamfunction Boussinesq), `lorenz_rhs`.

### `instability` (ch. 15, 17)
Shear instability eigenproblems, solved with matrix methods (pure
`numpy.linalg`): Rayleigh/Kuo barotropic `growth_rate(_curve)`, stratified
Taylor–Goldstein `taylor_goldstein_growth_rate(_curve)`, `hazel_profile`,
and `double_shear_layer` initial conditions for nonlinear roll-up.

### `baroclinic` (ch. 16)
The book's centerpiece: `eady_growth_rate(_curve)` (the Eady dispersion
relation) and the two-layer QG model — `invert_2layer` (coupled PV
inversion) and `rhs_2layer` (with β and bottom friction r).

### `symmetric` (ch. 17)
The Ri–Ro stability map: `ertel_pv_ratio`, `critical_ri`, `classify`
(stable / symmetric / inertial regions of a balanced front).

## Part VI — Turbulence and circulation

### `mixing` (ch. 19)
Passive tracer stirred by 2D turbulence: `rhs_coupled` (vorticity + tracer
with mean gradient Γ) and `effective_diffusivity` (K_eff from the eddy flux
⟨v′c′⟩/Γ).

### `circulation` (ch. 20)
Wind-driven gyres: `stommel_analytic` / `solve_stommel_fd` and
`munk_analytic` / `solve_munk_fd` (analytic solutions plus 2nd-order
finite-difference solvers, verified against each other in the tests), with
`wind_stress_curl` forcing.

### `overturning` (ch. 21)
Buoyancy-driven overturning: `abyssal_profile` (Munk's 1966
upwelling–diffusion balance, exact solution), `surface_heating`
(differential surface-concentrated forcing), `rhs_overturning`
(horizontal convection, reusing ch. 14's `ChannelGrid`).

## Adding a primitive

1. Implement it in the right module (or a new `gfdlib/<topic>.py`) — pure
   NumPy, docstring stating conventions and citing `NOTATION.md` symbols.
2. Add a test in `tests/test_gfdlib.py`: an analytic limit, a conservation
   law, or a convergence order — not just "it runs".
3. If it's a new module: add it to the import list and `__all__` in
   `gfdlib/__init__.py`.
4. `make test` must pass; CI re-checks it plus a WASM export of every
   notebook.
