import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 14 — Rayleigh-Benard Convection & the Lorenz-63 Model

        **Physical question.** Heat a fluid layer from below. Below a critical
        temperature contrast it stays at rest, conducting heat quietly; above
        it, convection rolls appear; push harder still and the flow turns
        chaotic. After this chapter you should be able to predict the onset
        Rayleigh number from first principles, and explain how Lorenz reduced
        this whole PDE to three ODEs that already contain the essence of
        deterministic chaos.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### Setup and the two dimensionless numbers

        Boussinesq convection between free-slip, fixed-temperature plates at
        $z=0,1$, periodic in $x$, nondimensionalized by the depth $H$, the
        thermal diffusion time $H^2/\kappa$, and the imposed contrast
        $\Delta T$ (symbols as in [NOTATION](../notation)):

        $$\zeta_t + J(\psi,\zeta) = Pr\,\nabla^2\zeta + Pr\,Ra\,\theta_x,
          \qquad
          \theta_t + J(\psi,\theta) = \nabla^2\theta + \psi_x,
          \qquad \nabla^2\psi=\zeta,$$

        with $u=-\psi_z,\ w=\psi_x$ and $\theta$ the temperature *deviation*
        from the linear conduction profile (so $\theta\equiv0$ is the
        no-motion state, and the $+\psi_x$ term is advection of that
        background profile by $w$). Two numbers control everything:

        > **Rayleigh number**
        > $\;Ra=\dfrac{g\alpha\Delta T\,H^3}{\nu\kappa}$ — *buoyancy forcing
        > versus the two diffusivities that resist it.* The numerator is how
        > hard hot fluid is pushed up; the denominator is how fast viscosity
        > kills the motion and conduction erases the temperature anomaly
        > driving it. Convection requires the push to win: $Ra>Ra_c$.
        >
        > **Prandtl number** $\;Pr=\dfrac{\nu}{\kappa}$ — *which diffuses
        > faster, momentum or heat.* Air: $Pr\approx0.7$; water:
        > $Pr\approx7$; Earth's mantle: $Pr\sim10^{23}$ (momentum diffuses
        > essentially instantly; the flow is slaved to the temperature).

        ### Linear onset: why there is a critical $Ra$, and why a preferred cell size

        Linearize about rest ($\zeta,\theta$ small; drop the Jacobians) and
        insert normal modes $\sim e^{\sigma t}\sin(kx)\sin(\pi z)$ — the
        $\sin(\pi z)$ satisfies both free-slip ($\zeta=0$) and
        fixed-temperature ($\theta=0$) conditions at the plates. Setting the
        growth rate $\sigma=0$ gives the neutral curve:

        $$Ra_c(k)=\frac{(k^2+\pi^2)^3}{k^2},\qquad
          Ra_c=\min_k Ra_c(k)=\frac{27\pi^4}{4}\approx657.5\ \text{at}\ k_c=\frac{\pi}{\sqrt2}.$$

        The *shape* of $Ra_c(k)$ is worth a minute of thought. Very wide
        cells ($k\to0$) are inefficient — fluid must travel a long
        horizontal path for each unit of vertical heat transport — so the
        required forcing diverges as $1/k^2$. Very narrow cells
        ($k\to\infty$) put hot and cold fluid so close together that
        diffusion short-circuits them; that end diverges as $k^4$. In
        between sits a preferred cell width, $k_c$: cells roughly
        $2\sqrt2$ as wide as the layer is deep. The fluid *chooses its own
        pattern scale* — nothing in the forcing (uniform heating!) picked it.

        ### The Lorenz (1963) truncation

        Keep only the gravest mode of $\psi$ and the two thermally-relevant
        modes of $\theta$ compatible with the boundary conditions:

        $$\psi\propto X(t)\sin(k_cx)\sin(\pi z),\quad
          \theta\propto Y(t)\cos(k_cx)\sin(\pi z)-Z(t)\sin(2\pi z),$$

        so $X$ is the roll circulation speed, $Y$ the temperature contrast
        between rising and sinking fluid, and $Z$ the distortion of the
        *mean* vertical temperature profile away from linear conduction.
        Substituting and discarding every harmonic these three modes
        generate gives

        $$\dot X=\sigma(Y-X),\qquad \dot Y=rX-Y-XZ,\qquad \dot Z=XY-bZ,$$

        with $\sigma=Pr$, $r=Ra/Ra_c$ (forcing measured in units of onset),
        and $b=8/3$ (a geometric factor from $k_c$). The nonlinearities
        $XZ$ and $XY$ are the two survivors of the Jacobians — advection of
        the mean profile by the roll, and construction of the mean-profile
        distortion by the heat flux. Fixed-point anatomy:

        - **Origin** $(0,0,0)$ = conduction; stable for $r<1$, loses
          stability at exactly $r=1$ — the truncation *inherits* the
          $Ra_c$ of the full PDE by construction.
        - **Pair** $(\pm\sqrt{b(r-1)},\pm\sqrt{b(r-1)},r-1)$ = steady rolls
          (the sign is the roll's rotation direction); stable for
          $1<r<r_H\approx24.74$ (at $\sigma=10$, $b=8/3$).
        - Beyond $r_H$, **all three** fixed points are unstable — via a
          *subcritical* Hopf bifurcation, so there is no stable orbit to
          land on — and trajectories wander forever on the chaotic
          **Lorenz attractor**, hopping irregularly between the two rolls'
          basins.
        """
    )
    return


@app.cell
async def _(mo):
    # --- vetted primitives from gfdlib (see CLAUDE.md) -------------------
    import sys

    if sys.platform == "emscripten":
        import micropip
        await micropip.install(
            str(mo.notebook_location() / "public" / "gfdlib-0.1.0-py3-none-any.whl")
        )
    else:
        sys.path.insert(0, str(mo.notebook_dir().parent))  # repo root

    import numpy as np
    import matplotlib.pyplot as plt
    from gfdlib import timestep, plotting, convection
    return convection, np, plotting, plt, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        **2D convection:** periodic in $x$ (FFT), Dirichlet in $z$ (finite
        difference, interior points only) — `gfdlib.convection.ChannelGrid`.
        Each Poisson solve $\nabla^2\psi=\zeta$ uses a Thomas (tridiagonal)
        algorithm per horizontal wavenumber, vectorized across all
        wavenumbers, with the forward-elimination coefficients precomputed
        once. Time-stepping is classical RK4 at the (restrictive) diffusive
        CFL $\Delta t\approx0.2\min(\Delta x,\Delta z)^2$
        (`gfdlib.convection.rhs_boussinesq`).

        **Lorenz-63:** `gfdlib.convection.lorenz_rhs` with the same
        `gfdlib.timestep.rk4` used everywhere else in this book — no special
        machinery needed for a 3-variable ODE.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part A — Onset and the full 2D simulation

        The curve below is the marginal-stability boundary $Ra_c(k)$: below
        it the conduction state is linearly stable at that wavenumber; above
        it, unstable. The domain width $L_x=2\sqrt2$ is chosen so the first
        periodic Fourier mode sits exactly at $k_c$.
        """
    )
    return


@app.cell(hide_code=True)
def _(convection, np, plt):
    # --- diagnostic: neutral stability curve, no simulation needed ---------
    _k = np.linspace(0.3, 5, 400)
    fig0, ax0 = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax0.plot(_k, convection.neutral_ra(_k), lw=2)
    ax0.axhline(convection.RA_C, color="crimson", ls="--", lw=1.2,
                label=f"$Ra_c={convection.RA_C:.1f}$")
    ax0.axvline(convection.K_C, color="gray", ls=":", lw=1)
    ax0.set_ylim(0, 3000)
    ax0.set_xlabel("wavenumber $k$"); ax0.set_ylabel("$Ra_c(k)$")
    ax0.legend(fontsize=9); ax0.grid(alpha=0.3)
    ax0.set_title("neutral stability curve")
    fig0
    return


@app.cell(hide_code=True)
def _(mo):
    rb_Ra = mo.ui.slider(300, 5000, step=100, value=2000,
                         label="$Ra$", show_value=True)
    rb_Pr = mo.ui.slider(0.1, 10.0, step=0.1, value=1.0,
                         label="$Pr$", show_value=True)
    rb_nx = mo.ui.dropdown(options={"32": 32, "64": 64, "128": 128}, value="64",
                           label="$n_x$")
    rb_nz = mo.ui.dropdown(options={"8": 8, "16": 16, "24": 24, "32": 32}, value="16",
                           label="$n_z$")
    rb_T = mo.ui.slider(0.2, 2.0, step=0.1, value=1.0,
                        label="run time $T$", show_value=True)
    mo.vstack([
        mo.hstack([rb_Ra, rb_Pr], justify="start"),
        mo.hstack([rb_nx, rb_nz, rb_T], justify="start"),
        mo.md(r"> $Ra<Ra_c\approx658$: conduction persists. $Ra>Ra_c$: rolls "
              r"grow and saturate. In the browser, start with $n_x=32,n_z=8$."),
    ])
    return rb_Pr, rb_Ra, rb_T, rb_nx, rb_nz


@app.cell(hide_code=True)
def _(mo):
    get_rb, set_rb = mo.state(None)
    return get_rb, set_rb


@app.cell(hide_code=True)
def _(mo):
    rb_btn = mo.ui.run_button(label="▶ Run 2D simulation")
    rb_btn
    return (rb_btn,)


@app.cell
def _(
    convection, mo, np, rb_Pr, rb_Ra, rb_T, rb_btn, rb_nx, rb_nz, set_rb,
    timestep,
):
    mo.stop(not rb_btn.value, mo.callout(
        mo.md("Press **▶ Run 2D simulation** above to start."), kind="info"
    ))

    _nx, _nz, _Ra, _Pr, _T = rb_nx.value, rb_nz.value, rb_Ra.value, rb_Pr.value, rb_T.value
    _Lx = 2 * np.sqrt(2)
    grid = convection.ChannelGrid(_nx, _nz, _Lx)

    _rng = np.random.default_rng(7)
    _theta = 1e-3 * _rng.standard_normal((_nx, _nz))
    _zeta = np.zeros((_nx, _nz))
    _state = np.stack([_zeta, _theta])

    _dt = 0.2 * min(grid.dx, grid.dz) ** 2
    _nsteps = max(20, int(_T / _dt))
    _nsub = max(1, _nsteps // 60)

    def _rhs(t, s):
        return convection.rhs_boussinesq(s, grid, _Ra, _Pr)

    _snaps_th, _snaps_psi, _t_hist = [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _snaps_th.append(_state[1].copy())
            _snaps_psi.append(grid.poisson_solve(_state[0]).copy())
            _t_hist.append(_t)
        if _s < _nsteps:
            _state = timestep.rk4(_rhs, _state, _dt)
            _t += _dt

    set_rb({
        "snaps_th": _snaps_th, "snaps_psi": _snaps_psi,
        "t": np.array(_t_hist), "grid": grid, "Ra": _Ra, "Pr": _Pr,
    })

    mo.callout(mo.md(
        f"**Done.** $Ra={_Ra}$, $Pr={_Pr}$, grid {_nx}×{_nz}, {_nsteps} steps "
        f"→ {len(_snaps_th)} frames. $Ra/Ra_c={_Ra/convection.RA_C:.2f}$."
    ), kind="success")
    return (grid,)


@app.cell
def _(get_rb, mo):
    rb_sim = get_rb()
    mo.stop(rb_sim is None, mo.callout(
        mo.md("2D results will appear here once the simulation has run."), kind="warn"
    ))
    return (rb_sim,)


@app.cell(hide_code=True)
def _(mo, rb_sim):
    rb_frame = mo.ui.slider(0, len(rb_sim["snaps_th"]) - 1, step=1,
                            value=len(rb_sim["snaps_th"]) - 1,
                            label="frame", show_value=True)
    return (rb_frame,)


@app.cell(hide_code=True)
def _(mo, np, plt, rb_frame, rb_sim):
    # --- diagnostic: snapshot browser, temperature + streamfunction --------
    _i = rb_frame.value
    _g = rb_sim["grid"]
    _th = rb_sim["snaps_th"][_i]; _psi = rb_sim["snaps_psi"][_i]
    _XX, _ZZ = np.meshgrid(_g.x, _g.z, indexing="ij")

    fig1, axs1 = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    _vm = max(float(np.abs(_th).max()), 1e-10)
    _im0 = axs1[0].pcolormesh(_XX, _ZZ, _th, cmap="RdBu_r", vmin=-_vm, vmax=_vm, shading="auto")
    plt.colorbar(_im0, ax=axs1[0], fraction=0.046, label=r"$\theta$")
    axs1[0].set_title(f"temperature  $t={rb_sim['t'][_i]:.3f}$")
    axs1[0].set_xlabel("$x$"); axs1[0].set_ylabel("$z$")

    _vm2 = max(float(np.abs(_psi).max()), 1e-10)
    _im1 = axs1[1].pcolormesh(_XX, _ZZ, _psi, cmap="PuOr", vmin=-_vm2, vmax=_vm2, shading="auto")
    plt.colorbar(_im1, ax=axs1[1], fraction=0.046, label=r"$\psi$")
    axs1[1].set_title(f"streamfunction  $t={rb_sim['t'][_i]:.3f}$")
    axs1[1].set_xlabel("$x$"); axs1[1].set_ylabel("$z$")
    mo.vstack([rb_frame, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part B — The Lorenz-63 truncation

        Try $\sigma=10,\ r=28,\ b=8/3$ for the classic butterfly attractor.
        The same $(\sigma,r,b)$ chosen here also drive the Lyapunov-separation
        and bifurcation diagnostics below.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    lz_sigma = mo.ui.slider(1.0, 20.0, step=0.5, value=10.0, label=r"$\sigma$", show_value=True)
    lz_r = mo.ui.slider(0.5, 60.0, step=0.5, value=28.0, label="$r$", show_value=True)
    lz_b = mo.ui.slider(0.5, 4.0, step=1 / 6, value=8 / 3, label="$b$", show_value=True)
    lz_T = mo.ui.slider(10, 100, step=5, value=50, label="integration time", show_value=True)
    lz_ntraj = mo.ui.slider(1, 5, step=1, value=3, label="trajectories", show_value=True)
    mo.vstack([
        mo.hstack([lz_sigma, lz_r, lz_b], justify="start"),
        mo.hstack([lz_T, lz_ntraj], justify="start"),
    ])
    return lz_T, lz_b, lz_ntraj, lz_r, lz_sigma


@app.cell(hide_code=True)
def _(mo):
    get_lz, set_lz = mo.state(None)
    return get_lz, set_lz


@app.cell(hide_code=True)
def _(mo):
    lz_btn = mo.ui.run_button(label="▶ Run Lorenz")
    lz_btn
    return (lz_btn,)


@app.cell
def _(
    convection, lz_T, lz_b, lz_btn, lz_ntraj, lz_r, lz_sigma, mo, np, set_lz,
    timestep,
):
    mo.stop(not lz_btn.value, mo.callout(
        mo.md("Press **▶ Run Lorenz** above to integrate."), kind="info"
    ))

    _sigma, _r, _b, _T, _n = lz_sigma.value, lz_r.value, lz_b.value, lz_T.value, lz_ntraj.value
    _dt = 0.01
    _nsteps = int(_T / _dt)
    _rng = np.random.default_rng(42)

    def _rhs(t, s):
        return convection.lorenz_rhs(s, _sigma, _r, _b)

    _trajs = []
    for _i in range(_n):
        _ic = np.array([1.0, 1.0, 1.0]) + _rng.normal(0, 0.5, 3)
        _traj = np.zeros((_nsteps + 1, 3))
        _traj[0] = _ic
        for _s in range(_nsteps):
            _traj[_s + 1] = timestep.rk4(_rhs, _traj[_s], _dt)
        _trajs.append(_traj)

    set_lz({
        "trajectories": _trajs, "t": np.linspace(0, _T, _nsteps + 1),
        "sigma": _sigma, "r": _r, "b": _b,
    })

    mo.callout(mo.md(
        f"**Done.** {_n} trajectory/ies, $t\\in[0,{_T}]$, "
        f"$\\sigma={_sigma}$, $r={_r}$, $b={_b:.2f}$."
    ), kind="success")
    return


@app.cell
def _(get_lz, mo):
    lz_sim = get_lz()
    mo.stop(lz_sim is None, mo.callout(
        mo.md("Lorenz results will appear here after running."), kind="warn"
    ))
    return (lz_sim,)


@app.cell(hide_code=True)
def _(lz_sim, plt):
    # --- diagnostic: phase portraits + time series --------------------------
    _colors = ["#2563eb", "crimson", "#16a34a", "#d97706", "#7c3aed"]
    fig2, axs2 = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for _idx, _traj in enumerate(lz_sim["trajectories"]):
        _c = _colors[_idx % len(_colors)]
        _skip = max(1, len(lz_sim["t"]) // 2000)
        axs2[0, 0].plot(_traj[::_skip, 0], _traj[::_skip, 1], lw=0.4, alpha=0.7, color=_c)
        axs2[0, 1].plot(_traj[::_skip, 0], _traj[::_skip, 2], lw=0.4, alpha=0.7, color=_c)
        axs2[0, 2].plot(_traj[::_skip, 1], _traj[::_skip, 2], lw=0.4, alpha=0.7, color=_c)
        if _idx == 0:
            axs2[1, 0].plot(lz_sim["t"], _traj[:, 0], lw=0.7, color=_c)
            axs2[1, 1].plot(lz_sim["t"], _traj[:, 1], lw=0.7, color=_c)
            axs2[1, 2].plot(lz_sim["t"], _traj[:, 2], lw=0.7, color=_c)
    axs2[0, 0].set_xlabel("$X$"); axs2[0, 0].set_ylabel("$Y$")
    axs2[0, 1].set_xlabel("$X$"); axs2[0, 1].set_ylabel("$Z$")
    axs2[0, 2].set_xlabel("$Y$"); axs2[0, 2].set_ylabel("$Z$")
    axs2[1, 0].set_xlabel("$t$"); axs2[1, 0].set_ylabel("$X$")
    axs2[1, 1].set_xlabel("$t$"); axs2[1, 1].set_ylabel("$Y$")
    axs2[1, 2].set_xlabel("$t$"); axs2[1, 2].set_ylabel("$Z$")
    for _ax in axs2.ravel():
        _ax.grid(alpha=0.2)
    fig2.suptitle(f"$\\sigma={lz_sim['sigma']},\\ r={lz_sim['r']},\\ b={lz_sim['b']:.2f}$")
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Sensitivity to initial conditions

        Two trajectories released $10^{-8}$ apart, at the same
        $(\sigma,r,b)$ chosen above. Both are first spun up for 15 time units
        from a fixed start so they land **on the attractor** before the clock
        starts — the transient relaxation onto the attractor has its own
        (much smaller) rate and is not the Lyapunov exponent.

        **How to read the right panel:** on a log axis, exponential
        separation $\|\delta u\|\sim\varepsilon e^{\lambda_1t}$ is a straight
        line whose slope is the leading **Lyapunov exponent** $\lambda_1$.
        The line must eventually *saturate* — the attractor is a bounded
        object, so two trajectories can never get farther apart than its
        diameter. The consequence Lorenz drew is arithmetic, and brutal:
        improving the initial condition by a factor of 10 buys only
        $\ln(10)/\lambda_1\approx2.6$ extra time units of predictability.
        Accuracy bought exponentially, predictability gained linearly —
        that asymmetry is why weather forecasts have a horizon (~2 weeks)
        that no conceivable observing system will push past by much.
        """
    )
    return


@app.cell(hide_code=True)
def _(convection, lz_sim, np, plt, timestep):
    # --- diagnostic: Lyapunov separation ------------------------------------
    _sigma, _r, _b = lz_sim["sigma"], lz_sim["r"], lz_sim["b"]
    _dt, _T2 = 0.01, 20.0
    _nsteps = int(_T2 / _dt)
    _eps = 1e-8

    def _rhs(t, s):
        return convection.lorenz_rhs(s, _sigma, _r, _b)

    # spin up onto the attractor first -- starting the clock from (1,1,1)
    # measures the transient relaxation onto the attractor, not the
    # asymptotic divergence rate, and badly underestimates lambda_1
    _s0 = np.array([1.0, 1.0, 1.0])
    for _ in range(int(15.0 / _dt)):
        _s0 = timestep.rk4(_rhs, _s0, _dt)

    _s1 = np.zeros((_nsteps + 1, 3)); _s2 = np.zeros((_nsteps + 1, 3))
    _s1[0] = _s0; _s2[0] = _s0 + [_eps, 0, 0]
    for _i in range(_nsteps):
        _s1[_i + 1] = timestep.rk4(_rhs, _s1[_i], _dt)
        _s2[_i + 1] = timestep.rk4(_rhs, _s2[_i], _dt)
    _t = np.linspace(0, _T2, _nsteps + 1)
    _sep = np.linalg.norm(_s1 - _s2, axis=1)

    fig3, axs3 = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axs3[0].plot(_t, _s1[:, 0], lw=1, label="trajectory 1")
    axs3[0].plot(_t, _s2[:, 0], lw=1, ls="--", alpha=0.8, label="trajectory 1 + $10^{-8}$")
    axs3[0].set_xlabel("$t$"); axs3[0].set_ylabel("$X$"); axs3[0].legend(fontsize=9)
    axs3[0].grid(alpha=0.3)

    _valid = _sep > 1e-15
    axs3[1].semilogy(_t[_valid], _sep[_valid], lw=1.5, color="crimson", label=r"$\|\delta u(t)\|$")
    _mask = (_t > 1) & (_t < 8) & _valid
    if _mask.sum() > 5:
        _lam = np.polyfit(_t[_mask], np.log(_sep[_mask]), 1)[0]
        axs3[1].semilogy(_t[_mask], np.exp(np.log(_sep[_mask][0]) + _lam * (_t[_mask] - _t[_mask][0])),
                         "k--", lw=1.2, label=f"$\\lambda_1\\approx{_lam:.2f}$")
    axs3[1].set_xlabel("$t$"); axs3[1].set_ylabel(r"separation $\|\delta u\|$")
    axs3[1].legend(fontsize=9); axs3[1].grid(alpha=0.2)
    fig3
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Bifurcation diagram

        Sweep $r$ at fixed $\sigma,b$ (from the sliders above): for each $r$,
        integrate past the initial transient, then record every local maximum
        of $Z(t)$. A single curve means a stable fixed point or simple
        oscillation; a vertical smear means chaos. The dashed line is the
        analytic fixed-point value $Z=r-1$.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    bif_rmax = mo.ui.slider(20.0, 60.0, step=2.0, value=60.0,
                            label="$r_{max}$", show_value=True)
    bif_nr = mo.ui.slider(100, 500, step=50, value=300,
                          label="number of $r$ values", show_value=True)
    bif_btn = mo.ui.run_button(label="▶ Run bifurcation sweep")
    mo.vstack([mo.hstack([bif_rmax, bif_nr], justify="start"), bif_btn])
    return bif_btn, bif_nr, bif_rmax


@app.cell(hide_code=True)
def _(mo):
    get_bif, set_bif = mo.state(None)
    return get_bif, set_bif


@app.cell
def _(
    bif_btn, bif_nr, bif_rmax, convection, lz_sim, mo, np, set_bif, timestep,
):
    mo.stop(not bif_btn.value, mo.callout(
        mo.md("Press **▶ Run bifurcation sweep** above to start."), kind="info"
    ))

    _sigma, _b = lz_sim["sigma"], lz_sim["b"]
    _r_vals = np.linspace(0.5, bif_rmax.value, bif_nr.value)
    _dt = 0.01

    def _rhs(t, s):
        return convection.lorenz_rhs(s, _sigma, _r_vals, _b)

    # vectorized across all r values at once: state shape (3, n_r)
    _state = np.ones((3, len(_r_vals))) + 0.1 * np.arange(len(_r_vals)) / len(_r_vals)

    for _ in range(int(20.0 / _dt)):        # discard transient
        _state = timestep.rk4(_rhs, _state, _dt)

    _maxima_r, _maxima_z = [], []
    _z_prev2 = _z_prev = None
    for _ in range(int(40.0 / _dt)):        # collect local maxima of Z
        _state = timestep.rk4(_rhs, _state, _dt)
        _z = _state[2]
        if _z_prev2 is not None:
            _is_max = (_z_prev > _z_prev2) & (_z_prev > _z)
            _maxima_r.append(_r_vals[_is_max])
            _maxima_z.append(_z_prev[_is_max])
        _z_prev2, _z_prev = _z_prev, _z

    set_bif({
        "r": _r_vals,
        "maxima_r": np.concatenate(_maxima_r) if _maxima_r else np.array([]),
        "maxima_z": np.concatenate(_maxima_z) if _maxima_z else np.array([]),
    })

    mo.callout(mo.md(f"**Done.** Swept {len(_r_vals)} values of $r$."), kind="success")
    return


@app.cell(hide_code=True)
def _(get_bif, mo, np, plt):
    _bif = get_bif()
    mo.stop(_bif is None, mo.callout(
        mo.md("The bifurcation diagram will appear here once the sweep has run."),
        kind="warn",
    ))
    fig4, ax4 = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax4.plot(_bif["maxima_r"], _bif["maxima_z"], ",", color="#2563eb", alpha=0.5)
    _r_line = np.linspace(1, _bif["r"].max(), 200)
    ax4.plot(_r_line, _r_line - 1, "k--", lw=1, label="$Z=r-1$ (fixed point)")
    ax4.axvline(24.74, color="crimson", ls=":", lw=1, label=r"$r_H\approx24.74$")
    ax4.set_xlabel("$r$"); ax4.set_ylabel("local maxima of $Z$")
    ax4.legend(fontsize=9); ax4.grid(alpha=0.3)
    fig4
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Find $Ra_c$ numerically.** In Part A, sweep $Ra$ near 658 (try
          600, 660, 750): watch $\theta$ decay, barely persist, or grow. The
          growth rate right at onset is *zero*, so near-critical runs take a
          long time to declare themselves — that slowness ("critical slowing
          down") is itself a universal signature of being near a
          bifurcation, used today to detect approaching tipping points in
          climate records.
        - **Prandtl number.** At fixed $Ra=2000$, compare $Pr=0.1$ (air-like
          order of magnitude) and $Pr=10$ (water-like): which develops
          convection faster? Note $Ra$ doesn't contain the growth *rate* —
          $Pr$ multiplies the buoyancy torque in the $\zeta$ equation, so it
          sets how quickly the roll spins up even though it cannot change
          *whether* it does.
        - **Read the bifurcation diagram.** Identify the single-branch region
          (steady rolls), locate $r_H$, and find at least one narrow window
          of $r>24.74$ where the smear briefly collapses back to a simple
          curve (a periodic window inside the chaos — chaos and order
          interleave at every scale of $r$).
        - **Doubling time.** From the fitted $\lambda_1$, compute how much
          longer you could predict the flow if your initial condition were
          10× more precise ($\Delta t=\ln 10/\lambda_1$). Now 100× more
          precise. Notice the returns diminishing.

        ### What you should have seen

        Below $Ra_c\approx658$ the 2D simulation's temperature perturbation
        decays to zero; above it, rolls grow and saturate into steady
        convection — exactly the transition the neutral curve predicts. The
        Lorenz truncation reproduces the same story in three variables: a
        stable origin for $r<1$, stable convective fixed points for
        $1<r<r_H\approx24.74$, and a chaotic attractor beyond, where nearby
        trajectories separate exponentially at rate $\lambda_1\approx0.9$.
        The bifurcation diagram makes the whole route from order to chaos
        visible in a single plot — this is the calculation, drawn from a real
        fluid-dynamics problem, that founded chaos theory.

        **Where this goes next.** This chapter is the book's template for
        *every* instability: identify a control parameter, find the
        critical value by linear theory, then watch nonlinearity decide
        what grows into. Chs. 15–17 replay exactly that script with shear
        (Rayleigh–Kuo), stratified shear (Taylor–Goldstein), and rotation
        (symmetric instability) as the antagonists; Ch. 16's baroclinic
        instability is the version that makes weather. And the machinery
        built here is reused literally: Ch. 21 drives this same
        `ChannelGrid` Boussinesq solver with *differential* surface heating
        to produce the ocean's overturning circulation. The predictability
        lesson, meanwhile, is the founding theorem of ensemble weather
        forecasting.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
