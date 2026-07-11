import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 6 — Rotating Shallow Water & Geostrophic Adjustment

        **Physical question.** Release a fluid at rest with a bump of height
        on it. Fast gravity waves radiate the disturbance outward, but
        rotation traps a balanced remnant behind. After this chapter you
        should be able to predict, from the ratio of the bump's size to the
        deformation radius, how much of the initial disturbance survives as
        organized geostrophic flow — and, on a $\beta$-plane, which way that
        remnant drifts.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### From shallow water to the linear system

        The rotating shallow-water (RSW) equations govern a thin layer of
        constant-density fluid with a free surface at $H+\eta$ (symbols as
        in [NOTATION](../notation)):

        $$\frac{D\mathbf u}{Dt}+f\hat{\mathbf z}\times\mathbf u=-g\nabla\eta,
          \qquad \eta_t+\nabla\cdot\big[(H+\eta)\mathbf u\big]=0.$$

        For small disturbances about rest ($|\mathbf u|$ small,
        $|\eta|\ll H$) drop every product of two small quantities: the
        advection $(\mathbf u\cdot\nabla)\mathbf u$ and the $\eta\mathbf u$
        flux both go, leaving a *linear* system. Now nondimensionalize —
        time by $f_0^{-1}$ (the rotation period is the clock), height by
        $H$, velocity by the gravity-wave speed $c=\sqrt{gH}$, and length by
        the one scale the system builds from its own constants:

        > **Rossby deformation radius**
        > $\;L_R=\dfrac{c}{f_0}=\dfrac{\sqrt{gH}}{f_0}$ — *the distance a
        > gravity wave travels in one rotation period; the scale at which
        > rotation and buoyancy contest control of the flow.* Atmosphere:
        > $L_R\sim1000\,$km. Ocean (first baroclinic mode): $L_R\sim
        > 30$–$50\,$km — which is why ocean "weather" (eddies) is 20× smaller
        > than atmospheric weather. Disturbances **wider** than $L_R$ feel
        > rotation before they can disperse; **narrower** ones disperse
        > before rotation matters. That single sentence is this whole
        > chapter.

        The nondimensional system (hats), on a $\beta$-plane:

        $$\hat\eta_t+\hat u_x+\hat v_y=0,\qquad
          \hat u_t-\hat f\hat v=-\hat\eta_x,\qquad
          \hat v_t+\hat f\hat u=-\hat\eta_y,
          \qquad \hat f=1+\hat\beta\,(y-y_0).$$

        ### The wave families (substitute and see)

        Try $e^{i(kx+ly-\omega t)}$ on the $f$-plane ($\hat\beta=0$). The
        $3\times3$ linear system has solutions only when its determinant
        vanishes, which factors into $\omega\big(\omega^2-(1+k^2+l^2)\big)=0$:

        - $\omega=0$: a **steady geostrophic mode** — balance is an *exact
          solution* of the linear equations, at every wavenumber. It doesn't
          oscillate; it just sits there.
        - $\omega=\pm\sqrt{1+k^2+l^2}$: fast **inertia–gravity waves**
          (dimensionally $\omega^2=f_0^2+c^2K^2$ — pure inertial
          oscillations at large scale, pure gravity waves at small scale).
          These are the radiators: they carry energy away.

        With $\hat\beta\neq0$ the $\omega=0$ mode is promoted to a slow,
        westward-propagating **Rossby wave**,
        $\omega_R=-\hat\beta k/(k^2+1)$ — the balanced remnant stops being
        exactly steady (Chs. 8–9 take that thread up properly).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Potential vorticity: the conserved skeleton

        Which combination of fields *doesn't* radiate away? Cross-differentiate
        the two momentum equations to get the vorticity equation, and combine
        it with the height equation ($\zeta=\hat v_x-\hat u_y$, $f$-plane):

        $$\zeta_t+(\hat u_x+\hat v_y)=0
          \quad\text{and}\quad
          \hat\eta_t+(\hat u_x+\hat v_y)=0
          \;\;\Longrightarrow\;\;
          \partial_t\underbrace{(\zeta-\hat\eta)}_{q}=0.$$

        The divergence — the one thing gravity waves are made of — cancels
        *identically*, so the **linear potential vorticity** $q=\zeta-\hat\eta$
        is conserved pointwise at every location, for all time, no matter how
        violent the wave transients. (On the $\beta$-plane it instead obeys
        $q_t=-\hat\beta\hat v$: meridional motion across the planetary
        vorticity gradient is the *only* thing that can change it.)

        This turns adjustment into a prediction machine. The end state is
        steady and geostrophic ($\hat u_{bal}=-\eta_{bal,y}$,
        $\hat v_{bal}=\eta_{bal,x}$, so $\zeta_{bal}=\nabla^2\eta_{bal}$),
        and it must carry the same $q$ the fluid started with:

        $$(\nabla^2-1)\,\eta_{bal}=q(t{=}0)=-\hat\eta_0
          \qquad\text{(since }\mathbf u=0\text{ initially).}$$

        The final balanced state is therefore fixed **before a single
        gravity wave has moved** — adjustment merely redistributes energy
        around a PV field it cannot touch. In Fourier space the inversion
        reads $\eta_{bal}=\hat\eta_0/(k^2+l^2+1)$: a **low-pass filter with
        cutoff at the deformation radius**. Scales $\gg L_R$ ($k\ll1$) pass
        through untouched; scales $\ll L_R$ are annihilated. Keep that
        filter in mind when you compare the presets below.
        """
    )
    return


@app.cell
async def _(mo):
    # --- vetted primitives from gfdlib (see CLAUDE.md) -------------------
    import sys

    if sys.platform == "emscripten":
        # running in the browser (Pyodide/WASM): install the gfdlib wheel
        # shipped alongside the exported notebook in its public/ folder
        import micropip
        await micropip.install(
            str(mo.notebook_location() / "public" / "gfdlib-0.1.0-py3-none-any.whl")
        )
    else:
        sys.path.insert(0, str(mo.notebook_dir().parent))  # repo root

    import numpy as np
    import matplotlib.pyplot as plt
    from gfdlib import spectral, timestep, plotting, shallowwater
    return np, plotting, plt, shallowwater, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Pseudo-spectral on a doubly-periodic $[0,L)^2$ domain, $L=20\,L_R$:
        derivatives are exact in Fourier space (`gfdlib.spectral.Grid`); the
        Coriolis terms $\hat f\hat u,\hat f\hat v$ are evaluated in physical
        space since $\hat f=\hat f(y)$. The system is exactly linear, so
        classical RK4 (`gfdlib.timestep.rk4`, whose stability region covers
        the imaginary axis) integrates the gravity- and Rossby-wave
        oscillations stably at a CFL-limited $\Delta t\approx0.4\,\Delta x$
        (`gfdlib.shallowwater.rhs`). There are no quadratic nonlinear terms,
        so no dealiasing is needed.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    scenario = mo.ui.dropdown(
        options={
            "Small-scale burst — waves win (σ ≪ L_R)": "small",
            "Large-scale surge — balance survives (σ ≫ L_R)": "large",
            "β-plane drift — Rossby propagation": "beta",
        },
        value="Small-scale burst — waves win (σ ≪ L_R)",
        label="Preset scenario",
    )
    scenario
    return (scenario,)


@app.cell(hide_code=True)
def _(mo, scenario):
    _presets = {
        "small": dict(n=128, sigma=0.5, amp=0.5, beta=0.0,  T=15.0),
        "large": dict(n=128, sigma=3.0, amp=0.5, beta=0.0,  T=15.0),
        "beta":  dict(n=128, sigma=2.0, amp=0.5, beta=0.16, T=15.0),
    }
    _p = _presets[scenario.value]

    n_ui = mo.ui.dropdown(
        options={"64": 64, "128": 128, "256": 256}, value=str(_p["n"]),
        label="resolution $n$",
    )
    sigma_ui = mo.ui.slider(0.25, 4.0, step=0.25, value=_p["sigma"],
                            label="perturbation width $\\sigma/L_R$", show_value=True)
    amp_ui = mo.ui.slider(0.1, 1.0, step=0.05, value=_p["amp"],
                          label="amplitude $\\eta_0/H$", show_value=True)
    beta_ui = mo.ui.slider(0.0, 0.5, step=0.02, value=_p["beta"],
                           label="$\\hat\\beta$", show_value=True)
    T_ui = mo.ui.slider(5.0, 30.0, step=2.5, value=_p["T"],
                        label="run time $T$ ($f_0^{-1}$)", show_value=True)

    mo.vstack([
        mo.hstack([n_ui, sigma_ui, amp_ui], justify="start"),
        mo.hstack([beta_ui, T_ui], justify="start"),
        mo.md(
            r"> **In the browser (WASM):** start with $n=64$. The domain is "
            r"$20\,L_R$ across — keep $T$ well under 20 so the radiating "
            r"gravity waves don't wrap around the periodic boundary."
        ),
    ])
    return amp_ui, beta_ui, n_ui, sigma_ui, T_ui


@app.cell(hide_code=True)
def _(mo):
    get_sim, set_sim = mo.state(None)   # results survive the run button resetting
    return get_sim, set_sim


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Run simulation")
    run_btn
    return (run_btn,)


@app.cell
def _(
    amp_ui, beta_ui, mo, n_ui, np, run_btn, set_sim, shallowwater, sigma_ui,
    spectral, T_ui, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run simulation** above to start."), kind="info"
    ))

    _n, _sigma, _amp, _beta, _T = (
        n_ui.value, sigma_ui.value, amp_ui.value, beta_ui.value, T_ui.value
    )
    _L = 20.0

    grid = spectral.Grid(_n, L=_L)               # build ONCE, reuse every step
    _f = shallowwater.coriolis(grid, _beta)

    _x0 = _y0 = 0.5 * _L
    _eta0 = _amp * np.exp(
        -((grid.x - _x0) ** 2 + (grid.y - _y0) ** 2) / (2 * _sigma ** 2)
    )
    _state = np.stack([_eta0, np.zeros_like(_eta0), np.zeros_like(_eta0)])
    _q0 = shallowwater.potential_vorticity(_state, grid)
    _eta_bal = shallowwater.invert_pv(_q0, grid)   # predicted end state, from t=0

    _dt = 0.4 * grid.dx
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 100)

    def _rhs(t, s):
        return shallowwater.rhs(s, grid, _f)

    _snaps, _t_hist, _KE, _PE = [], [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _snaps.append(_state.astype(np.float32))
            _t_hist.append(_t)
            _eta_s, _u_s, _v_s = _state
            _KE.append(0.5 * np.mean(_u_s ** 2 + _v_s ** 2))
            _PE.append(0.5 * np.mean(_eta_s ** 2))
        if _s < _nsteps:
            _state = timestep.rk4(_rhs, _state, _dt)
            _t += _dt

    _qf = shallowwater.potential_vorticity(_state, grid)

    set_sim({
        "snaps": _snaps, "t": np.array(_t_hist),
        "KE": np.array(_KE), "PE": np.array(_PE),
        "q0": _q0, "qf": _qf, "eta_bal": _eta_bal,
        "grid": grid, "beta": _beta, "sigma": _sigma,
    })

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_snaps)} frames."
    ), kind="success")
    return (grid,)


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the simulation has run."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(mo, sim):
    frame_ui = mo.ui.slider(0, len(sim["snaps"]) - 1, step=1, value=len(sim["snaps"]) - 1,
                            label="frame", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, mo, plotting, plt, sim):
    # --- diagnostic 1: snapshot browser -- height + velocity vectors -------
    _i = frame_ui.value
    _eta, _u, _v = sim["snaps"][_i]
    _grid = sim["grid"]
    _vmax = max(float(sim["snaps"][0][0].max()), 1e-12)

    fig1, ax1 = plt.subplots(figsize=(5, 4.5), constrained_layout=True)
    plotting.field(ax1, _eta, grid=_grid, signed=True, title=f"$\\eta$ at $t={sim['t'][_i]:.1f}$",
                   vmin=-_vmax, vmax=_vmax)
    _skip = max(1, _grid.n // 16)
    ax1.quiver(_grid.x[::_skip, ::_skip], _grid.y[::_skip, ::_skip],
               _u[::_skip, ::_skip], _v[::_skip, ::_skip],
               scale=4, width=0.004, alpha=0.7, color="k")
    mo.vstack([frame_ui, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Energy partition

        Initially all energy is potential ($\eta_0\neq0$, $\mathbf u=0$).
        Gravity-wave radiation converts it toward kinetic energy as the waves
        propagate outward; the balanced remnant left behind settles into a
        roughly fixed $KE$/$PE$ split.

        **A classical, slightly scandalous fact** (Rossby 1938; Gill §7.2):
        adjustment is energetically *wasteful*. For the textbook step-profile
        case, of the potential energy released by flattening the interface,
        only **one third** ends up as kinetic energy of the balanced state —
        the other two thirds is radiated away by gravity waves, gone for
        good. Balance is not an energy-conserving rearrangement; it is what
        remains after the fluid has paid a large wave tax. Watch for it
        below: the total-energy curve (dashed) stays flat because the
        *domain* is closed, but the energy near the disturbance drops as the
        wave front carries its share outward.
        """
    )
    return


@app.cell(hide_code=True)
def _(plt, sim):
    # --- diagnostic 2: energy partition time series -------------------------
    fig2, ax2 = plt.subplots(figsize=(6, 4), constrained_layout=True)
    _TE = sim["KE"] + sim["PE"]
    ax2.plot(sim["t"], sim["KE"] / _TE[0], lw=2, label="$KE/E_0$")
    ax2.plot(sim["t"], sim["PE"] / _TE[0], lw=2, color="crimson", label="$PE/E_0$")
    ax2.plot(sim["t"], _TE / _TE[0], lw=1.5, ls="--", color="gray", label="$TE/E_0$")
    ax2.set_xlabel("$t\\ (f_0^{-1})$"); ax2.set_ylabel("normalized energy")
    ax2.set_ylim(0, 1.15); ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Potential vorticity: what does and doesn't change

        Left/middle: $q=\zeta-\hat\eta$ at $t=0$ and at the final time — on
        the $f$-plane ($\hat\beta=0$) these should be visually identical (the
        max pointwise difference is printed below, and is at machine
        precision). Right: a slice through the domain center comparing the
        **predicted** balanced height $\eta_{bal}=$ `invert_pv(q0)` — computed
        from the $t=0$ field alone, before integrating a single step — against
        the **actual** simulated $\eta$ at the final time.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo, np, plotting, plt, sim):
    # --- diagnostic 3: PV conservation + PV-inversion prediction -----------
    _grid = sim["grid"]
    _q0, _qf, _eta_bal = sim["q0"], sim["qf"], sim["eta_bal"]
    _eta_f = sim["snaps"][-1][0]
    _mid = _grid.n // 2
    _dq = float(np.max(np.abs(_qf - _q0)))

    fig3, axs3 = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    _qmax = max(float(np.abs(_q0).max()), 1e-12)
    plotting.field(axs3[0], _q0, grid=_grid, signed=True, title="$q$ at $t=0$",
                   vmin=-_qmax, vmax=_qmax, cbar=False)
    plotting.field(axs3[1], _qf, grid=_grid, signed=True,
                   title=f"$q$ at $t={sim['t'][-1]:.1f}$", vmin=-_qmax, vmax=_qmax)

    _y = _grid.y[0, :]
    axs3[2].plot(_y, _eta_bal[_mid, :], lw=2, label="$\\eta_{bal}$ (predicted, from $q_0$)")
    axs3[2].plot(_y, _eta_f[_mid, :], lw=2, ls="--", color="crimson",
                 label="$\\eta$ (simulated, final $t$)")
    axs3[2].set_xlabel("$y$"); axs3[2].set_title("slice through domain center")
    axs3[2].legend(fontsize=8); axs3[2].grid(alpha=0.3)

    mo.vstack([
        mo.md(f"**max $|q(t_f)-q(0)|$ = {_dq:.2e}**  (max $|q_0|$ = {_qmax:.2f})"),
        fig3,
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Vary $\sigma/L_R$.** Compare the "small" and "large" presets: how
          much of $\eta_0$'s amplitude survives in $\eta_{bal}$ in each case?
          *Do the estimate first:* a Gaussian of width $\sigma$ projects onto
          wavenumbers $k\sim1/\sigma$, and the inversion multiplies each by
          $1/(k^2+1)$. For $\sigma=0.5$ that factor is $\sim1/5$; for
          $\sigma=3$ it is $\sim0.9$. Check the slice panel against these
          numbers.
        - **Let the waves finish radiating.** Increase $T$ in the small- or
          large-scale preset: does the residual between the predicted
          $\eta_{bal}$ and the simulated $\eta$ at the domain center shrink?
          It should — the prediction is only "wrong" by whatever wave energy
          hasn't left yet. (But don't push $T$ past $\sim20$: the domain is
          periodic, and the waves you radiated will wrap around and come
          back. That's an artifact of the box, not physics.)
        - **Break exact PV conservation.** Switch to the $\beta$-plane preset:
          $q$ at $t=0$ and at the final time are no longer identical. Where in
          the domain is $|q(t_f)-q(0)|$ largest? The $\beta$-plane source is
          $q_t=-\hat\beta\hat v$, so the change should concentrate where the
          *meridional velocity* has been strongest, not where $\eta$ is
          largest — check against the velocity arrows.
        - **Find the drift.** On the $\beta$-plane preset, track the location
          of the $\eta$ extremum on the snapshot browser as $t$ increases —
          confirm it moves westward (decreasing $x$), consistent with
          $\omega_R=-\hat\beta k/(k^2+1)$. Estimate its speed from the frames
          and compare with $\hat\beta/(k^2+1)$ using $k\sim1/\sigma$.

        ### What you should have seen

        A height bump released from rest splits into outward-radiating
        inertia-gravity waves plus a stationary geostrophic remnant. The
        remnant's amplitude relative to $\eta_0$ grows with $\sigma/L_R$: for
        $\sigma\ll L_R$ almost nothing survives (mostly gravity waves); for
        $\sigma\gg L_R$ most of $\eta_0$ survives basically unchanged. In both
        cases $q$ is exactly conserved and `invert_pv(q0)` — computed **before
        any time-stepping** — already predicts the balanced end state. Turning
        on $\hat\beta$ breaks pointwise PV conservation and sends the balanced
        remnant drifting slowly westward: geostrophic adjustment on a
        $\beta$-plane never truly comes to rest.

        **Where this goes next.** This chapter is the book's hinge. The
        "PV is conserved, PV is invertible, the rest is waves" logic you
        just watched is *the* organizing idea of large-scale GFD: Ch. 7
        makes invertibility the star, Ch. 8 derives the equation the
        balanced remnant obeys once it's free to evolve (quasi-geostrophy),
        and Ch. 9 follows the Rossby waves. When observational
        oceanographers see a mesoscale eddy that has survived for months,
        they are looking at the $\sigma\gtrsim L_R$ preset's endgame.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
