import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 23 — Conservation Laws & Wave Activity

        **Physical question.** Ch. 13 promised that waves deposit momentum
        in the mean flow "once dissipation or breaking enters" — implying
        that *without* them, waves deposit nothing. That implication is a
        theorem, and this chapter proves it and then verifies it against
        the full nonlinear model: for conservative, small-amplitude waves,
        the zonal-mean flow and the wave's **pseudomomentum** change in
        exact opposition, $\partial_t(\bar u+A)=0$ pointwise in latitude —
        the **non-acceleration theorem**. Waves are not a source of mean
        momentum; they are a *loan* of it, repaid in full unless the wave
        dies. After this chapter you should be able to derive $A$, state
        what breaks non-acceleration, and say why that breakdown is how
        the real stratosphere's winds are driven.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### Noether, applied to waves: pseudomomentum

        Ch. 22 read conservation laws off symmetries. The relevant symmetry
        here is **zonal translation**: a zonally-symmetric mean state looks
        the same shifted in $x$, so there is a conserved momentum-like
        charge associated with the *disturbance* riding on it — the
        **pseudomomentum**, or wave activity. For barotropic dynamics its
        small-amplitude density has a strikingly simple form (symbols as in
        [NOTATION](../notation)):

        $$A=\frac{\overline{q'^2}}{2\,\bar q_y},$$

        with $q'$ the PV disturbance, $\bar q_y=\beta-\bar U''$ the mean PV
        gradient, and the overbar a zonal average. Enstrophy over gradient
        — a measure of how far PV contours have been displaced (for a
        displacement $\eta$ of a PV contour, $q'\approx-\eta\,\bar q_y$, so
        $A\approx\tfrac12\bar q_y\overline{\eta^2}$: wave activity is the
        mean-square *sideways displacement* of the flow's PV contours,
        weighted by the gradient they cut across).

        ### Three short derivations, one theorem

        **(1) How $A$ changes.** Linearize PV conservation about the mean:
        $q'_t+\bar U q'_x+v'\bar q_y=0$. Multiply by $q'/\bar q_y$ and
        zonal-average (the $\bar Uq'_x$ term averages to zero by
        periodicity):

        $$\frac{\partial A}{\partial t}=-\,\overline{v'q'}.$$

        **(2) The Taylor identity.** The eddy PV flux is itself a total
        derivative: writing $q'=\nabla^2\psi'$ and zonal-averaging,

        $$\overline{v'q'}=\overline{\psi'_x\nabla^2\psi'}
          =-\frac{\partial}{\partial y}\,\overline{u'v'}$$

        (integrate by parts in $x$ twice; every pure-$x$ derivative of a
        product averages to zero on a periodic domain). The PV flux and
        the momentum-flux convergence are the *same object*.

        **(3) How $\bar u$ changes.** Zonal-average the momentum equation:
        $\bar u_t=-\partial_y\overline{u'v'}$, which by (2) equals
        $\overline{v'q'}$, which by (1) equals $-\partial A/\partial t$.
        Add:

        $$\boxed{\ \frac{\partial}{\partial t}\big(\bar u+A\big)=0\ }
        \qquad\text{pointwise in }y.$$

        This is the **non-acceleration theorem** (Charney-Drazin 1961;
        Eliassen-Palm; Andrews-McIntyre 1976): conservative,
        small-amplitude waves cannot change the mean flow *at all* except
        by lending it exactly what their own pseudomomentum gives up. A
        steady wave field ($A_t=0$) forces $\bar u_t=0$ no matter how
        vigorous the wave is.

        ### What breaks it — and why that is the interesting part

        The theorem has exactly two exits, and both are physics rather
        than fine print. **Transience**: while a wave packet grows through
        a latitude, $A$ rises and $\bar u$ dips (the loan); when the packet
        moves on, both return. **Dissipation or breaking**: if the wave
        dies at some latitude — Ch. 13's critical layer being the
        canonical place — its pseudomomentum is handed over *permanently*.
        The stratospheric quasi-biennial oscillation is this exit running
        as an engine: waves launched below break at descending critical
        lines and reverse the winds, over and over. "No acceleration
        without dissipation" is the single most consequential negative
        result in wave-mean theory.
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
    from gfdlib import spectral, timestep, plotting, wavemean
    return np, plotting, plt, spectral, timestep, wavemean


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        The theorem is *linear* theory; the honest test is against the
        **full nonlinear model**. Ch. 18's exact pseudo-spectral machinery
        (`gfdlib.spectral.Grid`, `gfdlib.timestep.ifrk4_step`, exact
        $\beta$ propagator, near-zero hyperviscosity) evolves a small
        Rossby-wave packet superposed on a sinusoidal shear
        $\bar U(y)=U_0\cos y$; `gfdlib.wavemean.wave_activity` computes
        $\bar u(y)$, $A(y)$, and $\bar q_y(y)$ from each snapshot. The
        theorem then predicts the two O(amplitude$^2$) changes —
        each individually tiny — cancel to the next order. The same check
        runs automatically in `tests/test_gfdlib.py`, where the residual
        comes out below 1% of either term. Note the design constraint
        baked into the defaults: $\beta>U_0$ keeps $\bar q_y>0$ everywhere,
        since $A$ is undefined where the mean PV gradient vanishes.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    a_ui = mo.ui.slider(0.01, 0.4, step=0.01, value=0.05,
                        label="packet amplitude $a$", show_value=True)
    k0_ui = mo.ui.slider(3, 10, step=1, value=6,
                         label="carrier wavenumber $k_0$", show_value=True)
    U0_ui = mo.ui.slider(0.0, 2.0, step=0.1, value=1.0,
                         label="shear amplitude $U_0$", show_value=True)
    beta_ui = mo.ui.slider(2.0, 10.0, step=0.5, value=5.0,
                           label="$\\beta$", show_value=True)
    nu_ui = mo.ui.slider(0.0, 2e-5, step=1e-6, value=0.0,
                         label="hyperviscosity $\\nu$ (0 = conservative)", show_value=True)
    T_ui = mo.ui.slider(1.0, 12.0, step=1.0, value=4.0,
                        label="run time $T$", show_value=True)
    n_ui = mo.ui.dropdown(options={"96": 96, "128": 128}, value="128",
                          label="resolution $n$")
    mo.vstack([
        mo.hstack([a_ui, k0_ui, U0_ui], justify="start"),
        mo.hstack([beta_ui, nu_ui], justify="start"),
        mo.hstack([T_ui, n_ui], justify="start"),
        mo.md(r"> Keep $\beta>U_0$ so $\bar q_y=\beta+U_0\cos y>0$ "
              r"everywhere. The $\nu$ slider is the theorem's OFF switch: "
              r"conservative at 0, dissipative above."),
    ])
    return T_ui, U0_ui, a_ui, beta_ui, k0_ui, n_ui, nu_ui


@app.cell(hide_code=True)
def _(mo):
    get_sim, set_sim = mo.state(None)   # results survive the run button resetting
    return get_sim, set_sim


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Run")
    run_btn
    return (run_btn,)


@app.cell
def _(
    T_ui, U0_ui, a_ui, beta_ui, k0_ui, mo, n_ui, np, nu_ui, run_btn, set_sim,
    spectral, timestep, wavemean,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run** above to start."), kind="info"
    ))

    _n, _a, _k0 = n_ui.value, a_ui.value, k0_ui.value
    _U0, _beta, _nu, _T = U0_ui.value, beta_ui.value, nu_ui.value, T_ui.value
    grid = spectral.Grid(_n)
    _x0 = _y0 = np.pi
    _sigma = 0.5

    _zeta_bar0 = _U0 * np.sin(grid.y)          # U(y) = U0*cos(y) => zeta = U0*sin(y)
    _env = np.exp(-((grid.x - _x0) ** 2 + (grid.y - _y0) ** 2) / (2 * _sigma ** 2))
    _zeta0 = _zeta_bar0 + _a * _env * np.cos(_k0 * (grid.x - _x0))
    _zh = grid.fft(_zeta0) * grid.dealias

    _L_op = -_nu * grid.k2 ** 2 + 1j * _beta * grid.kx * grid.k2_inv

    def _rhs_nl(t, F):
        return -grid.jacobian(grid.invert_laplacian(F), F)

    _dt = 0.004
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 60)

    _t_hist, _ubars, _As = [], [], []
    _snap0 = grid.ifft(_zh).astype(np.float32)
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _ub, _A, _ = wavemean.wave_activity(grid.ifft(_zh), grid, _beta)
            _t_hist.append(_t); _ubars.append(_ub); _As.append(_A)
        if _s < _nsteps:
            _zh = timestep.ifrk4_step(_zh, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    set_sim({
        "t": np.array(_t_hist),
        "ubar": np.array(_ubars),          # (frames, n)
        "A": np.array(_As),
        "zeta0": _snap0, "zeta1": grid.ifft(_zh).astype(np.float32),
        "grid": grid, "y": grid.y[0, :],
    })
    mo.callout(mo.md(f"**Done.** {_nsteps} steps on a {_n}×{_n} grid."), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the run has completed."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(mo, plotting, plt, sim):
    fig0, axs0 = plt.subplots(1, 2, figsize=(9.5, 4), constrained_layout=True)
    plotting.field(axs0[0], sim["zeta0"], grid=sim["grid"], signed=True,
                   title="$\\zeta$ at $t=0$: shear + packet", cbar=False)
    plotting.field(axs0[1], sim["zeta1"], grid=sim["grid"], signed=True,
                   title=f"$\\zeta$ at $t={sim['t'][-1]:.1f}$")
    fig0
    return


@app.cell(hide_code=True)
def _(np, plt, sim):
    # --- the theorem, pointwise in y ----------------------------------------
    _dU = sim["ubar"][-1] - sim["ubar"][0]
    _dA = sim["A"][-1] - sim["A"][0]

    fig1, axs1 = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axs1[0].plot(_dU, sim["y"], lw=2, color="#2563eb", label="$\\Delta\\bar u$")
    axs1[0].plot(-_dA, sim["y"], "--", lw=2, color="crimson", label="$-\\Delta A$")
    axs1[0].set_xlabel("change"); axs1[0].set_ylabel("$y$")
    axs1[0].set_title("mean-flow change vs. minus pseudomomentum change")
    axs1[0].legend(fontsize=9); axs1[0].grid(alpha=0.3)

    _resid = _dU + _dA
    axs1[1].plot(_resid, sim["y"], lw=2, color="k")
    axs1[1].set_xlim(axs1[0].get_xlim())
    axs1[1].set_xlabel("$\\Delta(\\bar u+A)$"); axs1[1].set_title(
        f"residual (same axis scale)  ·  max ratio "
        f"{np.abs(_resid).max() / max(np.abs(_dU).max(), 1e-30):.1%}")
    axs1[1].grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read it:** the left panel's two curves — computed from
        completely different fields (a velocity average vs. an
        enstrophy-weighted PV statistic) — should lie on top of each other
        at every latitude. The right panel plots their sum on the *same*
        axis scale: near-blank is the theorem. The percentage in its title
        is the residual relative to the signal; at small amplitude and
        $\nu=0$ it sits at the sub-percent level expected from the
        neglected $O(a^3)$ terms.
        """
    )
    return


@app.cell(hide_code=True)
def _(np, plt, sim):
    # --- Hovmoller: the loan being made and repaid --------------------------
    _dU_t = sim["ubar"] - sim["ubar"][0]
    _res_t = _dU_t + (sim["A"] - sim["A"][0])
    _vm = max(np.abs(_dU_t).max(), 1e-30)

    fig2, axs2 = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for _ax, _f, _ttl in [(axs2[0], _dU_t, "$\\Delta\\bar u(y,t)$"),
                          (axs2[1], _res_t, "$\\Delta(\\bar u+A)(y,t)$ — same color scale")]:
        _im = _ax.imshow(_f.T, origin="lower", aspect="auto", cmap="RdBu_r",
                         vmin=-_vm, vmax=_vm,
                         extent=[sim["t"][0], sim["t"][-1], sim["y"][0], sim["y"][-1]])
        _ax.set_xlabel("$t$"); _ax.set_ylabel("$y$"); _ax.set_title(_ttl)
    fig2.colorbar(_im, ax=axs2, shrink=0.85)
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Push the amplitude.** Raise $a$ toward 0.4 and watch the
          residual percentage climb: the theorem is an $O(a^2)$ statement
          and its error is $O(a^3)$ — so the *relative* residual should
          grow roughly linearly in $a$. Confirm that against two or three
          amplitudes.
        - **Throw the OFF switch.** Set $\nu$ to its maximum and rerun:
          the left panel's curves now visibly separate — dissipation
          destroys pseudomomentum without paying it back to $\bar u$
          (it hands the momentum to the small scales instead). This is
          the theorem failing *for the stated reason*, on demand.
        - **Watch the loan.** In the Hovmöller pair, the left panel shows
          $\Delta\bar u$ developing structure where the packet sits; the
          right panel stays blank *through* all of it — the compensation
          is not just an end-state accident but holds at every instant.
        - **Court the singularity.** Lower $\beta$ toward $U_0$: as
          $\min_y\bar q_y\to0$, the $A$ diagnostic develops spikes at the
          gradient's zero crossings and the residual deteriorates — not a
          bug, but the small-amplitude theory failing where its
          denominator vanishes (Ch. 13's critical-layer physics lives at
          exactly such places).

        ### What you should have seen

        Two curves computed from entirely different quantities agreeing
        pointwise, and their sum flat-lining through the whole run: the
        mean flow changed by *exactly* minus the wave's pseudomomentum
        change, latitude by latitude, instant by instant — until you
        switched on dissipation, at which point the books stopped
        balancing precisely as advertised. This closes the arc begun in
        Ch. 13: waves carry momentum as a loan (non-acceleration), and
        only their death — by breaking at critical layers, by radiative
        damping — converts the loan to a grant. That conversion drives
        the QBO, decelerates the winter stratospheric jet (sudden
        warmings), and shapes the ocean's eddy-driven circulation. Ch. 24
        closes the book by asking the complementary question: not how
        waves force the slow flow, but what "the slow flow" even is.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
