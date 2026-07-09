import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 8 — The Quasi-Geostrophic Approximation

        **Physical question.** Take Ch. 6's rotating shallow water and
        formally expand it for slow, balanced flow — what survives? After
        this chapter you should be able to write the single-layer QGPV
        equation from RSW, explain what the deformation radius does to an
        isolated vortex's velocity field and to the Rossby-wave dispersion
        relation, and watch a vortex drift westward while shedding a trailing
        wake of Rossby waves — a *beta-gyre*.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        Starting from Ch. 6's rotating shallow water and formally expanding
        in small Rossby number $Ro=U/(f_0L)$: the leading-order flow is
        geostrophic, and the $O(Ro)$ ageostrophic correction is exactly what
        advances the leading-order field forward in time. Carrying this
        through collapses RSW to a single evolution equation for the
        **quasi-geostrophic potential vorticity** (symbols as in
        [NOTATION](../notation)):

        $$q_t+J(\psi,q)=0,\qquad q=\nabla^2\psi-\frac{\psi}{L_R^2}+\beta y.$$

        Compare this to Ch. 7's barotropic $q=\zeta+\beta y$: QG adds a
        **stretching term** $-\psi/L_R^2$, from the free-surface deformation
        that RSW allows and pure 2D barotropic flow does not. Nondimensionalized
        by $L_R$ so $L_R=1$, the balance operator becomes
        $(\nabla^2-1)\psi=q-\beta y$ — **exactly** Ch. 6's shallow-water PV
        Helmholtz operator. QG inversion is nothing new: it's the same
        operator, applied to a field that's now free to evolve nonlinearly on
        a $\beta$-plane.
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
    from gfdlib import spectral, timestep, plotting, pv, qg
    return np, plotting, plt, pv, qg, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Reuses Ch. 18's exact pseudo-spectral machinery
        (`gfdlib.spectral.Grid`, `gfdlib.timestep.ifrk4_step`) with two
        changes: $\psi$ comes from `gfdlib.shallowwater.invert_pv` (the
        Helmholtz solve, not a pure Poisson solve), and the linear operator's
        Rossby-wave term uses `gfdlib.qg.dispersion_omega`'s $1/(k^2+1)$
        rather than $1/k^2$ — both are exact, so hyperviscosity and
        $\beta$-propagation remain unconditionally stable via the
        integrating factor. `gfdlib.pv.gaussian_blob` places the vortex.
        """
    )
    return


@app.cell(hide_code=True)
def _(qg, np, plt):
    # --- diagnostic 1: theory curve, no simulation needed -------------------
    _beta = 1.0
    _k = np.linspace(0.02, 5, 400)
    _omega_qg = qg.dispersion_omega(_k, 0.0, _beta)
    _omega_bt = -_beta * _k / _k ** 2   # barotropic (ch18), for comparison

    fig0, ax0 = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax0.plot(_k, _omega_qg, lw=2, label="QG: $-\\beta k/(k^2+1)$")
    ax0.plot(_k, _omega_bt, lw=1.5, ls="--", color="gray", label="barotropic: $-\\beta k/k^2$")
    ax0.axhline(-_beta / 2, color="crimson", ls=":", lw=1, label="$-\\beta/2$ (QG cap)")
    ax0.axvline(1.0, color="crimson", ls=":", lw=1)
    ax0.set_ylim(-2 * _beta, 0.2)
    ax0.set_xlabel("$k$ (units of $1/L_R$)"); ax0.set_ylabel("$\\omega$ ($\\beta=1$)")
    ax0.set_title("Rossby-wave dispersion: QG vs. barotropic")
    ax0.legend(fontsize=9); ax0.grid(alpha=0.3)
    fig0
    return


@app.cell(hide_code=True)
def _(mo):
    scen_md = mo.md(
        r"""
        Above: the barotropic dispersion relation (ch18) diverges as
        $k\to0$; the QG relation is capped at $\beta/2$, attained at
        $k=1/L_R$ — the deformation radius sets an intrinsic scale below
        which Rossby waves simply can't get any faster.
        """
    )
    scen_md
    return


@app.cell(hide_code=True)
def _(mo):
    amp_ui = mo.ui.slider(0.3, 2.0, step=0.1, value=1.0,
                          label="vortex amplitude", show_value=True)
    sigma_ui = mo.ui.slider(0.4, 2.5, step=0.1, value=1.0,
                            label="vortex radius (units of $L_R$)", show_value=True)
    beta_ui = mo.ui.slider(0.0, 1.0, step=0.05, value=0.3,
                           label="$\\beta$", show_value=True)
    n_ui = mo.ui.dropdown(options={"64": 64, "96": 96, "128": 128}, value="96",
                          label="resolution $n$")
    T_ui = mo.ui.slider(5.0, 60.0, step=5.0, value=30.0,
                        label="run time $T$", show_value=True)

    mo.vstack([
        mo.hstack([amp_ui, sigma_ui, beta_ui], justify="start"),
        mo.hstack([n_ui, T_ui], justify="start"),
        mo.md(r"> Domain is $20\,L_R$ across. Try $\sigma\ll L_R$ (a tight "
              r"vortex, strongly dispersive) vs. $\sigma\gg L_R$ (broad, "
              r"nearly non-dispersive)."),
    ])
    return T_ui, amp_ui, beta_ui, n_ui, sigma_ui


@app.cell(hide_code=True)
def _(mo):
    get_sim, set_sim = mo.state(None)   # results survive the run button resetting
    return get_sim, set_sim


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Run evolution")
    run_btn
    return (run_btn,)


@app.cell
def _(
    T_ui, amp_ui, beta_ui, mo, n_ui, np, pv, run_btn, set_sim, sigma_ui,
    spectral, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run evolution** above to start."), kind="info"
    ))

    _n, _amp, _sigma, _beta, _T = (
        n_ui.value, amp_ui.value, sigma_ui.value, beta_ui.value, T_ui.value
    )
    _L = 20.0
    grid = spectral.Grid(_n, L=_L)
    _xc = _yc = 0.5 * _L
    _q0 = pv.gaussian_blob(grid, _xc, _yc, _amp, _sigma)

    _nu, _nnu = 1e-4, 2
    _q_hat = grid.fft(_q0)
    _L_op = -_nu * grid.k2 ** _nnu + 1j * _beta * grid.kx / (grid.k2 + 1.0)

    def _rhs_nl(t, qh):
        # same operator as shallowwater.invert_pv, applied in spectral space
        # directly (qh is already spectral inside the RK4 substeps)
        _psi_h = -qh / (grid.k2 + 1.0)
        return -grid.jacobian(_psi_h, qh)

    _dt = 0.01
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 80)

    _snaps, _t_hist, _centroid = [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _q = grid.ifft(_q_hat)
            _snaps.append(_q.astype(np.float32))
            _t_hist.append(_t)
            _w = np.clip(_q, 0, None)
            if _w.sum() > 1e-6:
                _cx = (grid.x * _w).sum() / _w.sum()
                _cy = (grid.y * _w).sum() / _w.sum()
            else:
                _cx, _cy = np.nan, np.nan
            _centroid.append((_cx, _cy))
        if _s < _nsteps:
            _q_hat = timestep.ifrk4_step(_q_hat, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    set_sim({
        "snaps": _snaps, "t": np.array(_t_hist), "grid": grid,
        "centroid": np.array(_centroid), "xc": _xc, "yc": _yc,
    })

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_snaps)} frames."
    ), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the evolution has run."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(mo, sim):
    frame_ui = mo.ui.slider(0, len(sim["snaps"]) - 1, step=1, value=len(sim["snaps"]) - 1,
                            label="frame", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, mo, plotting, plt, sim):
    # --- diagnostic: snapshot browser of the evolving PV field --------------
    _i = frame_ui.value
    fig1, ax1 = plt.subplots(figsize=(5.5, 5), constrained_layout=True)
    plotting.field(ax1, sim["snaps"][_i], grid=sim["grid"], signed=True,
                   title=f"$q$ at $t={sim['t'][_i]:.1f}$")
    mo.vstack([frame_ui, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### The beta-drift

        Trajectory of the vortex core (its PV centroid) over the run. A
        westward drift is the beta-gyre signature: the asymmetric dipolar
        correction the vortex develops on a $\beta$-plane self-advects it,
        even though nothing in the flow was given any initial translation.
        """
    )
    return


@app.cell(hide_code=True)
def _(plt, sim):
    fig2, ax2 = plt.subplots(figsize=(5.5, 4), constrained_layout=True)
    _c = sim["centroid"]
    ax2.plot(_c[:, 0], _c[:, 1], "-o", ms=3, lw=1.5)
    ax2.plot(sim["xc"], sim["yc"], "k*", ms=12, label="start")
    ax2.set_xlabel("$x$"); ax2.set_ylabel("$y$")
    ax2.set_title("vortex-core trajectory")
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    ax2.set_aspect("equal")
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Turn off $\beta$.** Set $\beta=0$: the vortex should sit still
          and simply diffuse very slowly — no drift, no wake. Then turn
          $\beta$ back up and watch the wake reappear.
        - **Tight vs. broad vortex.** Compare $\sigma=0.4$ and $\sigma=2.0$
          at fixed $\beta$: which sheds a more visible Rossby wave wake in
          the same run time, and why does that match the dispersion-relation
          panel (which $k$'s does each vortex size project onto)?
        - **Read the cap.** At your chosen $\beta$, use the dispersion panel
          to read off $\beta/2$ — that's the fastest possible Rossby-wave
          frequency at any wavenumber. Does the wake in the snapshot browser
          ever seem to move faster than that bound would allow?

        ### What you should have seen

        An isolated QG vortex is not a permanent structure: on a $\beta$-plane
        it drifts (generically westward, with a small meridional component)
        while continuously shedding Rossby waves into a trailing wake — the
        *beta-gyre*. Both effects trace back to the same modification QG adds
        to Ch. 7's barotropic PV: the $-\psi/L_R^2$ stretching term, which
        reuses Ch. 6's Helmholtz inversion exactly and caps the Rossby-wave
        frequency at $\beta/2$. This is the mechanism behind the observed
        westward propagation of oceanic mesoscale eddies.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
