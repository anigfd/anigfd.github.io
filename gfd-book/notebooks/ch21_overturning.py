import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 21 — Buoyancy-Driven / Overturning Circulation

        **Physical question.** The Hadley cell and the ocean's meridional
        overturning circulation look nothing alike on a map, but they're
        driven by the same mechanism: heat a fluid unevenly along one
        boundary and it organizes into a single overturning cell — rising
        where it's heated, sinking where it's cooled. After this chapter
        you should be able to explain why that single-cell structure is
        the generic outcome (not convective rolls), and derive the classic
        "abyssal recipe" balance that sets the deep ocean's stratification.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        **Abyssal recipes (Munk 1966).** Deep water sinks at high
        latitudes and slowly upwells everywhere else; a steady 1D balance
        between that upwelling $w$ and diapycnal diffusion $\kappa$ sets
        the interior stratification (symbols as in
        [NOTATION](../notation)):

        $$w\frac{\partial T}{\partial z}=\kappa\frac{\partial^2T}{\partial z^2}.$$

        **The overturning cell.** Differential heating along one boundary
        of a 2D Boussinesq fluid — warm at one end, cool at the other —
        drives a vorticity source through the buoyancy term, exactly as in
        Ch. 14, but with the *sign* of the forcing now varying
        horizontally instead of only vertically:

        $$\frac{\partial\zeta}{\partial t}+J(\psi,\zeta)=Pr\nabla^2\zeta+Pr\!\cdot\!Ra\,
          \frac{\partial\theta}{\partial x},\qquad
          \frac{\partial\theta}{\partial t}+J(\psi,\theta)=\nabla^2\theta+Q(x,z),
          \qquad \nabla^2\psi=\zeta.$$
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
    from gfdlib import convection, timestep, overturning
    return convection, np, overturning, plt, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        **Part A** is a small linear ODE with constant coefficients,
        solved exactly (`gfdlib.overturning.abyssal_profile`) via a 2×2
        `numpy.linalg.solve` for the boundary-condition constants — cheap
        enough to be fully reactive. **Part B** reuses
        `gfdlib.convection.ChannelGrid`'s already-validated Boussinesq
        vorticity-streamfunction machinery unchanged (periodic $x$,
        Dirichlet-zero $z$, tridiagonal Poisson solve), but replaces
        Rayleigh-Bénard's uniform bottom heating with a differential,
        surface-concentrated heating field
        (`gfdlib.overturning.surface_heating`) — warm at one end of the
        domain, cool at the other — the classic "horizontal convection"
        mechanism (Rossby 1965) for a single overturning cell.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"## Part A — Abyssal recipes: the deep-ocean balance")
    return


@app.cell(hide_code=True)
def _(mo):
    w_ui = mo.ui.slider(0.05, 5.0, step=0.05, value=1.0, label="upwelling $w$", show_value=True)
    kappa_ui = mo.ui.slider(0.05, 1.0, step=0.05, value=0.3, label="diffusivity $\\kappa$",
                            show_value=True)
    mo.vstack([
        mo.hstack([w_ui, kappa_ui], justify="start"),
        mo.md(r"> $T(0)=0$ (deep/bottom water), $T(1)=1$ (surface water). "
              r"The ratio $w/\kappa$ is the only thing that matters — it's "
              r"the inverse of a **Péclet number** comparing advection to "
              r"diffusion."),
    ])
    return kappa_ui, w_ui


@app.cell(hide_code=True)
def _(kappa_ui, np, overturning, plt, w_ui):
    _z = np.linspace(0, 1, 200)
    _T = overturning.abyssal_profile(_z, w_ui.value, kappa_ui.value, H=1.0)
    _T_diffusive = overturning.abyssal_profile(_z, 1e-6, kappa_ui.value, H=1.0)

    fig0, ax0 = plt.subplots(figsize=(5, 5), constrained_layout=True)
    ax0.plot(_T_diffusive, _z, "k--", lw=1.2, label="pure diffusion ($w\\to0$)")
    ax0.plot(_T, _z, lw=2.5, color="crimson",
             label=f"$w/\\kappa$={w_ui.value / kappa_ui.value:.2f}")
    ax0.set_xlabel("$T(z)$"); ax0.set_ylabel("$z$ (0=bottom, 1=surface)")
    ax0.set_title("Munk's abyssal-recipe profile"); ax0.legend(fontsize=9); ax0.grid(alpha=0.3)
    fig0
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Push $w/\kappa$ up: the profile stays close to the **bottom**
        value through most of the water column, with the transition to
        surface properties compressed into a thin layer near $z=1$ — the
        real deep ocean's near-uniform abyssal temperature, capped by a
        thin thermocline. Munk (1966) ran this argument in reverse: given
        an *observed* stratification and a $\kappa$ measured from tracer-
        release experiments, the balance above tells you $w$ — the
        globally-averaged upwelling rate needed to close the ocean's mass
        budget, one of the first quantitative constraints on the deep MOC.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"## Part B — A minimal overturning cell")
    return


@app.cell(hide_code=True)
def _(mo):
    Ra_ui = mo.ui.slider(1e4, 5e5, step=1e4, value=2e5, label="$Ra$", show_value=True)
    Q0_ui = mo.ui.slider(0.5, 10.0, step=0.5, value=5.0, label="heating amplitude $Q_0$",
                         show_value=True)
    nx_ui = mo.ui.dropdown(options={"32": 32, "48": 48}, value="32", label="resolution $n_x$")
    T_ui = mo.ui.slider(0.05, 0.5, step=0.05, value=0.2, label="run time $T$", show_value=True)
    mo.vstack([
        mo.hstack([Ra_ui, Q0_ui], justify="start"),
        mo.hstack([nx_ui, T_ui], justify="start"),
    ])
    return Q0_ui, Ra_ui, T_ui, nx_ui


@app.cell(hide_code=True)
def _(mo):
    get_ot, set_ot = mo.state(None)
    return get_ot, set_ot


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Run overturning simulation")
    run_btn
    return (run_btn,)


@app.cell
def _(
    Q0_ui, Ra_ui, T_ui, convection, mo, np, nx_ui, overturning, run_btn,
    set_ot, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run overturning simulation** above to start."), kind="info"
    ))

    _nx = nx_ui.value
    _nz = _nx // 2
    _Lx = 2.0
    _Ra, _Q0, _T, _Pr = Ra_ui.value, Q0_ui.value, T_ui.value, 1.0

    grid = convection.ChannelGrid(_nx, _nz, _Lx)
    Q = overturning.surface_heating(grid, _Q0)

    _rng = np.random.default_rng(0)
    _zeta = np.zeros((_nx, _nz))
    _theta = 1e-3 * _rng.standard_normal((_nx, _nz))
    _state = np.stack([_zeta, _theta])

    _dt = 0.15 * min(grid.dx, grid.dz) ** 2
    _nsteps = max(20, int(_T / _dt))
    _nsub = max(1, _nsteps // 60)

    def _rhs(t, s):
        return overturning.rhs_overturning(s, grid, _Ra, _Pr, Q)

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

    set_ot({
        "snaps_th": _snaps_th, "snaps_psi": _snaps_psi,
        "t": np.array(_t_hist), "grid": grid, "Ra": _Ra, "Q0": _Q0,
    })

    mo.callout(mo.md(
        f"**Done.** $Ra={_Ra:.0f}$, $Q_0={_Q0}$, grid {_nx}×{_nz}, {_nsteps} steps "
        f"→ {len(_snaps_th)} frames."
    ), kind="success")
    return


@app.cell
def _(get_ot, mo):
    ot_sim = get_ot()
    mo.stop(ot_sim is None, mo.callout(
        mo.md("Results will appear here once the simulation has run."), kind="warn"
    ))
    return (ot_sim,)


@app.cell(hide_code=True)
def _(mo, ot_sim):
    ot_frame = mo.ui.slider(0, len(ot_sim["snaps_th"]) - 1, step=1,
                            value=len(ot_sim["snaps_th"]) - 1, label="frame", show_value=True)
    return (ot_frame,)


@app.cell(hide_code=True)
def _(mo, np, ot_frame, ot_sim, plt):
    _i = ot_frame.value
    _g = ot_sim["grid"]
    _th = ot_sim["snaps_th"][_i]; _psi = ot_sim["snaps_psi"][_i]
    _XX, _ZZ = np.meshgrid(_g.x, _g.z, indexing="ij")

    fig1, axs1 = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    _vm = max(float(np.abs(_th).max()), 1e-10)
    _im0 = axs1[0].pcolormesh(_XX, _ZZ, _th, cmap="RdBu_r", vmin=-_vm, vmax=_vm, shading="auto")
    plt.colorbar(_im0, ax=axs1[0], fraction=0.046, label=r"$\theta$")
    axs1[0].set_title(f"buoyancy  $t={ot_sim['t'][_i]:.3f}$")
    axs1[0].set_xlabel("$x$ (heated $\\to$ cooled $\\to$ heated)"); axs1[0].set_ylabel("$z$")

    _vm2 = max(float(np.abs(_psi).max()), 1e-10)
    _im1 = axs1[1].pcolormesh(_XX, _ZZ, _psi, cmap="PuOr", vmin=-_vm2, vmax=_vm2, shading="auto")
    plt.colorbar(_im1, ax=axs1[1], fraction=0.046, label=r"$\psi$")
    axs1[1].set_title(f"streamfunction  $t={ot_sim['t'][_i]:.3f}$")
    axs1[1].set_xlabel("$x$"); axs1[1].set_ylabel("$z$")
    mo.vstack([ot_frame, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        A single, coherent overturning cell — not the small convective
        rolls of Ch. 14 — rising over the heated column and sinking over
        the cooled one, exactly as in the abyssal-recipe intuition of
        Part A: warm, light fluid rises; cold, dense fluid sinks; the
        return flow closes the loop along the surface and the bottom.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Confirm the circulation sense.** At the frame slider's final
          position, the streamfunction's sign should show rising motion
          ($\partial\psi/\partial x>0$) directly under the heated column
          and sinking ($\partial\psi/\partial x<0$) under the cooled one —
          the same check `gfdlib`'s test suite makes automatically.
        - **Push $Q_0$.** Stronger differential heating drives a stronger
          cell (larger $|\psi|$) — does the RELATIONSHIP look linear, or
          does it saturate? This connects to the long-standing "does MOC
          strength scale linearly with the pole-to-equator buoyancy
          contrast" question in real ocean circulation theory.
        - **Connect Part A to Part B.** In Part A, larger $w/\kappa$ meant
          a thinner surface transition layer. In Part B's temperature
          field, does a stronger cell (larger $Q_0$ or $Ra$) similarly
          compress the buoyancy anomaly into a thinner layer near the
          surface?

        ### What you should have seen

        Two versions of the same idea: a 1D balance (Part A) showing that
        upwelling against diffusion sets a boundary-layer-like
        stratification, and a fully 2D, self-consistent circulation
        (Part B) showing WHY that upwelling exists in the first place —
        differential heating along a boundary organizes into one
        overturning cell, not many small ones, because the domain only
        has ONE large-scale buoyancy contrast to respond to. The same
        mechanism, at utterly different scales, drives the atmosphere's
        Hadley cell and the ocean's meridional overturning circulation —
        the two phenomena that, together with wind-driven gyres (Ch. 20),
        set the shape of the entire general circulation.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
