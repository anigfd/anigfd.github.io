import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 10 — Two-Layer QG & Available Potential Energy

        **Physical question.** Ch. 3 split a flow into "columnar" and
        "sheared" parts to talk about Taylor-Proudman. A two-layer ocean or
        atmosphere lets you make that split *exact*: every flow decomposes
        uniquely into a **barotropic mode** (the depth-mean, which behaves
        like Ch. 18's ordinary 2D turbulence) and a **baroclinic mode** (the
        inter-layer difference, which carries all the available potential
        energy and behaves like Ch. 8's single-layer QG). After this chapter
        you should be able to derive that decoupling, and watch what
        happens when a nonlinear disturbance starts purely in one mode: it
        does not stay there.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The exact mode split

        Ch. 16 wrote the 2-layer QGPV as (equal-thickness layers, symbols as
        in [NOTATION](../notation)):

        $$q_1=\nabla^2\psi_1+F(\psi_2-\psi_1),\qquad
          q_2=\nabla^2\psi_2+F(\psi_1-\psi_2)$$

        (perturbation PV; the mean $\beta y$ is carried separately). Define
        the **barotropic** and **baroclinic** modes

        $$\psi_{bt}=\frac{\psi_1+\psi_2}{2}\ \text{(depth mean)},\qquad
          \psi_{bc}=\frac{\psi_1-\psi_2}{2}\ \text{(half the shear)}.$$

        Add and subtract the two PV equations. The $F$ terms are
        antisymmetric ($F(\psi_2-\psi_1)$ vs. $F(\psi_1-\psi_2)$, exactly
        opposite), so they **cancel** in the sum and **double** in the
        difference:

        $$q_{bt}=\frac{q_1+q_2}{2}=\nabla^2\psi_{bt},
          \qquad
          q_{bc}=\frac{q_1-q_2}{2}=\nabla^2\psi_{bc}-2F\psi_{bc}.$$

        Read these literally: $q_{bt}$ is **exactly** Ch. 7/18's barotropic
        PV — no $F$, no memory of layering at all — and $q_{bc}$ is
        **exactly** Ch. 8's single-layer QGPV, with deformation-radius
        parameter $2F$ in place of $1$. The two-layer ocean, split this way,
        is *literally two decoupled copies of chapters you have already
        built* — verified directly in `gfdlib.baroclinic.to_modes`
        (not just asserted; see `tests/test_gfdlib.py`).

        ### Available potential energy, revisited

        Ch. 16 computed $EPE=\tfrac12F\langle(\psi_1-\psi_2)^2\rangle
        =2F\langle\psi_{bc}^2\rangle$ without naming it this way: APE lives
        **entirely in the baroclinic mode** — a flow with $\psi_{bc}\equiv0$
        (purely barotropic) has exactly zero APE, by construction, no matter
        how energetic its barotropic circulation is.

        ### Where the split BREAKS

        The *linear* operators decouple exactly. The Jacobian nonlinearity
        in the full 2-layer equations does not:
        $J(\psi_1,q_1)$ and $J(\psi_2,q_2)$, rewritten in mode variables,
        produce **cross terms** coupling $bt$ and $bc$ — an eddy's advection
        of its own PV mixes the two modes even when nothing else does. This
        chapter's numerical experiment is built entirely to expose that one
        fact: seed a disturbance purely in the baroclinic mode and watch a
        barotropic mode appear from nothing but nonlinear advection.
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
    from gfdlib import spectral, timestep, plotting, baroclinic, pv
    return baroclinic, np, plotting, plt, pv, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Everything below reuses ch16's exact `gfdlib.baroclinic` machinery
        (`invert_2layer`, `rhs_2layer`, `gfdlib.spectral.Grid`,
        `gfdlib.timestep.ifrk4_step`) with **no mean shear**
        ($U_1=U_2=0$) — this is an unforced initial-value problem (Ch. 6's
        adjustment logic, in the 2-layer setting), not an instability
        calculation. The initial condition is built directly in mode space:
        `gfdlib.baroclinic.from_modes` converts a chosen $(q_{bt},q_{bc})$
        pair into the $(q_1,q_2)$ the solver actually steps, using
        `gfdlib.pv.gaussian_blob` for the baroclinic-mode vortex shape.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    F_ui = mo.ui.slider(0.3, 3.0, step=0.1, value=1.0, label="$F$", show_value=True)
    amp_ui = mo.ui.slider(0.3, 2.0, step=0.1, value=1.0, label="vortex amplitude", show_value=True)
    sigma_ui = mo.ui.slider(0.5, 2.5, step=0.1, value=1.0,
                            label="vortex radius", show_value=True)
    beta_ui = mo.ui.slider(0.0, 1.0, step=0.05, value=0.2, label="$\\beta$", show_value=True)
    T_ui = mo.ui.slider(5.0, 60.0, step=5.0, value=30.0, label="run time $T$", show_value=True)
    n_ui = mo.ui.dropdown(options={"64": 64, "96": 96, "128": 128}, value="96", label="resolution $n$")
    mo.vstack([
        mo.hstack([F_ui, amp_ui, sigma_ui], justify="start"),
        mo.hstack([beta_ui, T_ui, n_ui], justify="start"),
        mo.md(r"> The initial condition is PURE baroclinic mode: "
              r"$q_{bt}(t{=}0)\equiv0$ everywhere, $q_{bc}(t{=}0)=$ a "
              r"single Gaussian blob. Any barotropic energy that appears "
              r"had to come from the nonlinear cross-coupling, not the "
              r"initial condition."),
    ])
    return F_ui, T_ui, amp_ui, beta_ui, n_ui, sigma_ui


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
    F_ui, T_ui, amp_ui, baroclinic, beta_ui, mo, n_ui, np, pv, run_btn,
    set_sim, sigma_ui, spectral, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run** above to start."), kind="info"
    ))

    _n, _F, _amp, _sigma, _beta, _T = (
        n_ui.value, F_ui.value, amp_ui.value, sigma_ui.value, beta_ui.value, T_ui.value
    )
    grid = spectral.Grid(_n)
    _xc = _yc = 0.5 * grid.L

    _q_bc0 = pv.gaussian_blob(grid, _xc, _yc, _amp, _sigma)
    _q_bt0 = np.zeros_like(_q_bc0)
    _q1_0, _q2_0 = baroclinic.from_modes(_q_bt0, _q_bc0)

    _nu, _nnu = 5e-4, 2
    _state_hat = grid.fft(np.stack([_q1_0, _q2_0])) * grid.dealias
    _L_op = -_nu * grid.k2 ** _nnu

    def _rhs_nl(t, sh):
        return baroclinic.rhs_2layer(sh, grid, _F, 0.0, 0.0, _beta, r=0.0)

    def _mode_energies(state_hat):
        q1_hat, q2_hat = state_hat
        psi1_hat, psi2_hat = baroclinic.invert_2layer(q1_hat, q2_hat, grid, _F)
        psi_bt_hat, psi_bc_hat = baroclinic.to_modes(psi1_hat, psi2_hat)
        u_bt = grid.ifft(-grid.ddy(psi_bt_hat)); v_bt = grid.ifft(grid.ddx(psi_bt_hat))
        u_bc = grid.ifft(-grid.ddy(psi_bc_hat)); v_bc = grid.ifft(grid.ddx(psi_bc_hat))
        ke_bt = 0.5 * np.mean(u_bt ** 2 + v_bt ** 2)
        ke_bc = 0.5 * np.mean(u_bc ** 2 + v_bc ** 2)
        ape = 2.0 * _F * np.mean(grid.ifft(psi_bc_hat) ** 2)
        return ke_bt, ke_bc, ape

    _dt = 0.01
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 80)

    _snaps_bt, _snaps_bc, _t_hist, _KEbt, _KEbc, _APE = [], [], [], [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _q1h, _q2h = _state_hat
            _qbt_h, _qbc_h = baroclinic.to_modes(_q1h, _q2h)
            _snaps_bt.append(grid.ifft(_qbt_h).astype(np.float32))
            _snaps_bc.append(grid.ifft(_qbc_h).astype(np.float32))
            _t_hist.append(_t)
            _ke_bt, _ke_bc, _ape = _mode_energies(_state_hat)
            _KEbt.append(_ke_bt); _KEbc.append(_ke_bc); _APE.append(_ape)
        if _s < _nsteps:
            _state_hat = timestep.ifrk4_step(_state_hat, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    set_sim({
        "snaps_bt": _snaps_bt, "snaps_bc": _snaps_bc, "t": np.array(_t_hist),
        "KEbt": np.array(_KEbt), "KEbc": np.array(_KEbc), "APE": np.array(_APE),
        "grid": grid,
    })

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_snaps_bt)} frames."
    ), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the run has completed."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(mo, sim):
    frame_ui = mo.ui.slider(0, len(sim["t"]) - 1, step=1, value=len(sim["t"]) - 1,
                            label="frame", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, mo, plotting, plt, sim):
    _i = frame_ui.value
    fig1, axs1 = plt.subplots(1, 2, figsize=(9.5, 4.5), constrained_layout=True)
    plotting.field(axs1[0], sim["snaps_bt"][_i], grid=sim["grid"], signed=True,
                   title=f"$q_{{bt}}$, $t={sim['t'][_i]:.1f}$")
    plotting.field(axs1[1], sim["snaps_bc"][_i], grid=sim["grid"], signed=True,
                   title=f"$q_{{bc}}$, $t={sim['t'][_i]:.1f}$")
    mo.vstack([frame_ui, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **Left panel, frame 0:** should be exactly zero everywhere (flat
        color) — that's the point. Any nonzero $q_{bt}$ in later frames was
        *generated*, not initialized.
        """
    )
    return


@app.cell(hide_code=True)
def _(plt, sim):
    fig2, ax2 = plt.subplots(figsize=(7, 4.3), constrained_layout=True)
    ax2.plot(sim["t"], sim["KEbc"], lw=2, color="crimson", label="baroclinic KE")
    ax2.plot(sim["t"], sim["APE"], lw=2, color="#d97706", label="APE")
    ax2.plot(sim["t"], sim["KEbt"], lw=2, color="#2563eb", label="barotropic KE")
    ax2.set_xlabel("$t$"); ax2.set_ylabel("energy")
    ax2.set_title("mode energetics"); ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Confirm the leak.** The barotropic KE curve should start
          essentially at zero (floating-point noise) and grow measurably
          once the vortex begins to distort — that growth has no source
          except the nonlinear cross-terms `to_modes`'s docstring predicts.
          Compare its final magnitude to the baroclinic mode's: does it stay
          a small correction, or become comparable?
        - **Push $\beta$.** Increase $\beta$: the baroclinic vortex should
          develop a beta-gyre wake exactly like Ch. 8's single-layer vortex
          (it obeys the same equation, with deformation parameter $2F$) —
          confirm the westward drift.
        - **Vary $F$.** Larger $F$ means a *smaller* effective deformation
          radius for the baroclinic mode ($L_R^2\propto1/(2F)$): does the
          vortex's Rossby-wave wake wavelength shrink as you increase $F$,
          matching Ch. 8's own $F$-dependence?
        - **Zero the vortex, keep $\beta$.** Set the vortex amplitude to its
          minimum: with no baroclinic disturbance to begin with, does any
          barotropic energy appear at all? (It shouldn't — there is nothing
          for the Jacobian to couple.)

        ### What you should have seen

        A disturbance placed purely in the baroclinic mode does not stay
        purely baroclinic: nonlinear self-advection leaks a small but
        measurable barotropic circulation into existence, even with zero
        mean shear and zero instability. This is the *elementary* version
        of a process central to real ocean and atmosphere dynamics —
        baroclinic eddies (generated by Ch. 16's instability) are the
        primary source that feeds Ch. 18's barotropic-dominated inverse
        cascade, and this leak, run for long enough at large scale, is a
        large part of *why* that upscale transfer happens. Two-layer QG
        is the crudest possible discretization of Ch. 11's continuous
        vertical-mode theory — one interior mode, one deformation radius
        — but it already contains the essential coupling that theory's
        higher modes only refine.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
