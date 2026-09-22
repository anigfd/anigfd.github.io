import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 19 — Eddy Transport & Mixing

        **Physical question.** Ch. 18's turbulent eddies don't just cascade
        energy and enstrophy — they also stir anything else riding along
        with the flow: heat, salt, chemical tracers, potential vorticity
        itself. After this chapter you should be able to explain what an
        "eddy diffusivity" actually measures, and compute one directly from
        a turbulent simulation rather than taking it as a free parameter.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### Why "eddy diffusivity" needs defending

        Molecular diffusion moves heat down its gradient because molecules
        genuinely random-walk. Eddies are not molecules — they are ordered,
        long-lived, and much bigger than the gradients they act on — so the
        claim that their net effect on a tracer looks like a (much larger)
        diffusivity is a *hypothesis*, one that every ocean and climate
        model bets on when it writes
        $\overline{v'c'}=-K\,\partial\overline C/\partial y$ for eddies it
        cannot resolve. This notebook's purpose is to test that bet in the
        one setting where the answer can be *measured* cleanly.

        ### The mean-gradient trick

        Split the tracer into an imposed, unbounded **mean gradient** and a
        doubly-periodic **perturbation**:
        $C_{total}(x,y,t)=\Gamma y+c'(x,y,t)$ (symbols as in
        [NOTATION](../notation)). Substituting into the advection–diffusion
        equation, $c'$ obeys

        $$\frac{\partial c'}{\partial t}+J(\psi,c')=-\Gamma v+\kappa\nabla^2c',$$

        stirred by the SAME vorticity field $\zeta=\nabla^2\psi$ evolving
        under its own barotropic dynamics (Ch. 18) — the tracer is passive:
        advected by the flow, but exerting no force back on it. The
        $-\Gamma v$ source term has a plain physical reading: a parcel
        moving *up* the mean gradient ($v>0$) arrives carrying less tracer
        than its new surroundings, i.e. a negative $c'$ — and that
        systematic correlation between $v'$ and $c'$ **is** the eddy flux.
        The turbulent flux defines an **effective diffusivity** through the
        flux–gradient closure:

        $$K_{eff}=-\frac{\overline{v'c'}}{\Gamma}.$$

        Dimensionally $K_{eff}\sim u_{rms}\,\ell_{mix}$ — velocity times the
        distance a parcel travels before its tracer identity is blended
        away (Prandtl's *mixing length*). For ocean mesoscale eddies
        ($u\sim0.1$ m/s, $\ell\sim50$ km) that gives
        $K_{eff}\sim10^3$–$10^4\,$m²/s, against a molecular
        $\kappa\sim10^{-7}$ m²/s for heat: stirring beats diffusion by ten
        orders of magnitude, which is why the parameterization question is
        existential for climate modeling.
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
    from gfdlib import spectral, timestep, diagnostics, plotting, mixing
    return diagnostics, mixing, np, plotting, plt, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        `gfdlib.mixing.rhs_coupled` stacks vorticity and tracer-perturbation
        into one state (exactly ch16's two-layer pattern), reusing ch18's
        exact turbulence machinery (`spectral.Grid`, McWilliams initial
        condition, hyperviscosity) for $\zeta$ unchanged, with the mean
        gradient folded in as an ordinary source term $-\Gamma v$ for $c'$
        — the same "split into mean + periodic perturbation" trick used for
        $\beta y$ (ch. 18) and mean shear (ch. 16). Because the $c'$
        equation is *linear*, $K_{eff}$ must be exactly independent of
        $\Gamma$ for a fixed flow — verified numerically, not assumed (see
        `gfdlib` tests and "Try this" below).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    n_ui = mo.ui.dropdown(options={"64": 64, "128": 128}, value="128", label="resolution $n$")
    k0_ui = mo.ui.slider(4, 16, step=1, value=8, label="initial peak $k_0$", show_value=True)
    nu_ui = mo.ui.slider(0.0, 1e-6, step=1e-8, value=1e-7, label="hyperviscosity $\\nu$",
                         show_value=True)
    Gamma_ui = mo.ui.slider(0.2, 3.0, step=0.2, value=1.0, label="mean gradient $\\Gamma$",
                            show_value=True)
    kappa_ui = mo.ui.slider(0.0, 1e-6, step=1e-8, value=1e-7,
                            label="molecular diffusivity $\\kappa$", show_value=True)
    nsteps_ui = mo.ui.slider(500, 3000, step=250, value=2000, label="steps", show_value=True)
    dt_ui = mo.ui.slider(0.005, 0.02, step=0.0025, value=0.01, label="$dt$", show_value=True)
    mo.vstack([
        mo.hstack([n_ui, k0_ui, nu_ui], justify="start"),
        mo.hstack([Gamma_ui, kappa_ui], justify="start"),
        mo.hstack([nsteps_ui, dt_ui], justify="start"),
    ])
    return Gamma_ui, dt_ui, k0_ui, kappa_ui, n_ui, nsteps_ui, nu_ui


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
    Gamma_ui, diagnostics, dt_ui, k0_ui, kappa_ui, mixing, mo, n_ui, np,
    nsteps_ui, nu_ui, run_btn, set_sim, spectral, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run simulation** above to start."), kind="info"
    ))

    _n, _dt, _nsteps = n_ui.value, dt_ui.value, nsteps_ui.value
    _k0, _nu, _kappa, _Gamma = k0_ui.value, nu_ui.value, kappa_ui.value, Gamma_ui.value
    _nnu = 2
    _E0 = 0.5
    _nsub = max(1, _nsteps // 100)

    grid = spectral.Grid(_n)

    # McWilliams (1984) initial condition, identical to ch18
    _rng = np.random.default_rng(1234)
    _K = np.where(grid.kmag == 0, 1e-10, grid.kmag)
    _amp = _K * np.sqrt(1.0 / (1.0 + (_K / _k0) ** 4))
    _zeta_hat = grid.dealias * _amp * np.exp(2j * np.pi * _rng.random(_K.shape))
    _ke, _ = diagnostics.energy_enstrophy(grid.invert_laplacian(_zeta_hat), grid)
    _zeta_hat = _zeta_hat * np.sqrt(_E0 / _ke)

    state_hat = np.stack([_zeta_hat, np.zeros_like(_zeta_hat)])
    _L_op = np.stack([-_nu * grid.k2 ** _nnu, -_kappa * grid.k2 ** _nnu])

    def _rhs_nl(t, sh):
        return mixing.rhs_coupled(sh, grid, _Gamma)

    _zeta_snaps, _c_snaps, _t_hist, _Keff_hist = [], [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _zh, _ch = state_hat
            _psi_hat = grid.invert_laplacian(_zh)
            _v = grid.ifft(grid.ddx(_psi_hat))
            _cprime = grid.ifft(_ch)
            _zeta_snaps.append(grid.ifft(_zh).astype(np.float32))
            _c_snaps.append(_cprime.astype(np.float32))
            _t_hist.append(_t)
            _Keff_hist.append(mixing.effective_diffusivity(_v, _cprime, _Gamma)
                              if _t > 0 else 0.0)
        if _s < _nsteps:
            state_hat = timestep.ifrk4_step(state_hat, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    set_sim({
        "zeta": _zeta_snaps, "c": _c_snaps, "t": np.array(_t_hist),
        "Keff": np.array(_Keff_hist), "Gamma": _Gamma, "grid": grid,
    })

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_t_hist)} frames."
    ), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the simulation has run."), kind="warn"
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
    plotting.field(axs1[0], sim["zeta"][_i], grid=sim["grid"], signed=True,
                   title=f"vorticity $\\zeta$, $t={sim['t'][_i]:.1f}$")
    plotting.field(axs1[1], sim["c"][_i], grid=sim["grid"], signed=True,
                   title=f"tracer perturbation $c'$, $t={sim['t'][_i]:.1f}$")
    mo.vstack([frame_ui, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### The effective diffusivity, measured directly

        $K_{eff}=-\overline{v'c'}/\Gamma$, computed at every saved frame.
        Watch it spin up from zero as filaments of tracer develop, then
        fluctuate around a well-defined magnitude that drifts only slowly
        compared to its initial rise — this is decaying turbulence, so the
        eddies doing the stirring are themselves still slowly evolving
        (merging into fewer, larger vortices), not a numerical artifact. A
        real, emergent turbulent diffusivity, not a parameter that was
        dialed in.
        """
    )
    return


@app.cell(hide_code=True)
def _(plt, sim):
    fig2, ax2 = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    ax2.plot(sim["t"], sim["Keff"], lw=2)
    ax2.axhline(0, color="k", lw=0.5)
    ax2.set_xlabel("$t$"); ax2.set_ylabel("$K_{eff}(t)$")
    ax2.set_title(f"effective diffusivity ($\\Gamma$={sim['Gamma']:.1f})")
    ax2.grid(alpha=0.3)
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Confirm $\Gamma$-independence.** Run once at $\Gamma=1$, note the
          late-time $K_{eff}$, then run again at $\Gamma=2$ with every other
          slider unchanged. Because the tracer equation is linear in $c'$,
          $K_{eff}$ should come out the same (to a fraction of a percent) —
          $K_{eff}$ is a property of the *flow*, not of how strong a
          gradient you chose to probe it with.
        - **Turn off the flow.** Set $\nu$ and $k_0$ however you like but
          imagine $\zeta\equiv0$: with no velocity there can be no $v'c'$
          correlation, so $K_{eff}\to0$ — confirmed in `gfdlib`'s test
          suite; you can see the same effect approximately by comparing
          early frames (before the flow has spun up any real velocity) to
          late ones.
        - **Compare to a mixing-length estimate.** Estimate
          $u_{rms}=\sqrt{2E}$ from the initial energy $E_0=0.5$ and the
          domain size $L=2\pi$; check that the late-time $K_{eff}$ you
          measure is the right order of magnitude for $u_{rms}\times L$.

        ### What you should have seen

        The tracer perturbation $c'$ develops long, thin filaments wrapped
        around the same eddies visible in $\zeta$ — turbulent stirring is
        literally the mechanism, visible directly in the field. $K_{eff}(t)$
        rises from zero as those filaments develop and then fluctuates
        around a clearly positive (down-gradient), well-defined magnitude —
        tracking the turbulence's own continued slow evolution (vortex
        mergers) rather than settling to a perfectly flat plateau, since
        this is freely-decaying, not statistically-steady, turbulence. A
        genuine eddy diffusivity, orders of magnitude larger than any
        molecular $\kappa$, and exactly the quantity ocean and atmosphere
        models parameterize when they can't afford to resolve the eddies
        themselves.

        **A caution worth carrying forward:** the flux-gradient closure
        worked here partly because the tracer was passive and the gradient
        imposed. For tracers that *feel back* on the flow — PV above all —
        eddy fluxes can be spatially inhomogeneous (near-zero inside
        Ch. 7's staircase risers, huge in the mixed treads) and even
        locally *up*-gradient. "The eddies act like diffusion" is a good
        first model and a famously dangerous last one.

        **Where this goes next.** Down-gradient PV mixing by exactly these
        eddies is what sharpens Ch. 18's jets; the diffusivity you measured
        is the quantity the Gent–McWilliams parameterization supplies to
        every non-eddy-resolving ocean climate model; and the microscale
        end of the same story — how the filaments' variance is finally
        destroyed — is Ch. 17's Kelvin–Helmholtz billows, feeding the
        abyssal $\kappa$ that Ch. 21's overturning balance runs on.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
