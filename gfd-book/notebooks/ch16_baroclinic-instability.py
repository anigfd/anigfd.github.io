import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 16 — Baroclinic Instability

        **Physical question.** Ch. 15's shear instability draws energy from
        horizontal shear of the mean flow's *kinetic* energy. Mid-latitude
        storms and ocean eddies draw energy from somewhere else entirely: the
        *potential* energy stored in a sloping temperature surface — the
        same thermal-wind shear from Ch. 5. After this chapter you should be
        able to compute the Eady growth rate for a given stratification and
        shear, and watch a baroclinic wave grow, break, and equilibrate in a
        two-layer model — the single most important instability in this
        book.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The energy source: available potential energy

        Ch. 5 showed that a meridional temperature gradient in thermal-wind
        balance means *sloping* buoyancy surfaces. Sloping surfaces store
        potential energy that flat ones don't:

        > **Available potential energy (APE)** — the fraction of a fluid's
        > potential energy that could be released by adiabatically
        > flattening its buoyancy surfaces. For the midlatitude atmosphere
        > it exceeds the kinetic energy of the winds by an order of
        > magnitude. Baroclinic instability is the mechanism that taps it.

        The release path is geometric — the **wedge of instability**: a
        parcel exchanged along a path *shallower than the buoyancy slope but
        steeper than horizontal* moves cold fluid down and warm fluid up,
        lowering the center of mass and freeing energy, all while staying
        (quasi-)balanced. Trajectories inside the wedge tilt against the
        shear — which is why growing baroclinic waves lean westward with
        height, and why that tilt is the observable signature of a storm
        still deepening.

        ### The Eady problem

        The cleanest quantitative model (Eady 1949): QG flow between rigid
        lids at $z=0,H$, uniform shear $U=\Lambda z$, constant $N$, no
        $\beta$ — chosen so the *interior* PV gradient is exactly **zero**.
        All the dynamics then lives on the two boundaries, where the
        advected temperature acts as a PV sheet: a warm anomaly on the
        ground behaves like a Rossby wave running one way, a temperature
        anomaly on the lid like one running the other way. Ch. 15's
        two-counter-propagating-waves resonance again — with the two waves
        now stacked *vertically* and coupled across the depth $H$.

        Normal modes $\psi'=\phi(z)e^{ik(x-ct)}$ satisfy
        $\phi''-k^2\phi=0$ in the interior (zero PV there!), matched to the
        boundary temperature equations $(U-c)\phi_z-\Lambda\phi=0$ at
        $z=0,H$ — a $2\times2$ eigenvalue problem for $c$, with growth rate
        $\sigma=\mu c_i$, $\mu=kNH/f_0$ (symbols as in
        [NOTATION](../notation)). Two robust predictions fall out:

        - a **short-wave cutoff** at $\mu\approx2.4$: waves much narrower
          than the deformation radius can't couple the two boundaries, so
          the resonance dies — precisely why the deformation radius
          $L_d=NH/f_0$ sets the size of storms and eddies;
        - a most-unstable scale $\mu\approx1.6$ with growth rate
          $\sigma\approx0.31\,\Lambda f_0/N$ — the "Eady timescale" used to
          map storm-track and eddy activity to this day.

        ### The 2-layer (Phillips) model

        The nonlinear, finite-amplitude analog: two QG layers with mean
        flows $U_1,U_2$ and PV

        $$q_1=\nabla^2\psi_1+F(\psi_2-\psi_1)+\beta y,\qquad
          q_2=\nabla^2\psi_2+F(\psi_1-\psi_2)+\beta y,$$

        with $F=f_0^2/(g'H)$ (the inverse deformation radius squared — the
        vertical coupling constant). The mean state hands each layer an
        *effective* PV gradient: $\beta_1=\beta+F\Delta U$,
        $\beta_2=\beta-F\Delta U$ ($\Delta U=U_1-U_2$). The **Phillips
        necessary condition** is that these can have opposite signs —
        i.e. $F\Delta U>\beta$ — a direct vertical analog of Ch. 15's
        Rayleigh–Kuo sign-change criterion, with the shear now supplying
        the negative gradient in the lower layer. Note what $\beta$ does
        here: it is purely *stabilizing*, setting a minimum shear below
        which no instability exists at all.
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
    from gfdlib import spectral, timestep, plotting, baroclinic
    return baroclinic, np, plotting, plt, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Both growth-rate problems are small ($2\times2$) eigenvalue
        problems, derived directly rather than taken from a memorized
        closed form, and solved with plain `numpy.linalg.eigvals` — no
        SciPy. The nonlinear life cycle reuses `gfdlib.spectral.Grid` and
        `gfdlib.timestep.ifrk4_step`, with `gfdlib.baroclinic.invert_2layer`
        (a closed-form $2\times2$ Helmholtz solve) providing $\psi_1,\psi_2$
        and hyperviscosity as the exact linear operator. A fixed mean shear
        with no other sink is an unlimited energy source, so
        `gfdlib.baroclinic.rhs_2layer` includes bottom Ekman drag on the
        lower layer — without it, growth continues past any physically
        reasonable amplitude regardless of hyperviscosity; with it, the
        classic growth-peak-decay-equilibrate life cycle emerges. A full
        life cycle needs ~$10^5$ time steps at the notebook's resolution —
        too heavy to integrate live in the browser (see performance notes in
        `CLAUDE.md`), so Part B's canonical run is precomputed offline and
        shipped as data for the widget below to scrub; a second, smaller
        panel lets you run a short, fully live simulation with your own
        parameters to test hypotheses about the growth phase.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"## Part A — The Eady growth-rate calculator")
    return


@app.cell(hide_code=True)
def _(baroclinic, np, plt):
    _mu = np.linspace(0.02, 3.2, 300)
    _growth = baroclinic.eady_growth_rate_curve(_mu)
    _mu_max = _mu[np.argmax(_growth)]

    fig0, ax0 = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax0.plot(_mu, _growth, lw=2)
    ax0.axvline(_mu_max, color="crimson", ls="--", lw=1,
               label=f"peak at $\\mu$={_mu_max:.2f}")
    ax0.set_xlabel("$\\mu=kNH/f_0$"); ax0.set_ylabel("growth rate $\\sigma=\\mu c_i$")
    ax0.set_title("Eady growth-rate curve"); ax0.legend(fontsize=9); ax0.grid(alpha=0.3)
    fig0
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Both fluids, same curve, different scales

        The nondimensional curve above is universal; converting
        $\mu_{max}\approx1.61$ to a physical wavelength
        $L=2\pi/(\mu_{max}/(N H/f_0))$ and e-folding time
        $\tau=1/\sigma_{max}$ for two very different settings shows why
        storms and ocean eddies look nothing alike, despite obeying the
        same equation.
        """
    )
    return


@app.cell(hide_code=True)
def _(baroclinic, mo, np):
    _mu_max = 1.606
    _sigma_max_nondim = 0.310

    def _report(name, N, f0, H, Lambda):
        Ld = N * H / f0
        k_max = _mu_max / Ld
        L_wave = 2 * np.pi / k_max
        sigma_dim = _sigma_max_nondim * Lambda * f0 / N   # sigma_nondim * (Lambda*f0/N)
        tau = 1.0 / sigma_dim
        return f"**{name}:** $L_d$={Ld/1e3:.0f} km, most-unstable wavelength $\\approx${L_wave/1e3:.0f} km, e-folding time $\\approx${tau/3600:.1f} hours"

    _atmos = _report("Atmosphere (storm track)", N=0.01, f0=1e-4, H=1e4, Lambda=3e-3)
    _ocean = _report("Ocean (mesoscale eddies)", N=0.005, f0=1e-4, H=1e3, Lambda=2e-3)
    mo.md(_atmos + "\n\n" + _ocean +
          "\n\n(Illustrative parameter choices, not a fit to any specific "
          "observed storm or eddy — but the ~10x scale separation and "
          "~10x timescale separation are the real, robust story.)")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part B — 2-layer nonlinear life cycle

        ### B1. The full life cycle (precomputed)

        This is the run PLAN.md's "life-cycle run" refers to: a
        small-amplitude perturbation seeded on a mean shear unstable by the
        Phillips criterion, integrated with `gfdlib.baroclinic.rhs_2layer`
        at the notebook's full $128\times128$ resolution out to $t=500$ —
        about $10^5$ time steps, precomputed offline (see the numerical-scheme
        note above) and shipped as data for the sliders below to scrub.
        """
    )
    return


@app.cell
async def _(mo, np):
    import sys as _sys

    if _sys.platform == "emscripten":
        import pyodide.http
        import io
        _resp = await pyodide.http.pyfetch(
            str(mo.notebook_location() / "public" / "ch16_baroclinic-instability.npz")
        )
        lifecycle = np.load(io.BytesIO(await _resp.bytes()))
    else:
        lifecycle = np.load(mo.notebook_dir() / "public" / "ch16_baroclinic-instability.npz")
    return (lifecycle,)


@app.cell(hide_code=True)
def _(lifecycle, mo):
    mo.md(
        f"**Precomputed run parameters:** $n$={int(lifecycle['n'])}, "
        f"$F$={float(lifecycle['F']):.2f}, $U_1$={float(lifecycle['U1']):.2f}, "
        f"$U_2$={float(lifecycle['U2']):.2f}, $\\beta$={float(lifecycle['beta']):.2f}, "
        f"Ekman drag $r$={float(lifecycle['r']):.2f}, seeded at wavenumber "
        f"$k$={int(lifecycle['kseed'])} — the same parameters used to derive "
        f"the life cycle described in `NOTATION.md`."
    )
    return


@app.cell(hide_code=True)
def _(lifecycle, mo):
    frame_ui = mo.ui.slider(0, len(lifecycle["t"]) - 1, step=1,
                            value=len(lifecycle["t"]) - 1, label="frame", show_value=True)
    frame_ui
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, lifecycle, np, plotting, plt):
    _i = frame_ui.value

    class _G:  # minimal shim: plotting.field only needs the domain length L
        L = 2 * np.pi

    fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(9.5, 4.5), constrained_layout=True)
    plotting.field(ax1a, lifecycle["q1"][_i], grid=_G(), signed=True,
                   title=f"upper layer $q_1'$, $t={lifecycle['t'][_i]:.0f}$")
    plotting.field(ax1b, lifecycle["q2"][_i], grid=_G(), signed=True,
                   title=f"lower layer $q_2'$, $t={lifecycle['t'][_i]:.0f}$")
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Energetics: the life cycle in one plot

        Eddy kinetic energy (barotropic-like) and eddy potential energy
        (proportional to layer-interface displacement squared) both grow
        exponentially, peak, and decay toward a statistical equilibrium.

        **How to read it:** during the linear phase both curves rise on
        parallel straight lines (log scale: same exponential rate — the
        growing normal mode has a fixed EPE/EKE partition). The energy
        *pathway* is mean APE $\to$ eddy PE $\to$ eddy KE: the wave first
        distorts the interface (extracting potential energy from the mean
        slope), then converts that displacement into swirling motion. The
        peak and decay is the wave breaking and the drag draining what the
        instability delivered; the mean shear — held fixed here — keeps
        feeding it, which is why the end state is a statistical equilibrium
        rather than rest. In the real atmosphere the same arc, growth to
        breaking, takes about a week and is called a *storm*.
        """
    )
    return


@app.cell(hide_code=True)
def _(lifecycle, plt):
    fig2, ax2 = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    ax2.plot(lifecycle["t"], lifecycle["EKE"], lw=2, label="EKE")
    ax2.plot(lifecycle["t"], lifecycle["EPE"], lw=2, color="crimson", label="EPE")
    ax2.plot(lifecycle["t"], lifecycle["EKE"] + lifecycle["EPE"], lw=1.5, ls="--",
             color="gray", label="total")
    ax2.set_yscale("log")
    ax2.set_xlabel("$t$"); ax2.set_ylabel("energy (log scale)")
    ax2.set_title("eddy energetics"); ax2.legend(fontsize=9); ax2.grid(alpha=0.3, which="both")
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### B2. Explore the mechanism yourself (live)

        A smaller, fully live version of the same model — low enough
        resolution and short enough run time to stay interactive in the
        browser — so you can test the hypotheses in "Try this" below with
        your own parameters. This panel only reaches the *growth phase*; for
        the full peak-decay-equilibration arc see B1 above.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    F_ui = mo.ui.slider(0.5, 2.0, step=0.1, value=1.0, label="$F$", show_value=True)
    dU_ui = mo.ui.slider(0.2, 1.2, step=0.1, value=0.6, label="$\\Delta U$", show_value=True)
    beta_ui = mo.ui.slider(0.0, 1.0, step=0.1, value=0.5, label="$\\beta$", show_value=True)
    r_ui = mo.ui.slider(0.0, 0.8, step=0.05, value=0.4, label="Ekman drag $r$", show_value=True)
    T_ui = mo.ui.slider(20.0, 100.0, step=10.0, value=60.0, label="run time $T$", show_value=True)
    mo.vstack([
        mo.hstack([F_ui, dU_ui, beta_ui], justify="start"),
        mo.hstack([r_ui, T_ui], justify="start"),
        mo.md(r"> Fixed at $n=48$ to stay interactive; this panel shows "
              r"growth/onset, not the full multi-hundred-time-unit life "
              r"cycle (that's what B1's precomputed run is for)."),
    ])
    return F_ui, T_ui, beta_ui, dU_ui, r_ui


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
    F_ui, T_ui, baroclinic, beta_ui, dU_ui, mo, np, r_ui, run_btn,
    set_sim, spectral, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run** above to start."), kind="info"
    ))

    _n = 48
    _F, _dU, _beta, _r, _T = F_ui.value, dU_ui.value, beta_ui.value, r_ui.value, T_ui.value
    _U1, _U2 = 0.5 * _dU, -0.5 * _dU
    grid = spectral.Grid(_n)

    _nu, _nnu = 5e-3, 4
    _kseed = 1
    _rng = np.random.default_rng(5)
    _eps = 0.01
    _q1 = _eps * np.cos(_kseed * grid.x) + 0.001 * _rng.standard_normal((_n, _n))
    _q2 = -_eps * np.cos(_kseed * grid.x) + 0.001 * _rng.standard_normal((_n, _n))
    _state_hat = grid.fft(np.stack([_q1, _q2])) * grid.dealias
    _L_op = -_nu * grid.k2 ** _nnu

    def _rhs_nl(t, sh):
        return baroclinic.rhs_2layer(sh, grid, _F, _U1, _U2, _beta, r=_r)

    def _energetics(state_hat):
        q1_hat, q2_hat = state_hat
        psi1_hat, psi2_hat = baroclinic.invert_2layer(q1_hat, q2_hat, grid, _F)
        u1 = grid.ifft(-grid.ddy(psi1_hat)); v1 = grid.ifft(grid.ddx(psi1_hat))
        u2 = grid.ifft(-grid.ddy(psi2_hat)); v2 = grid.ifft(grid.ddx(psi2_hat))
        eke = 0.5 * np.mean(u1 ** 2 + v1 ** 2 + u2 ** 2 + v2 ** 2)
        psi1, psi2 = grid.ifft(psi1_hat), grid.ifft(psi2_hat)
        epe = 0.5 * _F * np.mean((psi1 - psi2) ** 2)
        return eke, epe

    _dt = 0.005
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 60)

    _t_hist, _EKE, _EPE = [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _eke, _epe = _energetics(_state_hat)
            _EKE.append(_eke); _EPE.append(_epe); _t_hist.append(_t)
        if _s < _nsteps:
            _state_hat = timestep.ifrk4_step(_state_hat, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    set_sim({"t": np.array(_t_hist), "EKE": np.array(_EKE), "EPE": np.array(_EPE)})

    mo.callout(mo.md(f"**Done.** {_nsteps} steps on a {_n}×{_n} grid."), kind="success")
    return


@app.cell(hide_code=True)
def _(get_sim, mo, plt):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once you press Run."), kind="warn"
    ))

    fig3, ax3 = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    ax3.plot(sim["t"], sim["EKE"], lw=2, label="EKE")
    ax3.plot(sim["t"], sim["EPE"], lw=2, color="crimson", label="EPE")
    ax3.set_yscale("log")
    ax3.set_xlabel("$t$"); ax3.set_ylabel("energy (log scale)")
    ax3.set_title("your run's eddy energetics"); ax3.legend(fontsize=9)
    ax3.grid(alpha=0.3, which="both")
    fig3
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Kill the drag.** Set Ekman drag $r$ to $0$ in B2: the growth
          rate is unchanged (drag doesn't affect the *linear* instability),
          but nothing will ever cap it — B1's precomputed run shows what
          happens if you let that continue: without drag, energy grows past
          any physically reasonable amplitude, regardless of hyperviscosity.
          The general lesson: *what limits an instability is almost never
          what starts it.*
        - **Match Part A.** Using $F$ and $\Delta U$, estimate an
          Eady-equivalent $\Lambda\sim\Delta U/H$ and $N^2\sim f_0^2/F$,
          then check whether the seeded wavenumber ($k=1$) is close to
          Part A's $\mu_{max}$ for those parameters.
        - **Push $\beta$.** Increase $\beta$ in B2 until the flow no longer
          goes unstable (Phillips' necessary condition:
          $\beta_2=\beta-F\Delta U$ must be able to change sign, so the
          threshold is at $\beta=F\Delta U$ — compute it from your slider
          values *before* running) — find it experimentally by watching EKE
          stay flat instead of growing. This threshold is real physics:
          it's why the ocean's weakly-sheared interior is only marginally
          baroclinically unstable, while strongly-sheared western boundary
          current extensions are eddy factories.
        - **Watch the tilt.** In B1, pick a frame during the growth phase
          and compare the $q_1'$ and $q_2'$ patterns: the upper-layer wave
          should sit shifted *westward* (leftward) of the lower-layer wave.
          That phase tilt against the shear is the wedge-of-instability
          geometry made visible — and when the tilt vanishes near the
          energy peak, growth stops. Forecasters look for exactly this tilt
          in real soundings.

        ### What you should have seen

        Both halves of this notebook tell the same story at different
        levels: the Eady calculator predicts a most-unstable scale and an
        e-folding time from stratification and shear alone — length scales
        an order of magnitude apart for the atmosphere and ocean, exactly as
        observed — and the 2-layer model shows what that instability
        actually *does*: grow exponentially, extract available potential
        energy from the sloping mean state, break nonlinearly, and settle
        into a statistically steady turbulent equilibrium. This is the
        mechanism that generates essentially every mid-latitude storm and
        every ocean mesoscale eddy.

        **Where this goes next.** The equilibrated end state of this
        notebook — a soup of eddies stirred by an inexhaustible mean
        gradient — is the *starting point* of Part VI: Ch. 18 studies the
        eddy soup's own dynamics (cascades, jets), Ch. 19 asks what the
        eddies *transport*, and Chs. 20–21 build the mean circulations they
        feed on. Ch. 17 first finishes the instability survey with the
        cases this chapter's balanced framework can't see.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
