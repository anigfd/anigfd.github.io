import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 15 — Barotropic Instability

        **Physical question.** Some shear flows spontaneously roll up into
        chains of vortices; others just sit there, however hard they're
        perturbed. After this chapter you should be able to state Rayleigh's
        and Kuo's necessary condition for instability, compute a
        growth-rate-vs-wavenumber curve for a given profile, and watch the
        fastest-growing mode roll up nonlinearly into vortices.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### From the vorticity equation to an eigenvalue problem

        Linearize Ch. 7's $Dq/Dt=0$ about a parallel flow $U(y)$: write
        $\psi=-\!\int U\,dy+\psi'$, keep terms linear in primes, and note
        the background carries the vorticity gradient
        $\partial_yq_0=\beta-U''$:

        $$\Big(\partial_t+U\partial_x\Big)\nabla^2\psi'
          +(\beta-U'')\,\psi'_x=0.$$

        The coefficients depend on $y$ only, so Fourier in $x$ and $t$ but
        not $y$: substituting a normal mode $\psi'=\phi(y)\,e^{ik(x-ct)}$
        (each factor $\partial_t\to-ikc$, $\partial_x\to ik$; divide through
        by $ik$) gives the **Rayleigh–Kuo equation** (symbols as in
        [NOTATION](../notation)):

        $$(U-c)(\phi''-k^2\phi)+(\beta-U'')\phi=0.$$

        For each wavenumber $k$ this is an eigenvalue problem for the
        complex phase speed $c=c_r+ic_i$; since
        $\psi'\propto e^{ikc_it}\,e^{ik(x-c_rt)}$, a mode is unstable iff
        $c_i>0$, with growth rate $kc_i$.

        ### Rayleigh's theorem: instability needs a sign change

        Multiply the equation by $\phi^*/(U-c)$, integrate across the
        channel, and take the imaginary part. The $|\phi'|^2+k^2|\phi|^2$
        term is real and drops out, leaving

        $$c_i\int\frac{(\beta-U'')\,|\phi|^2}{|U-c|^2}\,dy=0.$$

        If $c_i\neq0$ the integral itself must vanish — and since
        $|\phi|^2/|U-c|^2>0$, that is only possible if $\beta-U''$
        **changes sign** somewhere in the domain. That is Rayleigh's
        inflection-point criterion ($\beta=0$) and Kuo's extension
        ($\beta\neq0$): a *necessary* condition delivered by three lines of
        integration, with no eigenfunction ever computed. It is necessary,
        not sufficient — plenty of profiles with an inflection point are
        still stable at any given $k$, which is what Part A's solver is for.

        **The mechanism behind the criterion:** a sign change in
        $\beta-U''$ means the flow supports Rossby-type waves (Ch. 9) riding
        on *oppositely-signed* PV gradients on either flank. Each wave
        propagates counter to the local flow, so the pair can become
        stationary relative to each other, **phase-lock**, and mutually
        amplify — each wave's induced velocity field pushes the other's PV
        contour further from equilibrium. Every shear instability in
        Part V is a version of this two-wave resonance; Ch. 16 replays it
        with the two waves stacked vertically instead of side by side.
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
    from gfdlib import spectral, timestep, plotting, instability
    return instability, np, plotting, plt, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        **Growth-rate solver** (new numerics for this book): discretize
        $y$ on a channel with rigid walls, and rewrite the Rayleigh-Kuo
        equation as a generalized eigenvalue problem
        $A\phi=cB\phi$ — converted to the standard form $B^{-1}A$ (since
        $B=\partial_y^2-k^2I$ is invertible) and solved with plain
        `numpy.linalg.eigvals` (`gfdlib.instability.growth_rate`) — no SciPy,
        no generalized eigensolver required. **Nonlinear roll-up** reuses
        Ch. 18's exact machinery (`gfdlib.spectral.Grid`,
        `gfdlib.timestep.ifrk4_step`) with a new initial condition,
        `gfdlib.instability.double_shear_layer`.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part A — Growth-rate solver

        Fully reactive (a 1D eigenvalue problem is cheap) — no Run button.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    profile_ui = mo.ui.dropdown(
        options={
            "tanh shear layer (unstable)": "tanh",
            "sech² jet (unstable)": "jet",
            "cosine bump (stable — no inflection point)": "stable",
        },
        value="tanh shear layer (unstable)",
        label="base-state profile $U(y)$",
    )
    delta_a_ui = mo.ui.slider(0.5, 3.0, step=0.25, value=1.0,
                              label="width $\\delta$", show_value=True)
    beta_a_ui = mo.ui.slider(0.0, 2.0, step=0.1, value=0.0,
                             label="$\\beta$", show_value=True)
    mo.vstack([
        mo.hstack([profile_ui, delta_a_ui, beta_a_ui], justify="start"),
        mo.md(r"> The stable profile is $U=-\cos(\pi y/(2Y_{max}))$ on the "
              r"same-size channel — $U''\le0$ everywhere in the interior, "
              r"so Rayleigh's criterion rules out instability entirely."),
    ])
    return beta_a_ui, delta_a_ui, profile_ui


@app.cell(hide_code=True)
def _(beta_a_ui, delta_a_ui, instability, np, profile_ui):
    _delta = delta_a_ui.value
    _Ymax = 15 * _delta
    _y = np.linspace(-_Ymax, _Ymax, 402)[1:-1]

    if profile_ui.value == "tanh":
        _U = np.tanh(_y / _delta)
    elif profile_ui.value == "jet":
        _U = 1.0 / np.cosh(_y / _delta) ** 2
    else:
        _U = -np.cos(np.pi * _y / (2 * _Ymax))

    _k_values = np.linspace(0.02, 3.0 / _delta, 60)
    _growth = instability.growth_rate_curve(_y, _U, _k_values, beta=beta_a_ui.value)
    profile_result = {"y": _y, "U": _U, "k": _k_values, "growth": _growth}
    return (profile_result,)


@app.cell(hide_code=True)
def _(np, plt, profile_result):
    fig1, axs1 = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axs1[0].plot(profile_result["U"], profile_result["y"], lw=2)
    axs1[0].set_xlabel("$U(y)$"); axs1[0].set_ylabel("$y$")
    axs1[0].set_title("base-state profile"); axs1[0].grid(alpha=0.3)

    axs1[1].plot(profile_result["k"], profile_result["growth"], lw=2)
    axs1[1].set_xlabel("$k$"); axs1[1].set_ylabel("growth rate $k c_i$")
    axs1[1].set_title("growth-rate curve"); axs1[1].grid(alpha=0.3)
    if profile_result["growth"].max() > 1e-8:
        _kmax = profile_result["k"][np.argmax(profile_result["growth"])]
        axs1[1].axvline(_kmax, color="crimson", ls="--", lw=1,
                        label=f"peak at $k$={_kmax:.2f}")
        axs1[1].legend(fontsize=9)
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part B — Nonlinear roll-up

        A periodic **double shear layer** (two opposite-signed tanh jets, so
        the profile closes up on a doubly-periodic domain), seeded with a
        small perturbation at zonal wavenumber $n$. Use Part A's growth-rate
        curve to pick a wavenumber near the peak before running.

        Linear theory (Part A) can only tell you *what grows and how fast
        while it is small*. What it cannot tell you is what the flow becomes:
        exponential growth ends when the perturbation velocity is comparable
        to the shear itself, and then the vorticity of the shear layer
        wraps up into discrete vortices — the nonlinear saturation. Watch
        for three stages in the snapshot browser: exponential amplification
        of the seeded wiggle (invisible at first), roll-up into cat's-eye
        vortices, then vortex merging and filamentation (Ch. 18's
        turbulence, arriving on schedule).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    delta_b_ui = mo.ui.slider(0.05, 0.3, step=0.01, value=0.16,
                              label="shear-layer width $\\delta$", show_value=True)
    npert_ui = mo.ui.slider(1, 8, step=1, value=3,
                            label="seeded wavenumber $n$", show_value=True)
    vpert_ui = mo.ui.slider(0.01, 0.1, step=0.01, value=0.05,
                            label="perturbation amplitude", show_value=True)
    beta_b_ui = mo.ui.slider(0.0, 10.0, step=1.0, value=0.0,
                             label="$\\beta$", show_value=True)
    n_ui = mo.ui.dropdown(options={"64": 64, "128": 128, "256": 256}, value="128",
                          label="resolution $n$")
    T_ui = mo.ui.slider(10.0, 40.0, step=5.0, value=30.0,
                        label="run time $T$", show_value=True)
    mo.vstack([
        mo.hstack([delta_b_ui, npert_ui, vpert_ui], justify="start"),
        mo.hstack([beta_b_ui, n_ui, T_ui], justify="start"),
    ])
    return T_ui, beta_b_ui, delta_b_ui, n_ui, npert_ui, vpert_ui


@app.cell(hide_code=True)
def _(mo):
    get_sim, set_sim = mo.state(None)   # results survive the run button resetting
    return get_sim, set_sim


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Run roll-up")
    run_btn
    return (run_btn,)


@app.cell
def _(
    T_ui, beta_b_ui, delta_b_ui, instability, mo, n_ui, np, npert_ui, run_btn,
    set_sim, spectral, timestep, vpert_ui,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run roll-up** above to start."), kind="info"
    ))

    _n, _delta, _npert, _vpert, _beta, _T = (
        n_ui.value, delta_b_ui.value, npert_ui.value, vpert_ui.value,
        beta_b_ui.value, T_ui.value,
    )
    grid = spectral.Grid(_n)
    _zeta0 = instability.double_shear_layer(grid, _delta, v_pert=_vpert, n_pert=_npert)

    _nu, _nnu = 1e-5, 2
    _zeta_hat = grid.fft(_zeta0) * grid.dealias
    _L_op = -_nu * grid.k2 ** _nnu + 1j * _beta * grid.kx * grid.k2_inv

    def _rhs_nl(t, zh):
        return -grid.jacobian(grid.invert_laplacian(zh), zh)

    _dt = 0.005
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 80)

    _snaps, _t_hist = [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _snaps.append(grid.ifft(_zeta_hat).astype(np.float32))
            _t_hist.append(_t)
        if _s < _nsteps:
            _zeta_hat = timestep.ifrk4_step(_zeta_hat, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    set_sim({"snaps": _snaps, "t": np.array(_t_hist), "grid": grid})

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_snaps)} frames."
    ), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the roll-up has run."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(mo, sim):
    frame_ui = mo.ui.slider(0, len(sim["snaps"]) - 1, step=1, value=len(sim["snaps"]) - 1,
                            label="frame", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, mo, plotting, plt, sim):
    fig2, ax2 = plt.subplots(figsize=(6, 5.5), constrained_layout=True)
    _i = frame_ui.value
    plotting.field(ax2, sim["snaps"][_i], grid=sim["grid"], signed=True,
                   title=f"$\\zeta$ at $t={sim['t'][_i]:.1f}$")
    mo.vstack([frame_ui, fig2])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Confirm the criterion.** In Part A, switch to the "stable"
          profile: the growth-rate curve should be exactly zero everywhere,
          for any $\delta$ or $\beta$. Rayleigh's criterion isn't just
          necessary here — it's decisive: with $U''$ one-signed in the
          interior, there are no counter-propagating waves to lock together.
        - **Miss the peak.** In Part B, set the seeded wavenumber $n$ well
          away from the peak $k$ you found in Part A (converting via
          $k=2\pi n/L$): does the roll-up still happen, and does it take
          longer? Watch *which* wavenumber actually wins — the seeded mode
          grows first, but background noise at the fastest-growing
          wavenumber is also being amplified the whole time, and given a
          long enough run it takes over. Nature seeds all wavenumbers at
          once; the growth-rate curve is a *selection* principle.
        - **Stabilize with $\beta$.** In Part A, increase $\beta$ for the
          jet profile: does the growth-rate curve's peak shrink, grow, or
          shift? $\beta>0$ everywhere pushes $\beta-U''$ toward one sign —
          for a strong enough $\beta$ the sign change (and with it the
          instability) disappears entirely. This is why planetary rotation
          gradients *stabilize* jets, and part of why zonal jets are the
          natural end state on a $\beta$-plane (Ch. 18).
        - **Count the vortices.** In Part B's final frame, count the
          vortex cores in one shear layer — it should match your chosen $n$
          (unless mergers have already begun; then count the earlier frames).

        ### What you should have seen

        A shear flow's fate is fixed by the curvature of its own velocity
        profile: no inflection point in $\beta-U''$ means guaranteed
        stability, confirmed by a growth-rate curve that's identically zero.
        Where an inflection point *does* exist, instability grows fastest at
        a specific wavenumber — and seeding exactly that wavenumber in the
        nonlinear simulation produces a clean chain of Kelvin-Helmholtz
        vortices, one per wavelength. This is the same mechanism (in a
        rotating, stratified form) behind meanders and eddy shedding in
        jet streams and western boundary currents.

        **Where this goes next.** This chapter's script — linearize, find a
        necessary criterion by an integral theorem, solve the eigenproblem
        for growth rates, then watch nonlinearity saturate the winner — is
        repeated verbatim in Ch. 16 (baroclinic instability: same theorem
        structure, with the PV gradient's sign change now in the *vertical*)
        and Ch. 17 (stratified and symmetric instabilities). Keep the
        two-counter-propagating-waves cartoon; it is the one mental model
        that survives every generalization.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
