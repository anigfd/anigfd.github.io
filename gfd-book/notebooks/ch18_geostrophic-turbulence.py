import marimo

__generated_with = "0.19.10"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 18 — Geostrophic (2D) Turbulence

        **Physical question.** Two-dimensional flow conserves *two* quadratic
        invariants — energy and enstrophy — and that second constraint reverses
        the direction of the energy cascade: energy flows **upscale**, small
        vortices merge into large ones. On a $\beta$-plane the upscale cascade
        runs into Rossby waves and is arrested into **zonal jets**. After this
        chapter you should be able to predict *when* vortices merge, *what*
        spectrum the enstrophy cascade produces, and *at what scale* $\beta$
        stops the cascade.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        Barotropic vorticity dynamics on a doubly-periodic $\beta$-plane
        (symbols as in [NOTATION](../notation)):

        $$\frac{\partial \zeta}{\partial t} + J(\psi,\zeta) + \beta\,v
          = -\nu\,(-\nabla^2)^{n_\nu}\,\zeta,
          \qquad \nabla^2\psi = \zeta,
          \qquad \mathbf{u} = (-\psi_y,\ \psi_x),$$

        so $v = \psi_x$ and $J(\psi,\zeta) = \psi_x\zeta_y - \psi_y\zeta_x$.
        $n_\nu = 1$ is Newtonian viscosity; $n_\nu > 1$ is *hyperviscosity*,
        which confines dissipation to the smallest scales.

        (Why study 2D turbulence in a GFD book at all? Because rotation and
        stratification make large-scale flow *behave* two-dimensionally:
        Ch. 8's QG dynamics is layerwise-2D by construction, and everything
        in this chapter transfers to it nearly verbatim. The atmosphere and
        ocean are, at large scales, the best 2D-turbulence laboratories in
        existence.)

        ### Two invariants, one inescapable conclusion

        With $\nu = 0$ and $\beta = 0$ the flow conserves both

        $$E = \tfrac{1}{2}\langle|\nabla\psi|^2\rangle
        \qquad\text{and}\qquad
        Z = \tfrac{1}{2}\langle\zeta^2\rangle .$$

        In spectral form $E=\int E(k)\,dk$ while $Z=\int k^2E(k)\,dk$ — the
        *same* spectrum, but enstrophy weights it by $k^2$. Now run
        **Fjørtoft's argument**: turbulence spreads the spectrum out
        (nonlinearity mixes scales — that is what turbulence *is*). But you
        cannot spread $E(k)$ while conserving both its area and its
        $k^2$-weighted area unless the *bulk of the energy* moves toward
        **small** $k$ while the *bulk of the enstrophy* moves toward
        **large** $k$. Try it with three wavenumbers and a pencil: move
        energy from $k$ to $2k$ and $k/2$; conservation forces most of it
        downscale in $k$. The conclusion is kinematic — no mechanism, no
        model, just the two conservation laws.

        This **dual cascade** is exactly backwards from 3D turbulence
        (where vortex stretching, forbidden in 2D, destroys enstrophy
        conservation and energy famously falls *down* the scales to
        dissipation). In 2D, friction at small scales is nearly harmless to
        the energy — the flow instead builds ever-larger structures.
        Kraichnan (1967) supplies the spectra: $E(k)\propto k^{-5/3}$ in
        the inverse-energy range and $E(k)\propto k^{-3}$ in the
        forward-enstrophy range.

        ### The arrest: where turbulence meets Rossby waves

        On a $\beta$-plane the upscale march does not continue forever.
        Compare timescales: an eddy of size $1/k$ turns over in
        $\tau_{turb}\sim1/(Uk)$, while a Rossby wave at that scale
        oscillates in $\tau_{wave}\sim k^2/(\beta k_x)$. At large $k$
        turbulence is faster and waves are irrelevant; at small $k$ the
        wave restoring wins and inhibits the nonlinear transfer. They match
        at

        > the **Rhines wavenumber** $\;k_\beta\simeq\sqrt{\beta/2U}$ — *the
        > scale where the inverse cascade runs into Rossby-wave stiffness
        > and stalls.* The arrest is anisotropic: purely zonal modes
        > ($k_x=0$) have **no** Rossby restoring at all, so energy funnels
        > into them — and the flow reorganizes into **zonal jets** of width
        > $\sim\pi/k_\beta$.
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
    from gfdlib import spectral, timestep, diagnostics, plotting
    return diagnostics, np, plotting, plt, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Pseudo-spectral on a doubly-periodic $[0,2\pi)^2$ grid: derivatives are
        exact in Fourier space; the Jacobian is evaluated in physical space and
        dealiased with the Orszag 2/3 rule. Both linear terms are *diagonal* in
        spectral space, so they are integrated **exactly** with an integrating
        factor,

        $$\hat L_{\mathbf k} = -\nu |\mathbf k|^{2n_\nu}
          + i\,\beta\,\frac{k_x}{|\mathbf k|^2},$$

        (the imaginary part is the Rossby-wave propagator,
        $\omega = -\beta k_x/|\mathbf k|^2$), while the advection term gets
        classical RK4 (`gfdlib.timestep.ifrk4_step`). The initial condition is
        the McWilliams (1984) random field with a spectrum peaked at $k_0$,
        rescaled to energy $E_0 = 1/2$.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    scenario = mo.ui.dropdown(
        options={
            "Decaying turbulence (McWilliams 1984)": "decay",
            "Hyperviscosity (n_nu=2) — sharper cascade": "hyper",
            "Beta-plane — zonal jets (Rhines arrest)": "beta",
        },
        value="Decaying turbulence (McWilliams 1984)",
        label="Preset scenario",
    )
    scenario
    return (scenario,)


@app.cell(hide_code=True)
def _(mo, scenario):
    _presets = {
        "decay": dict(n=128, dt=0.01,   nsteps=2000, k0=6,  nu=0.0,   nnu=1, beta=0.0),
        "hyper": dict(n=128, dt=0.01,   nsteps=2000, k0=6,  nu=1e-8,  nnu=2, beta=0.0),
        "beta":  dict(n=128, dt=0.0075, nsteps=3000, k0=10, nu=1e-7,  nnu=2, beta=20.0),
    }
    _p = _presets[scenario.value]

    n_ui = mo.ui.dropdown(
        options={"64": 64, "128": 128, "256": 256}, value=str(_p["n"]),
        label="resolution $n$",
    )
    dt_ui = mo.ui.slider(0.0025, 0.02, step=0.0025, value=_p["dt"],
                         label="time step $\\Delta t$", show_value=True)
    nsteps_ui = mo.ui.slider(500, 5000, step=500, value=_p["nsteps"],
                             label="steps", show_value=True)
    k0_ui = mo.ui.slider(2, 20, step=1, value=_p["k0"],
                         label="peak wavenumber $k_0$", show_value=True)
    nu_ui = mo.ui.slider(0.0, 1e-6, step=1e-8, value=_p["nu"],
                         label="viscosity $\\nu$", show_value=True)
    nnu_ui = mo.ui.slider(1, 4, step=1, value=_p["nnu"],
                          label="viscosity order $n_\\nu$", show_value=True)
    beta_ui = mo.ui.slider(0.0, 50.0, step=2.5, value=_p["beta"],
                           label="$\\beta$", show_value=True)

    mo.vstack([
        mo.hstack([n_ui, dt_ui, nsteps_ui], justify="start"),
        mo.hstack([k0_ui, nu_ui, nnu_ui, beta_ui], justify="start"),
        mo.md(
            r"> **In the browser (WASM):** Pyodide is ~3–10× slower than native. "
            r"Start with $n = 64$ and ≤ 2000 steps, then scale up."
        ),
    ])
    return beta_ui, dt_ui, k0_ui, n_ui, nnu_ui, nsteps_ui, nu_ui


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
    beta_ui, diagnostics, dt_ui, k0_ui, mo, n_ui, nnu_ui, np, nsteps_ui,
    nu_ui, run_btn, set_sim, spectral, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run simulation** above to start."), kind="info"
    ))

    _n, _dt, _nsteps = n_ui.value, dt_ui.value, nsteps_ui.value
    _k0, _nu, _nnu, _beta = k0_ui.value, nu_ui.value, nnu_ui.value, beta_ui.value
    _E0 = 0.5
    _nsub = max(1, _nsteps // 100)          # ~100 saved frames

    grid = spectral.Grid(_n)                 # build ONCE, reuse every step

    # McWilliams (1984) initial condition: |zeta_hat| ~ k / sqrt(1 + (k/k0)^4)
    _rng = np.random.default_rng(1234)
    _K = np.where(grid.kmag == 0, 1e-10, grid.kmag)
    _amp = _K * np.sqrt(1.0 / (1.0 + (_K / _k0) ** 4))
    zeta_hat = grid.dealias * _amp * np.exp(2j * np.pi * _rng.random(_K.shape))
    _ke, _ = diagnostics.energy_enstrophy(grid.invert_laplacian(zeta_hat), grid)
    zeta_hat = zeta_hat * np.sqrt(_E0 / _ke)

    # diagonal linear operator: (hyper)viscosity + exact Rossby propagator
    _L_op = -_nu * grid.k2 ** _nnu + 1j * _beta * grid.kx * grid.k2_inv

    def _rhs_nl(t, zh):                      # -J(psi, zeta), dealiased
        return -grid.jacobian(grid.invert_laplacian(zh), zh)

    def _energy_spectrum(zh):
        _kc, _S = diagnostics.isotropic_spectrum(
            grid.kmag * grid.invert_laplacian(zh), grid)   # |k psi_hat|^2 shells
        return _kc, 0.5 * _S                                # sums to E

    _snaps, _ubar, _t_snap, _E_hist, _Z_hist = [], [], [], [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _psi_hat = grid.invert_laplacian(zeta_hat)
            _e, _z = diagnostics.energy_enstrophy(_psi_hat, grid)
            _E_hist.append(_e); _Z_hist.append(_z); _t_snap.append(_t)
            _snaps.append(grid.ifft(zeta_hat).astype(np.float32))
            _ubar.append(grid.ifft(-grid.ddy(_psi_hat)).mean(axis=0))  # zonal-mean u(y)
        if _s < _nsteps:
            zeta_hat = timestep.ifrk4_step(zeta_hat, _rhs_nl, _dt, _L_op) * grid.dealias
            _t += _dt

    _kc0, _Ek0 = _energy_spectrum(grid.fft(_snaps[0].astype(float)))
    _kcf, _Ekf = _energy_spectrum(zeta_hat)

    set_sim({
        "snaps": _snaps, "t": np.array(_t_snap),
        "E": np.array(_E_hist), "Z": np.array(_Z_hist),
        "ubar": np.array(_ubar),                       # (frames, n)
        "k": _kcf, "Ek0": _Ek0, "Ekf": _Ekf,
        "grid": grid, "beta": _beta, "k0": _k0,
    })

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_snaps)} frames."
    ), kind="success")
    return grid, zeta_hat


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the simulation has run."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(plotting, plt, sim):
    # --- diagnostic 1: initial/final vorticity + invariants ---------------
    fig1, axs1 = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    plotting.field(axs1[0], sim["snaps"][0], grid=sim["grid"], signed=True,
                   title=f"$\\zeta$ at $t={sim['t'][0]:.1f}$", cbar=False)
    plotting.field(axs1[1], sim["snaps"][-1], grid=sim["grid"], signed=True,
                   title=f"$\\zeta$ at $t={sim['t'][-1]:.1f}$")
    axs1[2].plot(sim["t"], sim["E"] / sim["E"][0], lw=2, label="$E/E_0$")
    axs1[2].plot(sim["t"], sim["Z"] / sim["Z"][0], lw=2, color="crimson", label="$Z/Z_0$")
    axs1[2].set_xlabel("$t$"); axs1[2].set_ylim(0, 1.15)
    axs1[2].set_title("energy & enstrophy"); axs1[2].legend(); axs1[2].grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo, sim):
    frame_ui = mo.ui.slider(0, len(sim["snaps"]) - 1, step=1, value=len(sim["snaps"]) - 1,
                            label="frame", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, mo, np, plotting, plt, sim):
    # --- diagnostic 2: snapshot browser + instantaneous zonal-mean u(y) ---
    _i = frame_ui.value
    _vmax = max(float(np.abs(sim["snaps"][0]).max()), 1e-12)
    fig2, axs2 = plt.subplots(1, 2, figsize=(9.5, 4), constrained_layout=True,
                              width_ratios=[1.4, 1])
    plotting.field(axs2[0], sim["snaps"][_i], grid=sim["grid"], signed=True,
                   title=f"$\\zeta$ at $t={sim['t'][_i]:.2f}$",
                   vmin=-_vmax, vmax=_vmax)
    _y = np.arange(sim["grid"].n) * sim["grid"].dx
    axs2[1].plot(sim["ubar"][_i], _y, lw=2)
    axs2[1].axvline(0, color="k", lw=0.5)
    axs2[1].set_xlabel("$\\bar{u}(y)$"); axs2[1].set_ylabel("$y$")
    axs2[1].set_title("zonal-mean zonal flow"); axs2[1].grid(alpha=0.3)
    mo.vstack([frame_ui, fig2])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Spectrum and jet formation

        Left: the isotropic energy spectrum $E(k)$ against the Kraichnan slopes
        — $k^{-3}$ for the forward enstrophy cascade ($k > k_0$) and $k^{-5/3}$
        for the inverse energy range. Right: a Hovmöller diagram of
        $\bar u(y,t)$ — with $\beta > 0$, watch zonal stripes emerge and
        persist near the Rhines scale $k_\beta \simeq \sqrt{\beta/2U}$
        (dashed line on the spectrum).
        """
    )
    return


@app.cell(hide_code=True)
def _(np, plotting, plt, sim):
    # --- diagnostic 3: E(k) with reference slopes + ubar Hovmöller ---------
    fig3, axs3 = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)

    _k, _Ek0, _Ekf = sim["k"], sim["Ek0"], sim["Ekf"]
    _pos = _Ekf > 0
    axs3[0].loglog(_k, _Ek0, lw=1.5, alpha=0.5, label="initial")
    axs3[0].loglog(_k, _Ekf, lw=2, color="crimson", label="final")
    _kr = np.array([max(sim["k0"], 3), sim["grid"].n / 3])
    _E_at = np.interp(_kr[0], _k[_pos], _Ekf[_pos])
    axs3[0].loglog(_kr, _E_at * (_kr / _kr[0]) ** (-3.0), "k--", lw=1, label="$k^{-3}$")
    _kl = np.array([1.5, max(sim["k0"], 3)])
    axs3[0].loglog(_kl, _E_at * (_kl / _kr[0]) ** (-5.0 / 3.0), "k:", lw=1.2,
                   label="$k^{-5/3}$")
    if sim["beta"] > 0:
        _U = np.sqrt(2 * sim["E"][-1])
        axs3[0].axvline(np.sqrt(sim["beta"] / (2 * _U)), color="olive", ls="--",
                        lw=1.2, label="$k_\\beta$")
    axs3[0].set_xlabel("$k$"); axs3[0].set_ylabel("$E(k)$")
    axs3[0].set_title("isotropic energy spectrum")
    axs3[0].legend(fontsize=9); axs3[0].grid(which="both", alpha=0.2)

    _im = axs3[1].imshow(sim["ubar"].T, origin="lower", aspect="auto",
                         cmap=plotting.DIV,
                         vmin=-np.abs(sim["ubar"]).max(), vmax=np.abs(sim["ubar"]).max(),
                         extent=[sim["t"][0], sim["t"][-1], 0, sim["grid"].L])
    axs3[1].set_xlabel("$t$"); axs3[1].set_ylabel("$y$")
    axs3[1].set_title("$\\bar{u}(y,t)$ — Hovmöller")
    fig3.colorbar(_im, ax=axs3[1], shrink=0.85)
    fig3
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Selective decay.** Run the decaying preset with $\nu = 0$: $E$ stays
          flat while $Z$ falls. Why can the dealiasing filter remove enstrophy
          but almost no energy? *(Answer with Fjørtoft in hand: the cascade
          delivers enstrophy to the cutoff wavenumber, where the filter eats
          it, but the energy has gone the other way — by the time anything
          reaches the small scales, it carries lots of $\zeta^2$ and almost
          no $|\nabla\psi|^2$.)*
        - **Find the arrest.** In the β-plane preset, sweep $\beta$ from 0 to 50.
          At what $\beta$ do stripes first appear in the Hovmöller panel? Compare
          the jet spacing with $\pi/k_\beta$ (compute $U$ from the final
          energy, $U=\sqrt{2E}$). Also check the *sign* structure: are the
          eastward jets sharper than the westward ones? (They should be —
          recall Ch. 7's staircase: sharp PV risers make sharp *eastward*
          jets, while the westward flow spreads over the mixed treads.)
        - **Cascade sharpness.** Compare $n_\nu = 1$ and $n_\nu = 2$ at the same
          final $Z/Z_0$: which gives a longer $k^{-3}$ range, and why?
          (Hyperviscosity $\propto k^{2n_\nu}$ is negligible until very
          near the cutoff, so it leaves more of the inertial range
          untouched — that's its entire job.)
        - **Initial scale.** Move $k_0$ from 6 to 14: does the final vortex size
          care where the energy started? The inverse cascade erases its
          origins — a hint of why large-scale flows can have *universal*
          statistics despite wildly different forcing.

        ### What you should have seen

        Small vortices merge into a few large, long-lived coherent vortices —
        the inverse energy cascade; $E$ is nearly conserved while $Z$ decays
        (selective decay). The late-time spectrum steepens toward $k^{-3}$ above
        $k_0$. With $\beta$ on, the vorticity field elongates zonally, the
        Hovmöller panel develops persistent stripes (jets) whose spacing is set
        by the Rhines scale — the 2D cascade arrested by Rossby waves. This is
        the mechanism behind banded winds on giant planets, the multiple jets of
        the Southern Ocean, and the eddy-driven midlatitude jet.

        **Where this goes next.** Ch. 16's baroclinic instability is the
        *energy source* this chapter left unspecified — it injects eddy
        energy near the deformation radius, and the inverse cascade carries
        it upscale from there. Ch. 19 asks what this stirring does to
        anything carried by the flow (tracers, heat, PV — and PV mixing is
        what builds Ch. 7's staircase); Ch. 20 puts the turbulence in a
        basin with boundaries and gets the ocean gyres.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
