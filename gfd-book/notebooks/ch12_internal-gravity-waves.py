import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 12 — Internal Gravity Waves

        **Physical question.** Shake a stratified fluid at a single point and
        the energy doesn't spread out in circles — it concentrates into four
        narrow **beams** at an angle fixed by the forcing frequency, not by
        the size of the source (the "St. Andrew's Cross"). After this chapter
        you should be able to predict that beam angle from $\omega$, $N$, and
        $f$ alone, and explain why a beam **reflects** wherever it reaches a
        depth where the local stratification $N(z)$ drops below $\omega$.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The spring: buoyancy

        Displace a parcel upward by $\delta z$ in a stratified fluid and it
        finds itself denser than its new surroundings; the restoring
        acceleration is $\ddot{\delta z}=-N^2\,\delta z$ with

        > **buoyancy frequency** $\;N^2=-\dfrac{g}{\rho_0}\dfrac{d\rho}
        > {dz}=\dfrac{db}{dz}$ — *the frequency at which a vertically
        > displaced parcel oscillates about its rest level.* Ocean
        > thermocline: $N\sim10^{-2}\,$s$^{-1}$ (10-minute period);
        > atmosphere: $N\sim10^{-2}\,$s$^{-1}$ too. Compare
        > $f\sim10^{-4}\,$s$^{-1}$: buoyancy is a stiff spring, rotation a
        > soft one, and the two-decade gap between them is where internal
        > waves live.

        A parcel displaced *along a slope* at angle $\theta$ from vertical
        feels only the component of gravity along its path, and oscillates
        at $N\cos\theta$ — slower for steeper-from-vertical paths. Hold that
        thought: it *is* the dispersion relation.

        ### The wave equation

        Linear, incompressible motion in a vertical $(x,z)$ plane, uniform
        rotation $f$, buoyancy frequency $N(z)$ possibly varying with depth.
        Using the streamfunction $u=\psi_z,\ w=-\psi_x$ (symbols as in
        [NOTATION](../notation) — note this sign convention is specific to
        the vertical plane), take the curl of the momentum equations to get
        an equation for $q=\nabla^2\psi$, then eliminate the buoyancy (via
        $b_t=-N^2w$) and the rotation-coupled along-front velocity (via
        $v_t=-fu$) by taking one more time derivative:

        $$\frac{\partial^2 q}{\partial t^2} = -\Big(N^2(z)\,\psi_{xx} +
          f^2\,\psi_{zz}\Big) + S(x,z,t).$$

        Read the right side as two springs: buoyancy $N^2$ acts on
        *horizontal* wiggles ($\psi_{xx}$: tilted columns), rotation $f^2$
        on *vertical* wiggles ($\psi_{zz}$: sheared layers).

        ### Dispersion: frequency depends on angle, not size

        For **constant** $N$, substitute a plane wave
        $\sim e^{i(kx+mz-\omega t)}$ (so $\nabla^2\to-(k^2+m^2)$):

        $$\omega^2=\frac{N^2k^2+f^2m^2}{k^2+m^2}=N^2\cos^2\theta+f^2\sin^2\theta,
          \qquad \cos\theta=\sqrt{\frac{\omega^2-f^2}{N^2-f^2}},$$

        where $\theta$ is the angle of $\mathbf k=(k,m)$ from *vertical* —
        i.e. the angle of the wave **crests** (and the parcel motion, which
        is along crests for an incompressible transverse wave) from
        *horizontal*. This is the parcel-on-a-slope frequency from above,
        now with rotation stiffening the near-horizontal paths.

        Two consequences make internal waves genuinely strange:

        - $\omega$ depends **only on the angle** of $\mathbf k$, never its
          magnitude. Force at one frequency and you select a *direction*,
          not a wavelength — energy from a point source radiates along four
          **beams** at the fixed angle $\theta$, instead of spreading in
          rings the way every intuition from surface waves says it should.
        - Because $\omega$ is constant along any ray through $\mathbf k$-space
          origin, $\mathbf c_g=\nabla_{\mathbf k}\omega$ must be
          **perpendicular to $\mathbf k$**: group velocity is *along the
          crests*, phase velocity *across* them. In particular,
          downward-marching phase means **upward**-traveling energy — so
          when a mooring record shows phase lines descending, the energy
          source is *below*, not above. Every observational oceanographer
          learns to make this sign flip; in the snapshot browser you can
          watch it happen.

        If the geometric argument feels too slick, verify it by brute
        force. Differentiating $\omega^2=(N^2k^2+f^2m^2)/(k^2+m^2)$ with
        respect to each wavenumber component (write $K^2=k^2+m^2$ and use
        the quotient rule) gives

        $$\mathbf c_g=\left(\frac{\partial\omega}{\partial k},
          \frac{\partial\omega}{\partial m}\right)
          =\frac{1}{\omega K^2}\Big(k\,(N^2-\omega^2),\ m\,(f^2-\omega^2)\Big),$$

        and the dot product with $\mathbf k=(k,m)$ collapses using the
        dispersion relation itself:

        $$\mathbf c_g\cdot\mathbf k
          =\frac{k^2N^2+m^2f^2-\omega^2K^2}{\omega K^2}
          =\frac{\omega^2K^2-\omega^2K^2}{\omega K^2}=0
          \qquad\text{— exactly, not approximately.}$$

        The component form also shows *why* the vertical components
        oppose: for $f<\omega<N$, the factor $(N^2-\omega^2)$ is positive
        but $(f^2-\omega^2)$ is **negative** — so $c_{g,z}$ and the
        vertical phase speed $\omega m/K^2$ always carry opposite signs.
        The mooring-record sign flip is not folklore; it is that minus
        sign.

        Waves exist only where $f\le|\omega|\le N(z)$. Where $N(z)$ drops
        below $\omega$, the wave cannot propagate and **reflects** at that
        turning level — no boundary condition is imposed there; the wave
        equation does it on its own.
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
    from gfdlib import spectral, plotting, internalwaves
    return internalwaves, np, plotting, plt, spectral


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        The $(x,z)$ plane is doubly periodic and solved with the same
        `gfdlib.spectral.Grid` used in the horizontal chapters — its $y$-axis
        stands in for $z$ here. Each step inverts $\nabla^2\psi=q$ spectrally,
        evaluates $\psi_{xx},\psi_{zz}$ spectrally, then multiplies by
        $N^2(z)$ and $f^2$ **in physical space** so $N$ may vary with depth
        (`gfdlib.internalwaves.step_leapfrog`). Time-stepping is classical
        Stormer-Verlet (leapfrog) on $q$: since $\omega(k,m)$ never exceeds
        $N$ regardless of wavenumber (it is a weighted average of $N^2$ and
        $f^2$), the scheme is unconditionally stable at any resolution for a
        fixed $\Delta t\lesssim 1/N_{\max}$. A Gaussian point source
        oscillating at frequency $\omega$ forces the domain center.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    scenario = mo.ui.dropdown(
        options={
            "Uniform stratification — St. Andrew's Cross": "uniform",
            "Thermocline duct — reflection at a turning level": "duct",
        },
        value="Uniform stratification — St. Andrew's Cross",
        label="Preset scenario",
    )
    scenario
    return (scenario,)


@app.cell(hide_code=True)
def _(mo, scenario):
    _presets = {
        "uniform": dict(n=128, Nmean=1.0, Namp=0.0, f=0.1, omega=0.5, nT=10.0),
        "duct":    dict(n=128, Nmean=0.65, Namp=0.35, f=0.1, omega=0.5, nT=10.0),
    }
    _p = _presets[scenario.value]

    n_ui = mo.ui.dropdown(
        options={"64": 64, "128": 128, "256": 256}, value=str(_p["n"]),
        label="resolution $n$",
    )
    Nmean_ui = mo.ui.slider(0.4, 1.2, step=0.05, value=_p["Nmean"],
                            label="$\\hat N_{mean}$", show_value=True)
    Namp_ui = mo.ui.slider(0.0, 0.5, step=0.05, value=_p["Namp"],
                           label="$\\hat N_{amp}$ (0 = uniform)", show_value=True)
    f_ui = mo.ui.slider(0.0, 0.3, step=0.02, value=_p["f"],
                        label="$\\hat f$", show_value=True)
    omega_ui = mo.ui.slider(0.15, 0.9, step=0.05, value=_p["omega"],
                            label="forcing $\\hat\\omega$", show_value=True)
    nT_ui = mo.ui.slider(4.0, 20.0, step=1.0, value=_p["nT"],
                         label="run time (wave periods)", show_value=True)

    mo.vstack([
        mo.hstack([n_ui, Nmean_ui, Namp_ui], justify="start"),
        mo.hstack([f_ui, omega_ui, nT_ui], justify="start"),
        mo.md(
            r"> $\hat N(z)=\hat N_{mean}+\hat N_{amp}\cos(2\pi(z-z_c)/L)$ — "
            r"stratification peaks at the source depth $z_c$ and dips below "
            r"it. Need $\hat f<\hat\omega<\hat N_{mean}+\hat N_{amp}$ for "
            r"waves to exist at the source."
        ),
    ])
    return Namp_ui, Nmean_ui, f_ui, n_ui, nT_ui, omega_ui


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
    Namp_ui, Nmean_ui, f_ui, internalwaves, mo, n_ui, np, nT_ui, omega_ui,
    run_btn, set_sim, spectral,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run simulation** above to start."), kind="info"
    ))

    _n, _Nmean, _Namp, _f, _omega, _nT = (
        n_ui.value, Nmean_ui.value, Namp_ui.value, f_ui.value, omega_ui.value,
        nT_ui.value,
    )
    _L = 2 * np.pi

    grid = spectral.Grid(_n, L=_L)               # build ONCE, reuse every step
    _xc = _zc = 0.5 * _L
    _Nz = _Nmean + _Namp * np.cos(2 * np.pi * (grid.y - _zc) / _L)
    _N2 = _Nz ** 2

    _sigma_s = _L / 40
    _source_env = np.exp(
        -((grid.x - _xc) ** 2 + (grid.y - _zc) ** 2) / (2 * _sigma_s ** 2)
    )

    _Twave = 2 * np.pi / _omega
    _dt = _Twave / 40
    _nsteps = int(_nT * _Twave / _dt)
    _nsub = max(1, _nsteps // 80)

    _q, _q_prev = np.zeros_like(grid.x), np.zeros_like(grid.x)
    _snaps, _t_hist = [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _psi_hat = grid.invert_laplacian(grid.fft(_q))
            _w = grid.ifft(-grid.ddx(_psi_hat))
            _snaps.append(_w.astype(np.float32))
            _t_hist.append(_t / _Twave)
        if _s < _nsteps:
            _src = _source_env * np.sin(_omega * _t)
            _q, _q_prev = internalwaves.step_leapfrog(
                _q, _q_prev, grid, _N2, _f, _src, _dt
            ), _q
            _t += _dt

    _theta = internalwaves.beam_angle(_omega, _Nz.max(), _f)

    set_sim({
        "snaps": _snaps, "t": np.array(_t_hist), "Nz": _Nz[0, :],
        "grid": grid, "f": _f, "omega": _omega, "theta": _theta,
        "xc": _xc, "zc": _zc,
    })

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps on a {_n}×{_n} grid · {len(_snaps)} frames "
        f"· predicted beam angle from horizontal $\\theta={np.degrees(_theta):.1f}°$."
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
def _(mo):
    mo.md(
        r"""
        ### Theory: predicted frequency vs. angle

        Before looking at the field: this is what $\omega(\theta)$ looks like
        for the chosen $\hat N_{mean}$ and $\hat f$, with the forcing
        $\hat\omega$ and its predicted beam angle marked.
        """
    )
    return


@app.cell(hide_code=True)
def _(Nmean_ui, f_ui, np, omega_ui, plt):
    # --- diagnostic 1: theory curve, no simulation needed -------------------
    _N, _f, _omega = Nmean_ui.value, f_ui.value, omega_ui.value
    _theta_arr = np.linspace(0, np.pi / 2, 300)
    _omega_arr = np.sqrt(_N ** 2 * np.cos(_theta_arr) ** 2 + _f ** 2 * np.sin(_theta_arr) ** 2)

    fig1, ax1 = plt.subplots(figsize=(5.5, 4), constrained_layout=True)
    ax1.plot(np.degrees(_theta_arr), _omega_arr, lw=2)
    ax1.axhline(_N, color="crimson", ls="--", lw=1, label="$\\hat N_{mean}$")
    ax1.axhline(_f, color="olive", ls="--", lw=1, label="$\\hat f$")
    if _f < _omega < _N:
        _cos_th = np.sqrt(np.clip((_omega ** 2 - _f ** 2) / (_N ** 2 - _f ** 2), 0, 1))
        _th0 = np.degrees(np.arccos(_cos_th))
        ax1.axhline(_omega, color="k", lw=1)
        ax1.axvline(_th0, color="k", lw=1)
        ax1.plot([_th0], [_omega], "ko", ms=6)
    ax1.set_xlabel(r"$\theta$ from vertical (deg)"); ax1.set_ylabel(r"$\hat\omega$")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo, sim):
    frame_ui = mo.ui.slider(0, len(sim["snaps"]) - 1, step=1, value=len(sim["snaps"]) - 1,
                            label="frame (wave periods)", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(frame_ui, mo, np, plotting, plt, sim):
    # --- diagnostic 2: N(z) profile + turning levels, and the w(x,z) field -
    _grid = sim["grid"]
    _i = frame_ui.value
    _w = sim["snaps"][_i]
    _z = _grid.y[0, :]
    _Nz = sim["Nz"]

    fig2, axs2 = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True,
                              width_ratios=[1, 1.5])
    axs2[0].plot(_Nz, _z, lw=2)
    axs2[0].axvline(sim["omega"], color="k", ls="--", lw=1.2, label="$\\hat\\omega$")
    axs2[0].axhline(sim["zc"], color="gray", lw=0.8, ls=":", label="source depth")
    _crossings = _z[:-1][np.diff(np.sign(_Nz - sim["omega"])) != 0]
    for _zt in _crossings:
        axs2[0].axhline(_zt, color="crimson", lw=0.8, ls="--")
    axs2[0].set_xlabel("$\\hat N(z)$"); axs2[0].set_ylabel("$z$")
    axs2[0].set_title("stratification + turning levels")
    axs2[0].legend(fontsize=8); axs2[0].grid(alpha=0.3)

    _vmax = max(float(np.abs(sim["snaps"][-1]).max()), 1e-12)
    plotting.field(axs2[1], _w, grid=_grid, signed=True,
                   title=f"$w$ at $t={sim['t'][_i]:.1f}$ wave periods",
                   vmin=-_vmax, vmax=_vmax)
    axs2[1].plot(sim["xc"], sim["zc"], "k*", ms=10)

    mo.vstack([frame_ui, fig2])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Read off the angle.** In the uniform preset, use a ruler (or the
          pixel grid) on the final frame to estimate the beam angle from
          horizontal, and compare it to $\theta$ printed above the plots.
          Then move $\hat\omega$ and predict *before re-running*: higher
          frequency $\to$ steeper or shallower beams? (From
          $\cos\theta\propto\sqrt{\omega^2-f^2}$: higher $\omega$ means
          crests closer to vertical, i.e. beams closer to vertical too.)
        - **Watch phase vs. group.** Step through consecutive frames along
          one beam of the X: the *crests inside the beam* march across it
          (perpendicular to the beam) while the *beam envelope* extends
          along itself. You are watching $\mathbf c_g\perp\mathbf c_p$
          directly — the single weirdest verified prediction of this
          dispersion relation.
        - **Turn off rotation.** Set $\hat f=0$: the beams should sharpen
          slightly (check the theory curve — how much does $\theta$ actually
          change between $\hat f=0$ and $\hat f=0.3$ at fixed $\hat\omega$?
          Rotation only matters to near-inertial waves, $\omega\to f$).
        - **Find the duct.** In the "Thermocline duct" preset, read the two
          turning-level depths off the left panel (where $\hat N(z)=\hat
          \omega$) and confirm the wave field on the right stays confined
          between them. This is a waveguide: the real ocean's main
          thermocline ducts internal-tide energy across whole basins this
          way.
        - **Break the duct.** Increase $\hat\omega$ toward $\hat N_{mean}+
          \hat N_{amp}$: the turning levels should move apart until they
          leave the domain, and the beams should reach the periodic boundary
          and reconnect, as in the uniform case.

        ### What you should have seen

        A point source produces four beams forming an "X" (a St. Andrew's
        Cross), at an angle set only by $\hat\omega,\hat N,\hat f$ — not by
        the source size or the grid resolution. With $\hat N_{amp}>0$, the
        beams **curve** as they move into weaker stratification, and
        **reflect** genuinely — from the wave equation itself, with no ad hoc
        rule — at the turning level where the local $\hat N(z)$ drops to
        $\hat\omega$, trapping the wave field in a duct around the source
        depth. This is the same mechanism (a frequency-dependent turning
        point) that governs internal-tide generation at ocean ridges and
        mountain-wave reflection in the stratified atmosphere.

        **Where this goes next.** Internal waves are the fast, unbalanced
        motion that Part II's QG world filtered out — the two chapters are
        complementary halves of the same fluid. When these waves reach
        amplitudes where they overturn, or meet a **critical layer** where
        the background flow speed matches their phase speed, they break and
        deposit their momentum into the mean flow (Ch. 13) — the mechanism
        that drives the stratospheric QBO and much of the deep ocean's
        mixing (Ch. 21's $\kappa$ is largely *made* of broken internal
        waves). And the stability of stratified shear flow against
        overturning is exactly Ch. 17's Richardson-number story.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
