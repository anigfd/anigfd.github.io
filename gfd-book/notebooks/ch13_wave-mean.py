import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 13 — Wave-Mean-Flow Interaction (Intro)

        **Physical question.** Two ways a wave, even a perfectly linear one,
        leaves a permanent mark on the mean state it rides on. First: does
        a floating object end up exactly where it started once a wave train
        passes, or does it drift? (It drifts — **Stokes drift** — even
        though the *Eulerian*-mean velocity at any fixed point is exactly
        zero.) Second: what happens to a wave as it approaches a latitude
        where the mean flow matches its own phase speed? (It cannot get
        there — a **critical layer** — and in doing so it sets up exactly
        the conditions for depositing its momentum into the mean flow.)
        After this chapter you should be able to derive the Stokes-drift
        formula for a simple wave from scratch, and explain why a Rossby
        wave's wavenumber and amplitude both diverge as it approaches a
        critical layer.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part A — Stokes drift

        ### Governing equations

        Take the simplest possible wave: a single Fourier mode,
        $u'(x,t)=A\cos(kx-\omega t)$, with **zero Eulerian mean** at every
        fixed point ($\langle u'(x,t)\rangle_t=0$ trivially, by definition
        of a pure oscillation). Now ask a different question: not "what is
        the velocity at a fixed point," but "where does a marked fluid
        parcel actually end up." To leading order in the wave amplitude,
        a parcel released at $x_0$ has displacement
        $\xi(t)=\int_0^tu'(x_0,t')\,dt'=\dfrac{A}{\omega}\sin(\omega t-kx_0)$
        (symbols as in [NOTATION](../notation)). The parcel's TRUE velocity
        is the field evaluated at its ACTUAL position $x_0+\xi(t)$, not at
        $x_0$; Taylor-expanding to the next order in $\xi$,

        $$u'(x_0+\xi,t)\approx u'(x_0,t)+\xi\,\frac{\partial u'}{\partial x}\bigg|_{x_0,t}.$$

        Time-averaging: the first term is the Eulerian mean, exactly zero.
        The second term is new — it is a genuine, second-order-in-amplitude
        mean velocity that a parcel actually experiences, even though no
        term in the original equation looks like a "mean flow":

        $$u_S=\left\langle \xi\,\frac{\partial u'}{\partial x}\right\rangle
          =\frac{A^2k}{2\omega}.$$

        This is the **Stokes drift** — positive (in the direction of phase
        propagation) for $k,\omega>0$. The three-line derivation above is
        the entire content of "waves move fluid even when their Eulerian
        mean says they don't": it is a pure consequence of evaluating the
        velocity at the parcel's *displaced* position rather than its mean
        one, nothing more exotic.
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
    from gfdlib import wavemean, timestep
    return np, plt, timestep, wavemean


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Numerical scheme (Part A)

        `gfdlib.wavemean.wave_velocity` is the prescribed oscillatory
        field; a single particle is advected through it with
        `gfdlib.timestep.rk4` — exactly Ch. 1's $d\mathbf x/dt=\mathbf u$
        machinery, now applied to a genuinely time-dependent flow. The
        measured long-time drift rate is compared directly against
        `gfdlib.wavemean.stokes_drift`'s closed form. Fully reactive — cheap
        enough for no Run button.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    A_ui = mo.ui.slider(0.01, 0.15, step=0.01, value=0.06,
                        label="$A$ (wave velocity amplitude)", show_value=True)
    k_ui = mo.ui.slider(0.5, 3.0, step=0.25, value=1.0, label="$k$", show_value=True)
    omega_ui = mo.ui.slider(0.5, 4.0, step=0.25, value=2.0, label="$\\omega$", show_value=True)
    nperiods_ui = mo.ui.slider(10, 80, step=5, value=40, label="run time (periods)",
                               show_value=True)
    mo.vstack([mo.hstack([A_ui, k_ui], justify="start"),
               mo.hstack([omega_ui, nperiods_ui], justify="start")])
    return A_ui, k_ui, nperiods_ui, omega_ui


@app.cell(hide_code=True)
def _(A_ui, k_ui, nperiods_ui, np, omega_ui, timestep, wavemean):
    # --- advect one particle through the wave (reactive, cheap) ------------
    A, k, omega, n_periods = A_ui.value, k_ui.value, omega_ui.value, nperiods_ui.value
    T_period = 2 * np.pi / omega
    substeps = 40
    dt = T_period / substeps
    nsteps = n_periods * substeps

    def _rhs(t, s):
        return np.array([wavemean.wave_velocity(s[0], t, A, k, omega)])

    _state = np.array([0.0])
    xs, ts = [0.0], [0.0]
    _t = 0.0
    for _s in range(nsteps):
        _state = timestep.rk4(_rhs, _state, dt, _t)
        _t += dt
        xs.append(float(_state[0])); ts.append(_t)
    xs, ts = np.array(xs), np.array(ts)

    strobe_idx = np.arange(0, nsteps + 1, substeps)
    x_strobe, t_strobe = xs[strobe_idx], ts[strobe_idx]
    predicted = wavemean.stokes_drift(A, k, omega)
    measured = xs[-1] / ts[-1]
    return T_period, measured, predicted, t_strobe, ts, x_strobe, xs


@app.cell(hide_code=True)
def _(measured, mo, np, plt, predicted, t_strobe, ts, x_strobe, xs):
    fig1, axs1 = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    axs1[0].plot(ts, xs, lw=0.7, alpha=0.6, color="#2563eb", label="$x(t)$, every substep")
    axs1[0].plot(ts, predicted * ts, "k--", lw=1.5, label="predicted drift line")
    axs1[0].set_xlabel("$t$"); axs1[0].set_ylabel("$x$")
    axs1[0].set_title("raw trajectory"); axs1[0].legend(fontsize=8); axs1[0].grid(alpha=0.3)

    axs1[1].plot(t_strobe, x_strobe, "o-", ms=4, color="crimson",
                label="stroboscopic (once per period)")
    axs1[1].plot(ts, predicted * ts, "k--", lw=1.5, label="predicted drift line")
    axs1[1].set_xlabel("$t$"); axs1[1].set_ylabel("$x$")
    axs1[1].set_title("stroboscopic view (oscillation removed)")
    axs1[1].legend(fontsize=8); axs1[1].grid(alpha=0.3)
    mo.vstack([
        fig1,
        mo.md(f"**Measured drift rate:** {measured:.5f}  ·  **Predicted "
              f"($A^2k/2\\omega$):** {predicted:.5f}  ·  relative difference "
              f"{abs(measured-predicted)/predicted*100:.1f}%"),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read it:** the left panel is what you'd actually record
        with a tracked buoy — a fast oscillation riding on a barely-visible
        slow drift. The right panel samples the trajectory once every wave
        period, exactly removing the fast oscillation by construction: the
        remaining points should fall almost exactly on the predicted
        straight line. This "strobe" trick — sample once per cycle — is
        precisely how Stokes drift is separated from orbital motion in real
        wave-buoy data.

        ## Part B — The Rossby-wave critical layer

        ### Governing equations

        A Rossby wave riding on a background zonal shear $U(y)$ (not just
        planetary curvature, as in Ch. 9) has a **Doppler-shifted**
        dispersion relation:

        $$\omega(y,k,l)=U(y)\,k-\frac{\beta k}{k^2+l^2}.$$

        The medium is steady ($U,\beta$ don't depend on $t$), so $\omega$
        is exactly conserved along any ray (Ch. 9's WKB machinery, now with
        $U(y)$ added to the effective medium). Hamilton's ray equations for
        $(y,l)$ at fixed, conserved $k$:

        $$\frac{dy}{dt}=\frac{\partial\omega}{\partial l}=\frac{2\beta
          kl}{(k^2+l^2)^2},\qquad
          \frac{dl}{dt}=-\frac{\partial\omega}{\partial y}=-k\frac{dU}{dy}.$$

        Now the key argument. Suppose the ray approaches a latitude $y_c$
        where $U(y_c)=\omega/k$ — the **Doppler-shifted phase speed matches
        the mean flow**. As $y\to y_c$, the term $U(y)k$ alone would already
        equal $\omega$, so the remaining term $-\beta k/(k^2+l^2)$ must
        shrink to zero to keep $\omega$ fixed at its conserved value. That
        forces $k^2+l^2\to\infty$: the **meridional wavenumber diverges**.
        And since $dy/dt\propto l/(k^2+l^2)^2\to0$ as $l\to\infty$ faster
        than $l$ grows, the ray's own approach speed vanishes — it gets
        arbitrarily close to $y_c$ but, in this linear, inviscid theory,
        never arrives in finite time. This is the **critical layer**: a
        genuine mathematical singularity of linear ray theory, not a
        numerical artifact (verified directly below: $\omega$ stays exactly
        conserved along the ray even as $y$ and $l$ both change
        enormously).

        ### What diverges, exactly? The wave-action bookkeeping

        The wavenumber blowing up is only half the singularity — the wave's
        *amplitude* diverges too, and the argument is short enough to do.
        In a slowly-varying medium the conserved wave quantity is not
        energy (the shear can exchange energy with the wave) but **wave
        action** $A=E/\hat\omega$, the energy density divided by the
        *intrinsic* frequency $\hat\omega=\omega-U(y)k=-\beta k/(k^2+l^2)$
        — a WKB result previewed here and derived properly in Ch. 23. For
        a steady wave train, the meridional **flux of action** must be
        independent of $y$:

        $$c_{g,y}\,\frac{E}{\hat\omega}
          =\frac{2\beta kl/K^4}{-\beta k/K^2}\,E
          =-\frac{2l}{K^2}\,E=\text{const}
          \qquad\Longrightarrow\qquad
          E\propto\frac{K^2}{|l|}\sim|l|\ \to\ \infty.$$

        The energy *density* piles up without bound as the ray slows down —
        the same reason ocean swell steepens on a beach (slower group
        speed, same flux, so energy accumulates). And the divergence has a
        characteristic anatomy: since $E\sim\tfrac12K^2|\hat\psi|^2$, the
        streamfunction amplitude actually *shrinks*
        ($|\hat\psi|^2\sim E/K^2\propto1/|l|$), while the
        along-flow velocity $|u'|\sim|l||\hat\psi|$ and the vorticity
        $|\zeta'|\sim K^2|\hat\psi|$ grow
        like $|l|^{1/2}$ and $|l|^{3/2}$: the wave is being sheared into
        finer and finer, more and more intense filaments. No linear theory
        survives that indefinitely — either viscosity erases the filaments
        or the vorticity overturns nonlinearly, and in both cases the
        wave's momentum flux is deposited *at* the critical layer, which
        is the entire wave-mean interaction this chapter is named for.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Numerical scheme (Part B)

        `gfdlib.wavemean.ray_rhs_shear` evaluates the ray equations above
        directly (an exact symbolic derivative of the dispersion relation,
        unlike Ch. 9's central-difference approach — safe here since the
        formula is simple enough to differentiate by hand without risk).
        Integrated with `gfdlib.timestep.rk4`. The background flow is
        linear shear, $U(y)=\Lambda y$, chosen only because it gives a
        clean, exactly-locatable critical layer
        $y_c=\omega_0/(k\Lambda)$ to compare against — the qualitative
        behavior (wavenumber blow-up, vanishing approach speed) holds for
        any $U(y)$ with $dU/dy\neq0$ near $y_c$.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    Lambda_ui = mo.ui.slider(0.3, 2.0, step=0.1, value=1.0, label="$\\Lambda=dU/dy$",
                             show_value=True)
    k2_ui = mo.ui.slider(0.5, 2.0, step=0.1, value=1.0, label="$k$", show_value=True)
    beta2_ui = mo.ui.slider(0.1, 1.5, step=0.1, value=0.5, label="$\\beta$", show_value=True)
    l0_ui = mo.ui.slider(0.2, 3.0, step=0.1, value=1.0, label="launch $l_0$", show_value=True)
    T2_ui = mo.ui.slider(1.0, 15.0, step=1.0, value=8.0, label="run time $T$", show_value=True)
    mo.vstack([
        mo.hstack([Lambda_ui, k2_ui], justify="start"),
        mo.hstack([beta2_ui, l0_ui], justify="start"),
        mo.hstack([T2_ui], justify="start"),
        mo.md(r"> Launched at $y_0=0$. $y_c=\omega_0/(k\Lambda)$ is computed "
              r"from your own launch conditions and marked on every plot."),
    ])
    return Lambda_ui, T2_ui, beta2_ui, k2_ui, l0_ui


@app.cell(hide_code=True)
def _(Lambda_ui, T2_ui, beta2_ui, k2_ui, l0_ui, np, timestep, wavemean):
    # --- integrate the ray (reactive, cheap: a handful of ODEs) ------------
    Lam, k2, beta2, l0, T2 = (
        Lambda_ui.value, k2_ui.value, beta2_ui.value, l0_ui.value, T2_ui.value
    )

    def _dUdy(y):
        return Lam

    def _U(y):
        return Lam * y

    omega0 = wavemean.rossby_shear_dispersion(0.0, k2, l0, _U, beta2)
    y_c = omega0 / (k2 * Lam)

    _state = np.array([0.0, l0])
    dt2 = 1e-3
    nsteps2 = int(T2 / dt2)
    nsub2 = max(1, nsteps2 // 400)

    y_hist, l_hist, cgy_hist, t_hist2, omega_hist = [], [], [], [], []
    _t = 0.0
    for _s in range(nsteps2 + 1):
        if _s % nsub2 == 0:
            _y, _l = _state
            _cgy, _ = wavemean.ray_rhs_shear(_state, k2, beta2, _dUdy)
            y_hist.append(_y); l_hist.append(_l); cgy_hist.append(_cgy)
            t_hist2.append(_t)
            omega_hist.append(wavemean.rossby_shear_dispersion(_y, k2, _l, _U, beta2))
        if _s < nsteps2:
            _state = timestep.rk4(lambda t, s: wavemean.ray_rhs_shear(s, k2, beta2, _dUdy),
                                  _state, dt2)
            _t += dt2

    y_hist, l_hist, cgy_hist = np.array(y_hist), np.array(l_hist), np.array(cgy_hist)
    t_hist2, omega_hist = np.array(t_hist2), np.array(omega_hist)
    return cgy_hist, l_hist, omega0, omega_hist, t_hist2, y_c, y_hist


@app.cell(hide_code=True)
def _(cgy_hist, l_hist, mo, omega0, omega_hist, plt, t_hist2, y_c, y_hist):
    fig2, axs2 = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)

    axs2[0].plot(t_hist2, y_hist, lw=2, color="#2563eb")
    axs2[0].axhline(y_c, color="crimson", ls="--", lw=1.2, label=f"$y_c={y_c:.3f}$")
    axs2[0].set_xlabel("$t$"); axs2[0].set_ylabel("$y$")
    axs2[0].set_title("meridional position"); axs2[0].legend(fontsize=8); axs2[0].grid(alpha=0.3)

    axs2[1].plot(t_hist2, l_hist, lw=2, color="#d97706")
    axs2[1].set_xlabel("$t$"); axs2[1].set_ylabel("$l$")
    axs2[1].set_title("meridional wavenumber"); axs2[1].grid(alpha=0.3)

    axs2[2].plot(t_hist2, cgy_hist, lw=2, color="#16a34a")
    axs2[2].axhline(0, color="gray", lw=0.6)
    axs2[2].set_xlabel("$t$"); axs2[2].set_ylabel("$dy/dt$")
    axs2[2].set_title("meridional group velocity"); axs2[2].grid(alpha=0.3)

    _domega = float(np.max(np.abs(omega_hist - omega0)))
    mo.vstack([
        fig2,
        mo.md(f"**Max drift in $\\omega$ over the whole run:** {_domega:.2e} "
              f"(should be at the level of numerical roundoff/RK4 truncation "
              f"— $\\omega$ is exactly conserved by the ray equations, not "
              f"approximately)."),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read it:** with the default parameters the ray first
        drifts *away* from $y_c$ and reaches a **turning point** where
        $l=0$ and $dy/dt=0$ momentarily — exactly Ch. 9's "turning
        latitude" caustic, arising here for the same reason (a sign change
        in the meridional group velocity). Only after turning around does
        the ray head toward $y_c$, with $|l|$ growing steadily and $dy/dt$
        shrinking toward (but never quite reaching) zero as it approaches.
        Two different singular points, back to back, in one trajectory —
        do not confuse them: the turning point is a reflection (the ray
        continues, just in the other direction); the critical-layer
        approach is a genuine asymptote (the ray never gets there, in this
        theory).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Confirm exact $\omega$ conservation.** Read the number printed
          above the ray-tracing figure: it should be many orders of
          magnitude smaller than $\omega_0$ itself, confirming the ray
          equations conserve $\omega$ to numerical precision, not
          approximately.
        - **Locate the critical layer yourself.** Before running, compute
          $y_c=\omega_0/(k\Lambda)$ by hand from your slider values (you'll
          need $\omega_0$ from $U(0)k-\beta k/(k^2+l_0^2)$ first) and check
          it against the dashed line the notebook draws automatically.
        - **Push the run time.** Increase $T$ toward its maximum: does the
          ray ever cross $y_c$, or does it keep asymptoting closer without
          quite arriving? At what point does the wavenumber panel's growth
          start looking numerically suspicious (a hint that the underlying
          WKB ray theory itself — which assumes slowly-varying wavenumber —
          is starting to break down as $l$ grows very large, a limitation
          of the THEORY, not just the numerics)?
        - **Vary $\beta$.** Set $\beta$ near its minimum: does the ray reach
          closer to $y_c$ in the same run time, or not as close? Relate
          your answer to $dy/dt$'s explicit $\beta$ dependence.

        ### What you should have seen

        Two genuinely different mechanisms by which a linear wave leaves a
        permanent trace on the mean state: Stokes drift, a direct
        consequence of evaluating a wave field at a parcel's *true*
        (displaced) position rather than its mean one; and the critical
        layer, where a wave riding on shear approaches — but cannot reach —
        the latitude of matching phase speed, its wavenumber and (by
        wave-action arguments) its amplitude both diverging as it does.
        Neither mechanism required any dissipation or nonlinearity to
        derive — both are exact features of *linear* theory. What linear
        theory cannot tell you is what happens next: real critical layers,
        once any dissipation or finite amplitude is added, absorb the
        wave's momentum and genuinely accelerate or decelerate the mean
        flow there — the actual "wave-mean interaction" this chapter is
        named for, and the mechanism behind the stratosphere's
        quasi-biennial oscillation. Part VII's Ch. 23 returns to this with
        the proper wave-activity/pseudomomentum machinery needed to go
        beyond the linear, conservative theory built here.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
