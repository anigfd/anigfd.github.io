import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 7 — Vorticity and Potential Vorticity

        **Physical question.** Place a patch of spinning fluid anywhere you
        like — can you predict the flow it creates, without solving anything
        else? Drop a second patch nearby — can you predict what happens next?
        After this chapter you should be able to state the **invertibility
        principle** (potential vorticity, plus a balance condition, plus
        boundary conditions, fixes the entire flow) and use it to explain why
        jets sit exactly where PV is mixed into a **staircase**.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### Curl the momentum equation

        Take $\partial_x(\text{$v$-equation})-\partial_y(\text{$u$-equation})$
        of the 2D momentum equations on a $\beta$-plane. The pressure
        gradient — a pure gradient — has zero curl and vanishes; what
        survives, for a nondivergent flow, is an equation for the relative
        vorticity $\zeta=v_x-u_y$ alone:

        $$\frac{D\zeta}{Dt}+\beta v=0
          \qquad\Longleftrightarrow\qquad
          \frac{D}{Dt}\underbrace{(\zeta+\beta y)}_{q}=0.$$

        The $\beta v$ term just says a parcel moving north picks up
        planetary vorticity and must shed relative vorticity to compensate —
        so the *sum* $q=\zeta+\beta y$, the **potential vorticity**, is
        carried by each parcel like a dye (symbols as in
        [NOTATION](../notation)):

        $$\frac{Dq}{Dt}=q_t+J(\psi,q)=0,\qquad \zeta=\nabla^2\psi=q-\beta y.$$

        This is the simplest member of a family. The full Ertel PV,
        $q=(\boldsymbol\omega_a\cdot\nabla\theta)/\rho$, is materially
        conserved in *any* adiabatic, frictionless stratified fluid;
        Ch. 6's shallow-water $q=\zeta-\hat\eta$ and Ch. 8's QG PV are its
        thin-layer and balanced limits. One conservation law, worn three
        ways.

        ### Invertibility: vorticity acts at a distance

        Conservation alone is bookkeeping. The power move is the
        **invertibility principle**: given $q(x,y)$ everywhere, a balance
        condition (here, nondivergence: $\zeta=\nabla^2\psi$), and boundary
        conditions (periodic), the streamfunction $\psi$ — and therefore the
        entire velocity field $\mathbf u=(-\psi_y,\psi_x)$ — is **uniquely
        determined** by one Poisson solve.

        Note what kind of operator that is: $\nabla^{-2}$ is *nonlocal*. In
        an unbounded domain $\psi(\mathbf x)=\frac{1}{2\pi}\int
        \ln|\mathbf x-\mathbf x'|\,\zeta(\mathbf x')\,d^2x'$ — a vortex
        patch induces flow **everywhere**, decaying slowly with distance,
        exactly like a charge distribution induces an electrostatic
        potential. That is why two vortices that never touch can advect each
        other, and why "where the PV is" determines "what the whole fluid is
        doing".

        The complete dynamical loop, which this notebook runs over and over:
        **invert** $q\to\psi$ (balance), **advect** $q$ with the resulting
        flow (conservation), repeat. Invert once for a static picture, or
        invert every substep to watch the flow evolve — nothing else is
        needed.
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
    from gfdlib import spectral, timestep, plotting, pv
    return np, plotting, plt, pv, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Exactly the pseudo-spectral machinery from Ch. 18
        (`gfdlib.spectral.Grid`, `gfdlib.timestep.ifrk4_step`): a Poisson
        inversion (`invert_laplacian`) for the balance step, and an
        integrating-factor RK4 for the evolution (hyperviscosity and the
        $\beta$-term treated exactly, advection via RK4). `gfdlib.pv` builds
        the initial $q$ field for each scenario below.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    scenario = mo.ui.dropdown(
        options={
            "Like-signed vortex pair — merger": "merger",
            "Opposite-signed pair — self-propelling dipole": "dipole",
            "PV staircase — jets at the risers": "staircase",
        },
        value="Like-signed vortex pair — merger",
        label="Preset scenario",
    )
    scenario
    return (scenario,)


@app.cell(hide_code=True)
def _(mo, scenario):
    _presets = {
        "merger":     dict(n=128, sep=1.2, amp=1.0, sigma=0.4, n_steps=3, beta=0.0, T=8.0),
        "dipole":     dict(n=128, sep=1.0, amp=1.0, sigma=0.35, n_steps=3, beta=0.0, T=10.0),
        "staircase":  dict(n=128, sep=1.0, amp=1.0, sigma=0.4, n_steps=3, beta=0.0, T=15.0),
    }
    _p = _presets[scenario.value]

    n_ui = mo.ui.dropdown(options={"64": 64, "128": 128, "256": 256}, value=str(_p["n"]),
                          label="resolution $n$")
    sep_ui = mo.ui.slider(0.4, 2.5, step=0.1, value=_p["sep"],
                          label="separation (vortex presets)", show_value=True)
    amp_ui = mo.ui.slider(0.2, 2.0, step=0.1, value=_p["amp"],
                          label="amplitude", show_value=True)
    sigma_ui = mo.ui.slider(0.15, 0.8, step=0.05, value=_p["sigma"],
                            label="blob width (vortex presets)", show_value=True)
    nsteps_ui = mo.ui.slider(1, 6, step=1, value=_p["n_steps"],
                             label="bands (staircase preset)", show_value=True)
    beta_ui = mo.ui.slider(0.0, 20.0, step=1.0, value=_p["beta"],
                           label="$\\beta$", show_value=True)
    T_ui = mo.ui.slider(2.0, 30.0, step=1.0, value=_p["T"],
                        label="run time $T$", show_value=True)

    mo.vstack([
        mo.hstack([n_ui, amp_ui, beta_ui], justify="start"),
        mo.hstack([sep_ui, sigma_ui, nsteps_ui], justify="start"),
        mo.hstack([T_ui], justify="start"),
    ])
    return T_ui, amp_ui, beta_ui, n_ui, nsteps_ui, sep_ui, sigma_ui


@app.cell(hide_code=True)
def _(amp_ui, n_ui, nsteps_ui, np, pv, scenario, sep_ui, sigma_ui, spectral):
    # --- build the initial q field (reactive: no Run button needed) --------
    _n, _amp, _sep, _sigma, _nsteps = (
        n_ui.value, amp_ui.value, sep_ui.value, sigma_ui.value, nsteps_ui.value
    )
    grid0 = spectral.Grid(_n)
    _xc = _yc = 0.5 * grid0.L

    if scenario.value == "merger":
        q0 = (pv.gaussian_blob(grid0, _xc - _sep / 2, _yc, _amp, _sigma)
              + pv.gaussian_blob(grid0, _xc + _sep / 2, _yc, _amp, _sigma))
    elif scenario.value == "dipole":
        q0 = (pv.gaussian_blob(grid0, _xc, _yc - _sep / 2, _amp, _sigma)
              + pv.gaussian_blob(grid0, _xc, _yc + _sep / 2, -_amp, _sigma))
    else:
        q0 = pv.staircase_pv(grid0, n_steps=_nsteps, amp=_amp)

    return grid0, q0


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Invert, don't integrate

        Left: the PV field as placed. Middle: the streamfunction and
        velocity it **instantaneously** implies — computed with a single
        Poisson solve, no time-stepping. Right: the zonal-mean $q(y)$ and
        $u(y)$ — watch how, for the staircase preset, the wind extrema line
        up with the steepest part of the PV profile, not its flat plateaus.

        **Why jets live at the risers:** for a zonal-mean profile the
        inversion reads $\bar u=-\bar\psi_y$ with
        $\bar\psi_{yy}=\bar q-\beta y$, so $\bar u_y\sim-(\bar q-\beta y)$ —
        the *wind's curvature* tracks the PV anomaly, and the wind itself
        peaks where the PV *gradient* is concentrated. Mix PV flat in a band
        (a "tread") and you've killed the gradient there; all the gradient —
        and hence a sharp eastward jet — piles up at the "riser" between
        treads. Sharp PV gradients also *resist* mixing (they support strong
        Rossby restoring, Ch. 9), so the staircase is self-reinforcing.
        """
    )
    return


@app.cell(hide_code=True)
def _(grid0, np, plotting, plt, q0):
    fig1, axs1 = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)

    plotting.field(axs1[0], q0, grid=grid0, signed=True, title="PV $q(x,y)$", cbar=False)

    _psi_hat = grid0.invert_laplacian(grid0.fft(q0))
    _psi = grid0.ifft(_psi_hat)
    _u = grid0.ifft(-grid0.ddy(_psi_hat))
    _v = grid0.ifft(grid0.ddx(_psi_hat))
    plotting.field(axs1[1], _psi, grid=grid0, signed=True, title="$\\psi(x,y)$ (inverted)", cbar=False)
    _skip = max(1, grid0.n // 16)
    axs1[1].quiver(grid0.x[::_skip, ::_skip], grid0.y[::_skip, ::_skip],
                   _u[::_skip, ::_skip], _v[::_skip, ::_skip],
                   scale=None, width=0.004, alpha=0.7, color="k")

    _q_prof = q0.mean(axis=0)
    _u_prof = _u.mean(axis=0)
    _y = grid0.y[0, :]
    _ax2b = axs1[2].twiny()
    axs1[2].plot(_q_prof, _y, color="crimson", lw=2, label="$\\bar q(y)$")
    _ax2b.plot(_u_prof, _y, color="#2563eb", lw=2, label="$\\bar u(y)$")
    axs1[2].set_xlabel("$\\bar q$", color="crimson")
    _ax2b.set_xlabel("$\\bar u$", color="#2563eb")
    axs1[2].set_ylabel("$y$")
    axs1[2].set_title("zonal means")
    axs1[2].grid(alpha=0.3)
    fig1
    return


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
    T_ui, beta_ui, grid0, mo, np, q0, run_btn, set_sim, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run evolution** above to watch this PV field evolve."),
        kind="info",
    ))

    _grid, _beta, _T = grid0, beta_ui.value, T_ui.value
    _nu, _nnu = 2e-4, 2   # fixed mild hyperviscosity, not the point of this chapter

    _zeta_hat = _grid.fft(q0) * _grid.dealias
    _L_op = -_nu * _grid.k2 ** _nnu + 1j * _beta * _grid.kx * _grid.k2_inv

    def _rhs_nl(t, zh):
        return -_grid.jacobian(_grid.invert_laplacian(zh), zh)

    _dt = 0.01
    _nsteps = max(20, int(_T / _dt))
    _nsub = max(1, _nsteps // 60)

    _snaps, _t_hist = [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _snaps.append(_grid.ifft(_zeta_hat).astype(np.float32))
            _t_hist.append(_t)
        if _s < _nsteps:
            _zeta_hat = timestep.ifrk4_step(_zeta_hat, _rhs_nl, _dt, _L_op) * _grid.dealias
            _t += _dt

    set_sim({"snaps": _snaps, "t": np.array(_t_hist), "grid": _grid})

    mo.callout(mo.md(
        f"**Done.** {_nsteps} steps · {len(_snaps)} frames."
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
    # --- diagnostic: snapshot browser of the evolving PV field -------------
    _i = frame_ui.value
    fig2, ax2 = plt.subplots(figsize=(5, 4.5), constrained_layout=True)
    plotting.field(ax2, sim["snaps"][_i], grid=sim["grid"], signed=True,
                   title=f"$q$ at $t={sim['t'][_i]:.1f}$")
    mo.vstack([frame_ui, fig2])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Merger vs. no merger.** In the merger preset, increase
          separation until the two vortices stop merging within the run
          time and instead just co-rotate. Roughly how many blob-widths
          apart is that threshold? *Context:* for idealized equal vortex
          patches the critical separation is $\approx3.2$ radii — inside it,
          each vortex's strain field tears filaments off the other faster
          than rotation can protect it; outside, they orbit like point
          vortices essentially forever. See whether Gaussian blobs land near
          the same number.
        - **Dipole speed.** In the dipole preset, decrease separation:
          does the pair speed up or slow down? *Point-vortex estimate:* each
          vortex rides the other's induced velocity,
          $U\sim\Gamma/(2\pi d)$ for circulation $\Gamma$ and separation
          $d$ — so halving $d$ should roughly double the travel distance in
          the same $T$. This self-propulsion is why dipoles are the fluid's
          preferred way to *transport* vorticity (and heat, and tracers)
          across a domain.
        - **Break the staircase.** In the staircase preset, increase
          $\beta$ from 0: do the risers (and their jets) stay put, drift,
          or get smeared out by Rossby-wave radiation? Note which way any
          drift goes, and whether the jets *survive* — sharp PV gradients
          support fast Rossby waves, which make the risers elastic rather
          than fragile.
        - **Count the jets.** Change the number of bands and confirm the
          zonal-mean panel always shows $2\times$ that many alternating
          jets — one per riser. (Periodic in $y$: $n$ treads means $2n$
          sign-alternating risers.)

        ### What you should have seen

        Every picture in this notebook comes from the same two facts:
        $q$ determines $\psi$ uniquely (invertibility), and $q$ is carried
        unchanged by the flow it determines (conservation). Two like-signed
        patches of vorticity orbit each other and merge; two opposite-signed
        patches lock together and self-propel as a dipole; and a PV field
        mixed into alternating well-mixed bands inverts into a field of
        jets sitting exactly at the sharp risers between bands — not in the
        quiet, well-mixed plateaus. This last result is the PV-staircase
        picture believed to underlie multiple-jet systems from the
        atmosphere's storm track to Jupiter's banded winds.

        **Where this goes next.** The invert-then-advect loop is the entire
        computational content of quasi-geostrophy: Ch. 8 changes only the
        inversion operator (adding a stretching term), Ch. 9 linearizes the
        same equation to get Rossby waves, and Ch. 18 runs it at high
        amplitude to get geostrophic turbulence — which *produces* the
        staircase you placed here by hand.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
