import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 1 — Kinematics & the Material Derivative

        **Physical question.** Drop a patch of dye into a fluid near a vortex
        embedded in a straining current. Some of the dye survives as a
        coherent blob circling the vortex core; the rest is stretched into
        thin filaments and swept away. After this chapter you should be able
        to predict, from the local balance of **strain** and **rotation**
        alone, which fate a given patch of fluid meets — and to state exactly
        why "the flow at this instant" (a streamline) is not the same
        question as "where does this parcel go" (a trajectory).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### Trajectories vs. streamlines

        A fluid parcel's position $\mathbf x(t)$ obeys the simplest possible
        equation in all of fluid dynamics:

        $$\frac{d\mathbf x}{dt}=\mathbf u(\mathbf x,t).$$

        A **streamline** is a snapshot: the curve everywhere tangent to
        $\mathbf u$ at one fixed instant $t$. A **trajectory** (pathline) is
        the actual solution of the equation above, integrated through time.
        For a **steady** flow ($\partial\mathbf u/\partial t=0$) these
        coincide exactly — a parcel's path traces out a streamline. The
        moment the flow is unsteady, they generically diverge: a parcel can
        cross what was a streamline a moment ago, because the streamline
        itself has since moved. This notebook uses a flow that is steady in
        the Eulerian sense (fixed velocity field) but whose *material
        derivative* still does interesting things, because what a parcel
        experiences depends on where it travels — not on time alone.

        ### The material derivative

        For any field $\phi(\mathbf x,t)$ carried by the flow, the rate of
        change **following a parcel** is

        $$\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+(\mathbf u\cdot\nabla)\phi.$$

        The second term — advection — is why a *steady* velocity field can
        still produce unsteady-looking, complicated parcel histories: a
        parcel sampling different $\nabla\phi$ as it moves feels $\phi$
        changing even though $\partial\phi/\partial t\equiv0$ at any fixed
        point. Every equation in this book (vorticity, PV, buoyancy) is a
        statement about $D(\text{something})/Dt$; this chapter is about
        building the intuition for what that operator actually does to a
        parcel's neighborhood.

        ### Strain, rotation, and the Okubo-Weiss parameter

        Expand the velocity field linearly about any point: the local
        velocity-gradient tensor $\partial u_i/\partial x_j$ splits uniquely
        into a symmetric (strain) part and an antisymmetric (rotation) part.
        In 2D, three numbers capture it completely (symbols as in
        [NOTATION](../notation)):

        $$S_n=u_x-v_y\ \text{(normal strain)},\qquad
          S_s=v_x+u_y\ \text{(shear strain)},\qquad
          \zeta=v_x-u_y\ \text{(vorticity)}.$$

        $S_n$ stretches material lines along one axis and compresses them
        along the perpendicular axis; $S_s$ does the same but along axes
        tilted $45°$; $\zeta$ rotates the neighborhood rigidly, stretching
        nothing.

        ### Where $W$ actually comes from: an eigenvalue calculation

        The claim "strain wins means exponential separation" is not a
        metaphor — it is a two-line eigenvalue computation, worth doing
        once in full. Let $\delta\mathbf x$ be the separation between two
        nearby parcels. Both obey $d\mathbf x/dt=\mathbf u(\mathbf x)$;
        subtracting and Taylor-expanding $\mathbf u$ to first order in the
        separation gives a *linear* ODE with the velocity-gradient tensor
        as its (locally frozen) matrix:

        $$\frac{d(\delta\mathbf x)}{dt}=\mathsf A\,\delta\mathbf x,\qquad
          \mathsf A=\begin{pmatrix}u_x & u_y\\ v_x & v_y\end{pmatrix}.$$

        For incompressible 2D flow $u_x+v_y=0$, so $\operatorname{tr}
        \mathsf A=0$ and the eigenvalues satisfy simply
        $\lambda^2=-\det\mathsf A$. Writing the four entries in terms of
        $(S_n,S_s,\zeta)$ — incompressibility gives $u_x=S_n/2$,
        $v_y=-S_n/2$, and inverting the definitions gives
        $v_x=(S_s+\zeta)/2$, $u_y=(S_s-\zeta)/2$ — the determinant works
        out to

        $$\det\mathsf A=u_xv_y-u_yv_x
          =-\frac{S_n^2}{4}-\frac{S_s^2-\zeta^2}{4}
          =-\frac{W}{4},
          \qquad\text{so}\qquad
          \lambda=\pm\frac{\sqrt W}{2}.$$

        The **Okubo-Weiss parameter** $W=S_n^2+S_s^2-\zeta^2$ is therefore
        nothing but (four times) the discriminant of the local separation
        dynamics:

        - $W>0$: $\lambda$ is **real**, $\pm\sqrt W/2$. Nearby parcels
          separate exponentially at rate $\sqrt W/2$ along the unstable
          eigenvector (and compress along the other) — this is where
          material lines and dye patches are stretched into filaments.
        - $W<0$: $\lambda$ is **purely imaginary**, $\pm i\sqrt{|W|}/2$.
          The separation vector rotates at frequency $\sqrt{|W|}/2$
          without growing — nearby parcels orbit each other, a coherent
          vortex core.
        - $W=0$ is the (generically thin) boundary between the two
          regimes, where the local dynamics is degenerate (pure shear:
          linear, not exponential, growth).

        One honest caveat: the derivation froze $\mathsf A$ in time, which
        is only justified if the parcel doesn't move somewhere with a
        different $\mathsf A$ before the exponential behavior expresses
        itself — i.e., if $\mathsf A$ varies slowly along trajectories.
        Near sharp features this fails, and the Okubo-Weiss criterion is
        known to over-predict filamentation there; more sophisticated
        Lagrangian diagnostics (finite-time Lyapunov exponents) fix this
        at much greater cost. For this notebook's smooth flow, $W$ is an
        excellent guide.

        This notebook's flow is a steady superposition of uniform strain
        (rate $\alpha$) and a smooth (Lamb-Oseen-like) vortex of circulation
        $\Gamma$ and core radius $\sigma$ — a classic minimal model for
        exactly this competition, used throughout the vortex-dynamics
        literature (`gfdlib.kinematics`).
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
    from gfdlib import kinematics, timestep, plotting
    return kinematics, np, plotting, plt, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        `gfdlib.kinematics.velocity_field` evaluates the (steady) strain +
        vortex velocity in closed form; `gfdlib.kinematics.okubo_weiss`
        differentiates it with central differences (deliberately not a
        hand-derived analytic Okubo-Weiss formula for the vortex term — the
        same safety pattern `gfdlib.rossby` uses for its ray equations).
        Particle trajectories are integrated with `gfdlib.timestep.rk4`,
        stepping every particle's $(x,y)$ simultaneously as one flat array —
        no new integrator needed, this is exactly $d\mathbf x/dt=\mathbf u$.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    alpha_ui = mo.ui.slider(0.0, 0.5, step=0.02, value=0.15,
                            label="strain rate $\\alpha$", show_value=True)
    Gamma_ui = mo.ui.slider(0.5, 6.0, step=0.25, value=3.0,
                            label="vortex circulation $\\Gamma$", show_value=True)
    sigma_ui = mo.ui.slider(0.2, 1.5, step=0.1, value=0.6,
                            label="vortex core radius $\\sigma$", show_value=True)
    T_ui = mo.ui.slider(2.0, 20.0, step=1.0, value=10.0,
                        label="advection time $T$", show_value=True)
    mo.vstack([
        mo.hstack([alpha_ui, Gamma_ui], justify="start"),
        mo.hstack([sigma_ui, T_ui], justify="start"),
        mo.md(r"> A ring of tracer particles is released around the vortex "
              r"core at $t=0$ and advected by the flow above. The "
              r"Okubo-Weiss field is EXACTLY steady (the flow doesn't "
              r"change), so it's shown once, underneath the moving "
              r"particles."),
    ])
    return Gamma_ui, T_ui, alpha_ui, sigma_ui


@app.cell(hide_code=True)
def _(Gamma_ui, alpha_ui, kinematics, np, sigma_ui):
    # --- Okubo-Weiss field on a fixed background grid (reactive) -----------
    _L = 4.0
    xg = np.linspace(-_L, _L, 220)
    yg = np.linspace(-_L, _L, 220)
    XG, YG = np.meshgrid(xg, yg, indexing="ij")
    W_field = kinematics.okubo_weiss(XG, YG, alpha_ui.value, Gamma_ui.value, sigma_ui.value)
    return W_field, XG, YG


@app.cell(hide_code=True)
def _(mo):
    get_traj, set_traj = mo.state(None)   # results survive the run button resetting
    return get_traj, set_traj


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Release tracer ring & advect")
    run_btn
    return (run_btn,)


@app.cell
def _(
    Gamma_ui, T_ui, alpha_ui, kinematics, mo, np, run_btn, set_traj, sigma_ui,
    timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Release tracer ring & advect** above to start."), kind="info"
    ))

    _alpha, _Gamma, _sigma, _T = alpha_ui.value, Gamma_ui.value, sigma_ui.value, T_ui.value
    _n_particles = 240
    _r0 = 1.3 * _sigma
    _theta0 = np.linspace(0, 2 * np.pi, _n_particles, endpoint=False)
    _x0 = _r0 * np.cos(_theta0)
    _y0 = _r0 * np.sin(_theta0)
    _state = np.stack([_x0, _y0])

    def _rhs(t, s):
        _x, _y = s
        _u, _v = kinematics.velocity_field(_x, _y, _alpha, _Gamma, _sigma)
        return np.stack([_u, _v])

    _dt = 0.02
    _nsteps = max(20, int(_T / _dt))
    _nsub = max(1, _nsteps // 120)

    _snaps, _t_hist = [], []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _snaps.append(_state.copy())
            _t_hist.append(_t)
        if _s < _nsteps:
            _state = timestep.rk4(_rhs, _state, _dt)
            _t += _dt

    set_traj({"snaps": _snaps, "t": np.array(_t_hist), "r0": _r0})

    mo.callout(mo.md(
        f"**Done.** {_n_particles} particles, {_nsteps} steps → {len(_snaps)} frames."
    ), kind="success")
    return


@app.cell
def _(get_traj, mo):
    traj = get_traj()
    mo.stop(traj is None, mo.callout(
        mo.md("Results will appear here once the ring has been advected."), kind="warn"
    ))
    return (traj,)


@app.cell(hide_code=True)
def _(mo, traj):
    frame_ui = mo.ui.slider(0, len(traj["snaps"]) - 1, step=1, value=len(traj["snaps"]) - 1,
                            label="frame", show_value=True)
    return (frame_ui,)


@app.cell(hide_code=True)
def _(W_field, XG, YG, frame_ui, mo, plt, traj):
    # --- diagnostic: Okubo-Weiss field + advected tracer ring --------------
    _i = frame_ui.value
    _x, _y = traj["snaps"][_i]

    fig1, ax1 = plt.subplots(figsize=(6, 5.5), constrained_layout=True)
    _wm = max(float(abs(W_field).max()), 1e-8)
    _im = ax1.pcolormesh(XG, YG, W_field, cmap="RdBu_r", vmin=-_wm, vmax=_wm, shading="auto")
    ax1.contour(XG, YG, W_field, levels=[0.0], colors="k", linewidths=1.2)
    fig1.colorbar(_im, ax=ax1, shrink=0.85, label="Okubo-Weiss $W$")
    ax1.plot(_x, _y, "o", ms=3, color="k", alpha=0.85)
    ax1.set_xlabel("$x$"); ax1.set_ylabel("$y$")
    ax1.set_aspect("equal")
    ax1.set_title(f"$t={traj['t'][_i]:.1f}$  ·  black curve = $W$=0")
    mo.vstack([frame_ui, fig1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read it:** the black contour is $W=0$ — inside it (blue,
        $W<0$) vorticity dominates and the vortex core rigidly rotates any
        particle it encloses; outside it (red, $W>0$) strain dominates and
        particles separate exponentially along the strain axis. Particles
        that started inside the $W<0$ region should stay roughly together,
        circling the core; particles that started (or wander) into $W>0$
        territory should stretch apart into a thin filament aligned with the
        local strain axis.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Check the rate, not just the sign.** Far from the vortex the
          flow is nearly pure strain, so the eigenvalue calculation
          predicts nearby particles separate like $e^{\sqrt W t/2}
          =e^{\alpha t}$ (using $W=4\alpha^2$ for pure strain). Pick two
          adjacent particles that start in the red region, estimate how
          many time units it takes their separation to grow by $e\approx
          2.7$, and compare with $1/\alpha$ from your slider.
        - **Shrink the vortex.** Decrease $\Gamma$ (or increase $\alpha$)
          until the $W<0$ region disappears entirely: no radius is safe any
          more, and the whole ring should stretch into filaments regardless
          of where it started. This is literally the mechanism behind
          **vortex stripping** — a vortex embedded in strong enough ambient
          strain cannot survive, however coherent it started.
        - **Grow the vortex.** Increase $\Gamma$ at fixed $\alpha$: the
          $W<0$ region should grow, and more of the ring should survive as
          a coherent, slowly-deforming blob circling the core rather than
          filamenting.
        - **Watch the transition in real time.** Set the ring radius near
          the $W=0$ contour (adjust $\sigma$ so $1.3\sigma$ sits close to
          the boundary) and step through frames: do particles that start
          exactly on the boundary behave more like the inside group or the
          outside group as $t$ grows? (Neither, cleanly — the boundary
          itself is not a material curve, so particles starting there
          generically end up on one side or the other.)
        - **Streamline vs. trajectory.** In your head (or on paper), sketch
          the *instantaneous* streamlines of this flow near the vortex edge
          at $t=0$ (they are simply closed curves circling the core, exactly
          like the $W$ field's symmetry suggests). Now compare to the
          particle *trajectories* you just watched. Since the flow is
          steady, these should coincide — confirm that a filamenting
          particle's path really does trace out one connected streamline,
          just one that happens to pass through the strain-dominated region.

        ### What you should have seen

        The same steady flow does two completely different things to
        nearby parcels depending only on local geometry: a rigid, coherent
        rotation inside the Okubo-Weiss boundary, and an exponential,
        filament-producing stretch outside it. Nothing about the *governing
        equation* changed between the two regimes — $d\mathbf x/dt=\mathbf
        u(\mathbf x)$ throughout — only the *local structure* of $\mathbf u$
        did. This strain-vs-rotation competition, quantified by $W$, is the
        single most-used diagnostic for finding coherent vortices (ocean
        eddies, atmospheric storms) in real velocity data, and it is the
        seed of every later chapter's vorticity dynamics: Ch. 7's
        invertibility principle is precisely the statement that $\zeta$
        (one of the three numbers making up $W$) determines the whole flow.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
