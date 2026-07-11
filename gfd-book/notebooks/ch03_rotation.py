import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 3 — Effects of Rotation

        **Physical question.** Two questions, both about what rotation alone
        does to a fluid. First: give a parcel a push and remove every other
        force — what does it do? (Not travel in a straight line.) Second:
        force a rapidly-rotating fluid at one level only — does the rest of
        the column notice? (It does, and it shouldn't.) After this chapter
        you should be able to derive the exact trajectory of an unforced
        parcel on an $f$-plane, and state precisely why a homogeneous,
        rapidly-rotating fluid organizes into vertically-rigid **Taylor
        columns**.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part A — Inertial oscillations

        ### Governing equations

        Remove every force from the horizontal momentum equation except
        Coriolis (symbols as in [NOTATION](../notation)):

        $$\frac{du}{dt}=fv,\qquad \frac{dv}{dt}=-fu.$$

        Write $w=u+iv$: the pair collapses to a single complex ODE,
        $dw/dt=-if w$, with exact solution $w(t)=w_0e^{-ift}$ — a vector of
        *constant length* rotating clockwise ($f>0$, Northern Hemisphere) at
        constant angular frequency $f$. Integrating once more gives the
        parcel's position: a **circle** of radius $|w_0|/f$, traversed once
        every **inertial period** $T=2\pi/f$.

        This is the simplest possible consequence of rotation, and it is
        directly observed: current-meter records from the upper ocean, once
        the tides and mean flow are filtered out, are dominated almost
        everywhere by a spectral peak at exactly the local inertial
        frequency $f$ — parcels genuinely do coast around these circles,
        excited by any transient push (a storm's wind stress, most
        commonly) and only slowly damped by friction.
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
    from gfdlib import rotation, timestep, balance, plotting
    return balance, np, plotting, plt, rotation, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Numerical scheme

        `gfdlib.rotation.inertial_trajectory` is the closed-form solution
        above — cheap enough to be fully reactive, no time-stepping needed
        for the plot itself. As an independent check, `gfdlib.rotation.
        inertial_rhs` is also stepped with `gfdlib.timestep.rk4` and
        compared against the analytic curve directly on the plot.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    f_ui = mo.ui.slider(0.3, 3.0, step=0.1, value=1.0,
                        label="$f$ (Coriolis parameter)", show_value=True)
    u0_ui = mo.ui.slider(0.0, 2.0, step=0.1, value=1.0,
                         label="$u_0$ (initial eastward speed)", show_value=True)
    v0_ui = mo.ui.slider(-2.0, 2.0, step=0.1, value=0.0,
                         label="$v_0$ (initial northward speed)", show_value=True)
    nT_ui = mo.ui.slider(0.5, 5.0, step=0.5, value=2.0,
                         label="run time (inertial periods)", show_value=True)
    mo.vstack([
        mo.hstack([f_ui, nT_ui], justify="start"),
        mo.hstack([u0_ui, v0_ui], justify="start"),
    ])
    return f_ui, nT_ui, u0_ui, v0_ui


@app.cell(hide_code=True)
def _(f_ui, np, nT_ui, rotation, timestep, u0_ui, v0_ui):
    # --- analytic trajectory + independent RK4 check (both cheap, reactive) -
    _f, _u0, _v0, _nT = f_ui.value, u0_ui.value, v0_ui.value, nT_ui.value
    T_period = 2 * np.pi / _f
    _t = np.linspace(0, _nT * T_period, 400)
    x_a, y_a, u_a, v_a = rotation.inertial_trajectory(_t, _u0, _v0, _f)

    _dt = T_period / 200
    _nsteps = int(round(_nT * T_period / _dt))
    _rk4_state = np.array([_u0, _v0])
    x_n, y_n = [0.0], [0.0]
    for _s in range(_nsteps):
        _rk4_state = timestep.rk4(lambda t, s: rotation.inertial_rhs(s, _f), _rk4_state, _dt)
        x_n.append(x_n[-1] + _rk4_state[0] * _dt)
        y_n.append(y_n[-1] + _rk4_state[1] * _dt)
    x_n, y_n = np.array(x_n), np.array(y_n)
    return T_period, u_a, v_a, x_a, x_n, y_a, y_n


@app.cell(hide_code=True)
def _(T_period, np, plt, u0_ui, u_a, v0_ui, v_a, x_a, x_n, y_a, y_n):
    fig1, axs1 = plt.subplots(1, 2, figsize=(10.5, 4.5), constrained_layout=True)
    axs1[0].plot(x_a, y_a, lw=2.5, color="crimson", label="analytic")
    axs1[0].plot(x_n, y_n, "--", lw=1.2, color="k", alpha=0.7, label="RK4 (independent check)")
    axs1[0].plot(0, 0, "k*", ms=10)
    axs1[0].set_xlabel("$x$"); axs1[0].set_ylabel("$y$")
    axs1[0].set_title("trajectory"); axs1[0].set_aspect("equal"); axs1[0].legend(fontsize=8)
    axs1[0].grid(alpha=0.3)

    axs1[1].plot(u_a, v_a, lw=2.5, color="#2563eb")
    axs1[1].plot(u0_ui.value, v0_ui.value, "ko", ms=6)
    axs1[1].axhline(0, color="gray", lw=0.5); axs1[1].axvline(0, color="gray", lw=0.5)
    axs1[1].set_xlabel("$u$"); axs1[1].set_ylabel("$v$")
    axs1[1].set_title(f"hodograph  ·  $T=2\\pi/f={T_period:.2f}$")
    axs1[1].set_aspect("equal"); axs1[1].grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this (Part A)

        - **Read the period.** Set run time to exactly 1 inertial period and
          confirm the trajectory closes back on its starting point — then
          verify $T=2\pi/f$ against the value printed on the hodograph panel.
        - **Radius scales with speed, not direction.** Compare $(u_0,v_0)=
          (1,0)$ and $(0,1)$: same radius $|w_0|/f$, just rotated — the
          circle's size depends only on speed.
        - **Speed up rotation.** Increase $f$ at fixed $u_0,v_0$: the circle
          should shrink (radius $\propto1/f$) and spin faster (period
          $\propto1/f$) simultaneously — both from the same $f$.

        ## Part B — The Taylor-Proudman limit

        ### Governing equations

        Ch. 5's thermal-wind relation, $f\,u_{g,z}=-b_y$, is not a special
        case of Taylor-Proudman — it *contains* it. Set the buoyancy
        perturbation to zero ($b\equiv$ const, a perfectly homogeneous
        fluid) and the relation degenerates immediately:

        $$b_y\equiv0\ \Longrightarrow\ f\,u_{g,z}=0\ \Longrightarrow\
          \frac{\partial u_g}{\partial z}=\frac{\partial v_g}{\partial z}=0.$$

        This is the **Taylor-Proudman theorem**: a steady, geostrophically
        balanced, homogeneous, rotating flow cannot vary along the rotation
        axis. Physically: with no buoyancy to support a pressure that
        depends on $z$, hydrostatic balance ($\phi_z=b=0$) forces $\phi$ —
        and therefore the geostrophic velocity it determines — to be
        *literally the same at every height*. A fluid column moves as a
        single rigid unit, a **Taylor column**, from bottom to top.

        Below, `gfdlib.balance`'s exact frontal buoyancy field and
        thermal-wind integral (identical code to Ch. 5) are evaluated with
        the buoyancy contrast $\Delta b$ swept from a normal frontal value
        down toward zero — watch the vertical shear disappear continuously,
        not as a separate calculation, but as the $\Delta b\to0$ limit of
        the exact same formula.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    db_ui = mo.ui.slider(0.0, 1.0, step=0.02, value=0.5,
                         label="$\\Delta b$ (buoyancy contrast — 0 = homogeneous)",
                         show_value=True)
    Ly_b_ui = mo.ui.slider(0.3, 2.5, step=0.1, value=1.0,
                           label="$L_y$ (front half-width)", show_value=True)
    Htrop_b_ui = mo.ui.slider(0.3, 2.0, step=0.1, value=1.0,
                              label="$H_{trop}$", show_value=True)
    f_b_ui = mo.ui.slider(0.2, 2.0, step=0.1, value=1.0,
                          label="$f$", show_value=True)
    mo.vstack([
        mo.hstack([db_ui, Ly_b_ui], justify="start"),
        mo.hstack([Htrop_b_ui, f_b_ui], justify="start"),
        mo.md(r"> $N^2$ is fixed at 0 in this section — any residual $z$-tilt "
              r"of the isotherms would itself carry no $y$-gradient and so "
              r"cannot restore any shear (Ch. 5's $N^2$ slider showed exactly "
              r"this). Only $\Delta b$ (the FRONTAL, $y$-varying part) matters "
              r"here."),
    ])
    return Htrop_b_ui, Ly_b_ui, db_ui, f_b_ui


@app.cell(hide_code=True)
def _(Htrop_b_ui, Ly_b_ui, balance, db_ui, f_b_ui, np):
    # --- u_g(y=0, z) profile and its z-range as a function of db -----------
    _Htrop = Htrop_b_ui.value
    z_prof = np.linspace(1e-3, 2 * _Htrop, 300)
    u_profile = balance.thermal_wind_u(0.0, z_prof, db_ui.value, Ly_b_ui.value, _Htrop, f_b_ui.value)

    db_sweep = np.linspace(0.0, 1.0, 40)
    _shear_range = []
    for _db in db_sweep:
        _u = balance.thermal_wind_u(0.0, z_prof, _db, Ly_b_ui.value, _Htrop, f_b_ui.value)
        _shear_range.append(float(_u.max() - _u.min()))
    shear_range = np.array(_shear_range)
    return z_prof, db_ui, shear_range, u_profile, db_sweep


@app.cell(hide_code=True)
def _(db_sweep, z_prof, db_ui, mo, plt, shear_range, u_profile):
    fig2, axs2 = plt.subplots(1, 2, figsize=(10, 4.3), constrained_layout=True)
    axs2[0].plot(u_profile, z_prof, lw=2.5, color="#2563eb")
    axs2[0].axvline(0, color="gray", lw=0.6)
    axs2[0].set_xlabel("$u_g(y{=}0,z)$"); axs2[0].set_ylabel("$z$")
    axs2[0].set_title(f"vertical profile · $\\Delta b={db_ui.value:.2f}$")
    axs2[0].grid(alpha=0.3)

    axs2[1].plot(db_sweep, shear_range, lw=2, color="crimson")
    axs2[1].plot(db_ui.value,
                 float(u_profile.max() - u_profile.min()), "ko", ms=8, zorder=5)
    axs2[1].set_xlabel("$\\Delta b$"); axs2[1].set_ylabel("vertical shear range")
    axs2[1].set_title("shear range vs. buoyancy contrast")
    axs2[1].grid(alpha=0.3)
    mo.vstack([
        fig2,
        mo.md("Left: the profile flattens toward a vertical line — exactly "
              "columnar flow — as $\\Delta b\\to0$. Right: the shear range "
              "is **exactly linear** in $\\Delta b$ (the black dot marks "
              "your current slider value on that line) — Taylor-Proudman "
              "isn't approached asymptotically here, it's approached "
              "linearly, all the way to a literal zero."),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this (Part B)

        - **Drive $\Delta b$ to exactly zero.** Confirm the left panel's
          profile becomes a perfectly vertical line — $u_g$ genuinely
          independent of $z$, not just "small." Taylor-Proudman is exact
          in this limit, not an approximation.
        - **Check the linearity.** The right panel's line should pass
          through the origin. Why must it? (`thermal_wind_u` is linear in
          $\Delta b$ by inspection of the closed-form solution in
          `gfdlib.balance` — the shear range inherits that linearity
          exactly, with no fitting involved.)
        - **Connect to Part A.** Taylor columns and inertial oscillations
          are usually taught as unrelated facts about rotation. They are
          not: both are statements about what a *homogeneous* rotating
          fluid does when nothing else (buoyancy, friction) is present —
          Part A says an unforced parcel moves in circles; Part B says a
          forced column moves rigidly. Ch. 4's regime map will show both
          living in the same corner of parameter space ($Bu\ll1$).

        ### What you should have seen

        Two faces of the same fact: rotation alone, with no other physics
        to lean on, produces motion with a very specific character — closed
        circular orbits for a free parcel (Part A), and rigid vertical
        coherence for a balanced flow (Part B). Both traced directly back
        to $f$ alone: no buoyancy needed for the second, no forcing at all
        for the first. Real Taylor columns are visible directly in
        laboratory experiments (a rotating tank with a bump on the bottom
        casts a "shadow" of blocked flow all the way to the free surface)
        and in the deep, weakly-stratified abyssal ocean, where currents
        measurably follow bottom topography contours far more than they
        follow the overlying wind or buoyancy forcing.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
