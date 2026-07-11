import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 9 — Rossby Waves

        **Physical question.** A Rossby wave's group velocity depends on its
        wavenumber, and its wavenumber refracts as it travels through a
        changing background flow — so a wave packet doesn't travel in a
        straight line, it follows a *ray*. After this chapter you should be
        able to write the Rossby-wave dispersion relation and its group
        velocity, set up the WKB ray equations, and explain why — on a
        sphere, with the right background flow — those rays trace out exact
        **great circles**.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The restoring mechanism, then the algebra

        Rossby waves need no surface and no stratification — their spring is
        the planetary vorticity gradient itself. Displace a line of parcels
        northward: by Ch. 7's $D(\zeta+\beta y)/Dt=0$ they acquire
        $\zeta<0$ (clockwise spin), and that spin advects their neighbors —
        the western neighbor northward, the eastern one southward. The
        displacement pattern therefore slides **westward**. That one
        geometric argument fixes the sign of everything below.

        Now the algebra. Linearize the barotropic vorticity equation about
        rest: with $q'=\nabla^2\psi'$, the equation
        $q'_t+\beta\psi'_x=0$ admits plane waves
        $\psi'\sim e^{i(kx+ly-\omega t)}$ ($\nabla^2\to-K^2$,
        $K^2=k^2+l^2$), giving (symbols as in [NOTATION](../notation)):

        $$\omega=\frac{-\beta k}{k^2+l^2},\qquad
          \mathbf c_g=\Big(\frac{\partial\omega}{\partial k},\frac{\partial\omega}{\partial l}\Big)
          =\frac{\beta}{(k^2+l^2)^2}\big(k^2-l^2,\ 2kl\big).$$

        Read the asymmetry carefully: **phase** always moves westward
        ($\omega/k<0$ whenever $k\neq0$), but the **group velocity** — the
        velocity of energy, and the velocity that matters for this chapter —
        can point *east* (short zonal waves, $k^2>l^2$) or *west* (long
        waves). Crests and energy routinely travel in opposite directions;
        watching crests will actively mislead you about where the energy is
        going.

        ### Rays: waves in a slowly varying medium

        A slowly-varying background — a jet, or (as here) planetary
        curvature itself — makes $\omega$ a function of position as well as
        wavenumber, $\omega(\mathbf x,\mathbf k)$. Provided the background
        varies on scales much larger than a wavelength (**the WKB
        condition**), the packet keeps a well-defined local $\mathbf k$, and
        $(\mathbf x,\mathbf k)$ evolve as canonically conjugate variables —
        exactly position and momentum in classical mechanics, with $\omega$
        as the Hamiltonian:

        $$\frac{d\mathbf x}{dt}=\frac{\partial\omega}{\partial\mathbf k},
          \qquad\frac{d\mathbf k}{dt}=-\frac{\partial\omega}{\partial\mathbf x}.$$

        Two conservation laws come free, and both get used below: $\omega$
        is constant along a ray (time-independent medium), and any
        wavenumber component conjugate to a symmetry of the background is
        conserved (here: zonal symmetry $\Rightarrow$ fixed zonal
        wavenumber index $n$).

        ### On the sphere: why great circles

        Take a background solid-body super-rotation
        $U(\phi)=\Omega_s a\cos\phi$ on a planet rotating at $\Omega$. A
        **stationary** wave ($\omega=0$, the kind a mountain range or a
        persistent tropical heat source forces) must satisfy
        $K^2=K_s^2$ with

        $$K_s^2=\frac{2(\Omega+\Omega_s)}{\Omega_s a^2}
          \quad\text{— the same at every latitude,}$$

        because both the advecting flow $U$ and the meridional gradient of
        absolute vorticity scale by the *same* factor $(\Omega+\Omega_s)$
        for solid-body flow. A latitude-independent $|\mathbf k|$ is
        precisely the condition under which Hamilton's equations integrate
        to **exact great circles** (Hoskins & Karoly 1981) — the ray bends
        neither toward nor away from the equator, just as a geodesic on a
        homogeneous sphere. The diagnostic below verifies "exact" directly,
        with no wiggle room.
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
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)
    from gfdlib import timestep, rossby
    return np, plt, rossby, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        `gfdlib.rossby` deliberately does **not** hand-derive the ray
        equations symbolically: `ray_rhs` evaluates
        $\partial\omega/\partial\mathbf k,\ \partial\omega/\partial\mathbf x$
        by central differences of a single scalar function,
        `dispersion_omega` — correctness then rests only on that one
        formula, which the great-circle check below verifies independently.
        Integration is `gfdlib.timestep.rk4` on the 3-state ray
        $(\lambda,\phi,m)$ (the zonal wavenumber index $n$ is exactly
        conserved by the background's zonal symmetry, so it's a fixed
        parameter, not a fourth evolved variable). This is cheap — a handful
        of ODEs, no PDE — so the notebook is fully reactive: no Run button.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    Omega_s_ui = mo.ui.slider(0.05, 1.0, step=0.05, value=0.3,
                              label="$\\Omega_s$ (background rotation)", show_value=True)
    phi0_ui = mo.ui.slider(-60.0, 60.0, step=5.0, value=20.0,
                           label="launch latitude $\\phi_0$ (deg)", show_value=True)
    nrays_ui = mo.ui.slider(1, 7, step=1, value=5,
                            label="number of rays", show_value=True)
    T_ui = mo.ui.slider(5.0, 40.0, step=5.0, value=25.0,
                        label="run time $T$", show_value=True)

    mo.vstack([
        mo.hstack([Omega_s_ui, phi0_ui], justify="start"),
        mo.hstack([nrays_ui, T_ui], justify="start"),
        mo.md(r"> Rays fan out eastward-to-westward (launch angles $0$ to "
              r"$160°$ from local east) from the same starting point. "
              r"$\Omega_s>0$ is required — stationary Rossby waves need a "
              r"background **eastward** flow to exist at all."),
    ])
    return Omega_s_ui, T_ui, nrays_ui, phi0_ui


@app.cell(hide_code=True)
def _(Omega_s_ui, T_ui, np, nrays_ui, phi0_ui, rossby, timestep):
    # --- integrate each ray (reactive: cheap, no Run button needed) --------
    _Omega_s, _phi0, _T = Omega_s_ui.value, np.radians(phi0_ui.value), T_ui.value
    _angles = np.radians(np.linspace(0, 160, nrays_ui.value))
    _dt = 0.01
    _nsteps = max(50, int(_T / _dt))

    rays = []
    for _alpha0 in _angles:
        _state, _n = rossby.launch_state(0.0, _phi0, _alpha0, _Omega_s)
        _traj = [_state.copy()]
        for _ in range(_nsteps):
            _state = timestep.rk4(lambda t, s: rossby.ray_rhs(s, _n, _Omega_s), _state, _dt)
            _traj.append(_state.copy())
        _traj = np.array(_traj)
        _dev = rossby.great_circle_deviation(_traj[:, 0], _traj[:, 1])
        rays.append({"lam": _traj[:, 0], "phi": _traj[:, 1], "alpha0": _alpha0,
                     "max_dev": float(np.max(_dev))})
    return (rays,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### A fan of great circles

        Left: the rays on the globe. Right: the same rays on an
        equirectangular map — great circles look like sinusoidal S-curves
        here, exactly like real-world flight-path maps.
        """
    )
    return


@app.cell(hide_code=True)
def _(np, plt, rays):
    fig1 = plt.figure(figsize=(11, 5))
    ax3d = fig1.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig1.add_subplot(1, 2, 2)

    _u, _v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
    ax3d.plot_wireframe(np.cos(_u) * np.sin(_v), np.sin(_u) * np.sin(_v), np.cos(_v),
                        color="lightgray", linewidth=0.3, alpha=0.5)

    _colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(rays)))
    for _ray, _c in zip(rays, _colors):
        _lam, _phi = _ray["lam"], _ray["phi"]
        ax3d.plot(np.cos(_phi) * np.cos(_lam), np.cos(_phi) * np.sin(_lam), np.sin(_phi),
                  color=_c, lw=1.5)
        _lam_wrap = np.degrees(((_lam + np.pi) % (2 * np.pi)) - np.pi)
        ax2d.plot(_lam_wrap, np.degrees(_phi), ".", ms=1.5, color=_c)

    ax3d.set_box_aspect([1, 1, 1]); ax3d.set_title("globe")
    ax2d.set_xlim(-180, 180); ax2d.set_ylim(-90, 90)
    ax2d.set_xlabel("longitude"); ax2d.set_ylabel("latitude")
    ax2d.set_title("equirectangular map"); ax2d.grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("### How exact is “exact”?")
    return


@app.cell(hide_code=True)
def _(mo, np, rays):
    _lines = [
        f"- ray at launch angle {np.degrees(_r['alpha0']):5.0f}°: "
        f"max deviation from its fitted great-circle plane = **{_r['max_dev']:.2e}**"
        for _r in rays
    ]
    mo.md(
        "Each ray's position, projected onto the plane through the sphere's "
        "center fit from its own trajectory, should be zero for an exact "
        "great circle:\n\n" + "\n".join(_lines) +
        "\n\nThese are floating-point roundoff, not approximation error — "
        "the numerical-differentiation ray equations reproduce the "
        "Hoskins & Karoly theorem essentially exactly."
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Break the great circle.** Ray-trace with a background flow that
          is *not* solid-body (e.g. edit the notebook to use
          $U=U_0\cos^2\phi$ instead of $\cos\phi$) — the stationary
          wavenumber $K_s(\phi)$ will no longer be latitude-independent.
          Does the great-circle deviation diagnostic stay near zero, or
          grow? *This is the point of the exercise:* on the real Earth,
          $K_s(\phi)$ varies, rays refract toward latitudes of larger
          $K_s$, and regions where $K_s^2<0$ (easterlies!) are forbidden —
          which is why observed wave trains are *bent* great circles, and
          why stationary Rossby waves cannot propagate through the tropical
          easterlies.
        - **Turning latitudes.** Launch a ray nearly due north/south (large
          $\alpha_0$ near $90°$) and watch its latitude in the map view: it
          should reach a maximum, turn around, and head back — the ray's
          own version of a great circle's inclination-limited latitude
          range. In the refraction language: at the turning latitude the
          meridional wavenumber $m$ passes through zero, all of $K_s$ is in
          the zonal component, and the ray reflects.
        - **Vary $\Omega_s$.** At fixed launch angle, how does increasing
          $\Omega_s$ change $K_s$ (read the formula: $K_s^2\propto
          1+\Omega/\Omega_s$, so *stronger* westerlies mean *smaller*
          $K_s$ — longer stationary waves) and, visually, how "tightly
          wound" the great circle looks on the map?

        ### What you should have seen

        Every ray, launched in a different direction from the same point,
        traces a *different* great circle — but each one is exact, to
        floating-point precision, verified by a check that never assumed
        the answer. This is the same ray theory that predicts how energy
        from tropical heating anomalies (like El Niño) propagates into
        midlatitudes as quasi-stationary Rossby wave trains, arcing
        poleward-then-equatorward across the globe — real atmospheric
        teleconnection patterns are, to good approximation, segments of
        exactly this kind of great-circle ray.

        **Where this goes next.** This closes the book's linear-wave
        thread for balanced flow: Ch. 6 met these waves as the slow mode of
        adjustment, Ch. 8 gave them their finite-$L_R$ dispersion, and here
        they carry energy along rays. Part IV plays the same
        game — dispersion relation, group velocity, rays — for *internal
        gravity* waves in a stratified fluid (Ch. 12), where the geometry
        gets stranger: energy propagates *along* wave crests, not across
        them. And when Rossby waves grow too big to stay linear, they break
        and mix PV — which is where Ch. 18's turbulence and Ch. 7's
        staircases come from.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
