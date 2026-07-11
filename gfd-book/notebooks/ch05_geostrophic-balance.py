import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 5 — Geostrophic & Hydrostatic Balance

        **Physical question.** Draw a temperature field — warm air to the
        south, cold air to the north — and rotation turns it into a wind,
        without you specifying any wind directly. After this chapter you
        should be able to reconstruct the vertical shear of the geostrophic
        wind from a horizontal temperature (buoyancy) gradient alone, and
        explain why the strongest winds in the atmosphere sit right at the
        tropopause.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The scaling argument, done honestly

        Start from the horizontal momentum and hydrostatic equations for a
        rotating, Boussinesq fluid (symbols as in [NOTATION](../notation)),
        with $\phi=p/\rho_0$ the kinematic pressure and $b=-g\rho'/\rho_0$
        the buoyancy:

        $$\frac{Du}{Dt}-fv=-\phi_x,\qquad \frac{Dv}{Dt}+fu=-\phi_y,\qquad
          \phi_z=b.$$

        (The third equation — hydrostatic balance — is itself a scaling
        result: for motions much wider than they are deep, $L\gg H$, the
        vertical acceleration $Dw/Dt$ is smaller than $b$ by a factor
        $\sim(H/L)^2$ and drops out. Synoptic weather systems have
        $H/L\sim10^{-2}$, so this is excellent.)

        Now put sizes on the horizontal terms. Let $U$ be a typical velocity
        and $L$ a typical horizontal scale, and let the flow evolve on its
        own advective time $T\sim L/U$. Then, term by term:

        $$\underbrace{\frac{Du}{Dt}}_{\sim\,U^2/L}\;-\;
          \underbrace{fv}_{\sim\,fU}\;=\;-\phi_x .$$

        The ratio of acceleration to Coriolis defines the chapter's one
        dimensionless number:

        > **Rossby number** $\;Ro=\dfrac{U^2/L}{fU}=\dfrac{U}{fL}$ — *how
        > much the flow accelerates in the time rotation takes to turn it.*
        > Midlatitude weather: $U\sim10\,$m/s, $L\sim1000\,$km,
        > $f\sim10^{-4}\,$s$^{-1}$ $\Rightarrow Ro\sim0.1$. Ocean mesoscale
        > eddies: $U\sim0.1\,$m/s, $L\sim100\,$km $\Rightarrow Ro\sim0.01$.
        > A bathtub vortex: $Ro\sim10^{6}$ — rotation of the Earth is
        > irrelevant there, which is why this chapter is about planets, not
        > bathtubs.

        At $Ro\ll1$ the accelerations are negligible at leading order and
        the pressure gradient can only be balanced by the Coriolis force —
        **geostrophic balance**:

        $$u_g=-\frac{\phi_y}{f},\qquad v_g=\frac{\phi_x}{f}.$$

        Two things to internalize before moving on:

        - The wind blows **along** isobars, not down the pressure gradient —
          a fluid parcel pushed toward low pressure is deflected sideways by
          rotation until the two forces balance. This is why weather maps
          are useful: contours of pressure *are* streamlines.
        - Geostrophic flow is **automatically nondivergent**
          ($u_{g,x}+v_{g,y}=0$ for constant $f$), and comparing with the
          streamfunction convention in NOTATION shows $\psi=\phi/f$: the
          pressure field, rescaled, *is* the streamfunction. Every balanced
          model later in the book (Chs. 7, 8, 16, 18) leans on this.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Thermal wind: the whole derivation is one cross-derivative

        Geostrophy fixes the wind at each level from the pressure at that
        level, but pressure at different levels is tied together by
        hydrostatic balance. Differentiate $u_g=-\phi_y/f$ in $z$, swap the
        order of derivatives, and substitute $\phi_z=b$:

        $$f\,u_{g,z}=-\phi_{yz}=-(\phi_z)_y=-b_y,
          \qquad\text{and likewise}\qquad f\,v_{g,z}=b_x.$$

        That is the **thermal wind** relation, and the derivation really is
        that short — the physics is in reading it correctly:

        - The **vertical shear** of the wind is set entirely by the
          **horizontal** buoyancy (temperature) gradient. Cold air to the
          north ($b_y<0$) forces the westerly wind to *increase* with
          height ($u_{g,z}>0$) — the midlatitude jet stream in one line.
        - It fixes only the *shear*: adding any depth-independent
          (barotropic) flow leaves the relation untouched. Temperature
          alone cannot tell you the absolute wind — you need the wind at
          one reference level, which is why this notebook integrates
          upward from $u_g=0$ at the ground.

        This notebook reconstructs $u_g(y,z)$ by integrating
        $u_{g,z}=-b_y/f$ upward from an idealized frontal buoyancy field

        $$b(y,z)=N^2z+\Delta b\,\tanh(y/L_y)\cos\!\Big(\frac{\pi z}{2H_{trop}}\Big),$$

        a warm–cold contrast of strength $\Delta b$ concentrated in a front
        of half-width $L_y$, whose sign reverses at $z=H_{trop}$ (an
        idealized tropopause: above it, the meridional temperature gradient
        flips, as it does in the real lower stratosphere). All variables are
        nondimensional — think of $y,z$ in units of front width and
        tropopause height, and the wind in units of $\Delta b\,H_{trop}/(fL_y)$.
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
    from gfdlib import plotting, balance
    return balance, np, plotting, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Every field here is a **closed-form expression**
        (`gfdlib.balance`): the buoyancy field, its exact $y$-derivative, and
        the thermal-wind integral are all evaluated directly — no
        discretization, no time-stepping, no `Grid`. This is purely
        diagnostic calculus, so the notebook recomputes live as you move the
        sliders below; there is no "Run" button.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    db_ui = mo.ui.slider(0.1, 2.0, step=0.1, value=1.0,
                         label="$\\Delta b$ (buoyancy contrast)", show_value=True)
    Ly_ui = mo.ui.slider(0.3, 2.5, step=0.1, value=1.0,
                         label="$L_y$ (front half-width)", show_value=True)
    Htrop_ui = mo.ui.slider(0.3, 2.0, step=0.1, value=1.0,
                            label="$H_{trop}$ (tropopause height)", show_value=True)
    f_ui = mo.ui.slider(0.2, 2.0, step=0.1, value=1.0,
                        label="$f$ (Coriolis parameter)", show_value=True)
    N2_ui = mo.ui.slider(0.0, 1.0, step=0.1, value=0.3,
                         label="$N^2$ (background stratification)", show_value=True)
    mo.vstack([
        mo.hstack([db_ui, Ly_ui, Htrop_ui], justify="start"),
        mo.hstack([f_ui, N2_ui], justify="start"),
        mo.md(r"> $N^2$ tilts the isotherms for realism but carries no "
              r"meridional gradient, so it has **no effect** on the wind — "
              r"watch the $u_g$ panel to confirm."),
    ])
    return N2_ui, Htrop_ui, Ly_ui, db_ui, f_ui


@app.cell(hide_code=True)
def _(Htrop_ui, N2_ui, Ly_ui, balance, db_ui, f_ui, np):
    # --- evaluate the closed-form fields on a (y,z) grid --------------------
    db, Ly, Htrop, f, N2 = db_ui.value, Ly_ui.value, Htrop_ui.value, f_ui.value, N2_ui.value
    y = np.linspace(-4.0, 4.0, 240)
    z = np.linspace(1e-3, 2 * Htrop, 200)
    YY, ZZ = np.meshgrid(y, z, indexing="ij")

    b_field = balance.buoyancy_field(YY, ZZ, db, Ly, Htrop, N2=N2)
    u_field = balance.thermal_wind_u(YY, ZZ, db, Ly, Htrop, f)
    return Htrop, YY, ZZ, b_field, u_field


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### The temperature field and the wind it implies

        Left: the buoyancy (temperature) field, with isotherms sloping most
        steeply through the frontal zone. Right: the geostrophic wind $u_g$
        reconstructed purely from that field's meridional gradient — no
        wind was specified anywhere.

        **How to read the pair:** wherever isotherms are *packed* in $y$
        (large $|b_y|$), the wind changes rapidly with height; where they
        are flat, the shear vanishes. The wind panel is, in a precise sense,
        the *vertical integral* of the isotherm slope in the buoyancy panel.
        Forecasters use this reading constantly — a tight thermal gradient
        on an upper-air chart *implies* a jet above it, before any wind
        observation is consulted.
        """
    )
    return


@app.cell(hide_code=True)
def _(YY, ZZ, b_field, np, plotting, plt, u_field):
    fig1, axs1 = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    _vm = max(float(np.abs(b_field).max()), 1e-12)
    _im0 = axs1[0].pcolormesh(YY, ZZ, b_field, cmap=plotting.DIV, vmin=-_vm, vmax=_vm, shading="auto")
    axs1[0].contour(YY, ZZ, b_field, levels=12, colors="k", linewidths=0.4)
    fig1.colorbar(_im0, ax=axs1[0], shrink=0.85)
    axs1[0].set_title("buoyancy $b(y,z)$")
    axs1[0].set_xlabel("$y$"); axs1[0].set_ylabel("$z$")

    _vm2 = max(float(np.abs(u_field).max()), 1e-12)
    _im1 = axs1[1].pcolormesh(YY, ZZ, u_field, cmap=plotting.DIV, vmin=-_vm2, vmax=_vm2, shading="auto")
    axs1[1].contour(YY, ZZ, u_field, levels=12, colors="k", linewidths=0.4)
    fig1.colorbar(_im1, ax=axs1[1], shrink=0.85)
    axs1[1].set_title("geostrophic wind $u_g(y,z)$")
    axs1[1].set_xlabel("$y$"); axs1[1].set_ylabel("$z$")
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Where does the jet core sit?

        A vertical profile of $u_g$ straight through the front ($y=0$, where
        the meridional gradient is strongest). The dashed line marks
        $z=H_{trop}$ — where the buoyancy front's sign reverses.

        **Why the maximum is exactly there:** below $H_{trop}$ the shear
        $u_{g,z}=-b_y/f$ is positive (cold air poleward), so $u_g$ grows
        with height; above it $b_y$ flips sign, the shear turns negative,
        and $u_g$ decays. The wind maximum sits precisely where the
        *horizontal temperature gradient reverses* — which is the actual
        reason the observed jet stream cores at the tropopause
        ($\sim$200 hPa, $\sim$11 km), not somewhere in the mid-troposphere
        where the wind itself might seem "busiest".
        """
    )
    return


@app.cell(hide_code=True)
def _(Htrop, Ly_ui, balance, db_ui, f_ui, np, plt):
    _z = np.linspace(1e-3, 2 * Htrop, 400)
    _u_center = balance.thermal_wind_u(0.0, _z, db_ui.value, Ly_ui.value, Htrop, f_ui.value)

    fig2, ax2 = plt.subplots(figsize=(5, 4.5), constrained_layout=True)
    ax2.plot(_u_center, _z, lw=2)
    ax2.axhline(Htrop, color="crimson", ls="--", lw=1.2, label="$z=H_{trop}$")
    ax2.set_xlabel("$u_g(y{=}0,z)$"); ax2.set_ylabel("$z$")
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    _Ro = float(np.abs(_u_center).max()) / (f_ui.value * Ly_ui.value)
    ax2.set_title(f"jet core at $y=0$  ·  $Ro\\approx{_Ro:.2f}$")
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Isolate the effect of $N^2$.** Move the $N^2$ slider from 0 to 1:
          the isotherms in the left panel tilt more, but confirm the $u_g$
          panel and the jet-core plot are completely unchanged. *Why it must
          be so:* $N^2z$ contributes to $b_z$ but not to $b_y$, and thermal
          wind reads only $b_y$. Stratification controls how *stable* the
          fluid is (Chs. 11–12), not how it is *balanced*.
        - **Narrow the front.** Decrease $L_y$: does the jet get faster or
          slower? *Reasoning to check yourself against:* the shear is
          $-b_y/f\sim\Delta b/(fL_y)$, and the integration depth
          $\sim H_{trop}$ is unchanged, so the peak wind should scale like
          $\Delta b\,H_{trop}/(fL_y)$ — halving $L_y$ should roughly double
          the jet.
        - **Raise the tropopause.** Increase $H_{trop}$: the core should
          track $z=H_{trop}$ exactly (it is pinned to the sign reversal of
          $b_y$, nothing else), and the peak speed should grow $\propto
          H_{trop}$ — a deeper column of one-signed shear to integrate over.
        - **Check the balance.** The title of the jet-core plot reports an
          estimated $Ro=U_{jet}/(fL_y)$. Push $f$ down until $Ro$ approaches
          1. Nothing "breaks" in this notebook — the formulas are happy to
          evaluate — but the *premise* fails: at $Ro\sim1$ the neglected
          accelerations are as large as the terms we kept, so the computed
          $u_g$ is no longer a trustworthy estimate of the actual wind.
          Diagnostic relations fail silently; knowing their domain of
          validity is on you.

        ### What you should have seen

        A purely diagnostic calculation — no time integration anywhere — that
        still produces the single most important structure in
        midlatitude dynamics: a jet whose core sits exactly where the
        meridional temperature gradient changes sign, with speed set by how
        sharp the front is ($1/L_y$) and how strong the rotation is ($1/f$).
        Changing $N^2$ changes the picture but not the wind, because thermal
        wind only cares about the *horizontal* buoyancy gradient.

        **Where this goes next.** The balance you just used diagnostically
        becomes *dynamics* in the following chapters: Ch. 6 asks how a fluid
        that starts *out* of geostrophic balance gets into it (adjustment);
        Chs. 7–8 evolve balanced flow in time via potential vorticity; and
        in Ch. 16 the thermal-wind shear you built here — a front in
        balance — turns out to store available potential energy that
        baroclinic instability converts into weather. The jet you drew is
        stable in this chapter only because nothing here is allowed to move.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
