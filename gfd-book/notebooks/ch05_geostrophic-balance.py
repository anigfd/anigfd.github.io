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

        Start from the horizontal momentum and hydrostatic equations
        (symbols as in [NOTATION](../notation)), with $\phi=p/\rho_0$ and
        $b=-g\rho'/\rho_0$:

        $$\frac{Du}{Dt}-fv=-\phi_x,\qquad \frac{Dv}{Dt}+fu=-\phi_y,\qquad
          \phi_z=b.$$

        At low Rossby number $Ro=U/(fL)\ll1$, the acceleration terms
        $Du/Dt,Dv/Dt$ are $O(Ro)$ smaller than the Coriolis and pressure
        terms and drop out at leading order, leaving **geostrophic balance**:

        $$u_g=-\frac{\phi_y}{f},\qquad v_g=\frac{\phi_x}{f}.$$

        This is exactly the streamfunction relation from NOTATION with
        $\psi=\phi/f$ — geostrophic flow is automatically nondivergent, and
        $\phi/f$ *is* its streamfunction. Differentiating $u_g$ in $z$ and
        using hydrostatic balance to swap in $b$ gives the **thermal wind**
        relation:

        $$f\,u_{g,z}=-b_y,\qquad f\,v_{g,z}=b_x.$$

        Vertical shear of the geostrophic wind is fixed entirely by the
        *horizontal* buoyancy gradient. This notebook reconstructs $u_g(y,z)$
        by integrating that relation upward from an idealized frontal
        buoyancy field
        $$b(y,z)=N^2z+\Delta b\,\tanh(y/L_y)\cos\!\Big(\frac{\pi z}{2H_{trop}}\Big),$$
        whose sign reverses at $z=H_{trop}$ (an idealized tropopause).
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
          panel and the jet-core plot are completely unchanged.
        - **Narrow the front.** Decrease $L_y$: does the jet get faster or
          slower? Relate this to the $1/L_y$ in `dbdy_frontal`.
        - **Raise the tropopause.** Increase $H_{trop}$: does the jet core
          move, and does its peak speed change?
        - **Check the balance.** The title of the jet-core plot reports an
          estimated $Ro=U_{jet}/(fL_y)$. Push $f$ down until $Ro$ approaches
          1 — at that point, is geostrophic balance still a good
          approximation?

        ### What you should have seen

        A purely diagnostic calculation — no time integration anywhere — that
        still produces the single most important structure in
        midlatitude dynamics: a jet whose core sits exactly where the
        meridional temperature gradient changes sign, with speed set by how
        sharp the front is ($1/L_y$) and how strong the rotation is ($1/f$).
        Changing $N^2$ changes the picture but not the wind, because thermal
        wind only cares about the *horizontal* buoyancy gradient. This is the
        same balance, expressed as $\psi=\phi/f$, that every later chapter's
        geostrophic flows build on.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
