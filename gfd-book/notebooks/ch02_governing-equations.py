import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 2 — Governing Equations & the Boussinesq Approximation

        **Physical question.** The full equations of motion for a rotating,
        stratified fluid have a dozen terms — acceleration, Coriolis,
        pressure gradient, buoyancy, friction, diffusion. At any *given*
        scale, most of them are negligible. After this chapter you should be
        able to take a length scale, a velocity scale, and a latitude, and
        say *in advance* — without solving anything — which terms survive
        and which reduced equation set (this book's actual subject) applies.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The full system (Boussinesq)

        For a fluid whose density varies only weakly about a reference
        $\rho_0$ (the Boussinesq approximation — valid whenever the depth of
        the fluid is much less than a density scale height, true for both
        the ocean and, less obviously, most of the troposphere), mass,
        momentum, and buoyancy obey (symbols as in
        [NOTATION](../notation)):

        $$\nabla\cdot\mathbf u=0,\qquad
          \frac{D\mathbf u}{Dt}+f\hat{\mathbf z}\times\mathbf u=-\nabla\phi+b\hat{\mathbf z}+\nu\nabla^2\mathbf u,
          \qquad \frac{Db}{Dt}=\kappa\nabla^2 b.$$

        Every later chapter in this book is what remains of this system
        after dropping the terms that don't matter at that chapter's scale.

        ### Nondimensionalize, and the coefficients ARE the answer

        Scale length by $L$, velocity by $U$, time advectively by $L/U$,
        and pressure by $f_0UL$ (chosen so the Coriolis and pressure terms
        both come out $O(1)$ — the geostrophic scaling). The horizontal
        momentum equation becomes (hats denote nondimensional variables,
        $Ro=U/(f_0L)$):

        $$Ro\,\frac{D\hat{\mathbf u}}{D\hat t}+\hat{\mathbf z}\times\hat{\mathbf u}
          =-\nabla\hat\phi+\cdots$$

        Every term in the *original* equation now carries an explicit
        dimensionless coefficient in front of it — and that coefficient
        *is* the answer to "does this term matter here?" A coefficient
        $\ll1$ means: drop that term, at leading order, and you get one of
        this book's balanced or wave equations for free. This is the whole
        method of the book, stated once, explicitly, so every later
        derivation you already saw (Ch. 5's geostrophy from $Ro\to0$,
        Ch. 12's stratified-wave equations, Ch. 14's Boussinesq
        convection) can be recognized as an instance of it.

        | Number | Definition | Multiplies |
        |---|---|---|
        | $Ro=U/(f_0L)$ | acceleration / Coriolis | $D\mathbf u/Dt$ |
        | $Ek=\nu/(f_0H^2)$ | friction / Coriolis (vertical/boundary-layer scale $H$) | $\nu\partial^2\mathbf u/\partial z^2$ |
        | $Fr=U/(NH)$ | inertia / buoyancy restoring | vertical motion vs. stratification |
        | $Re=UL/\nu$ | inertia / viscosity | $\nu\nabla^2\mathbf u$ |
        | $Bu=(L_R/L)^2$, $L_R=NH/f_0$ | (rotation length / your length)$^2$ | relates $Ro$ and $Fr$: $Bu=(Fr/Ro)^2$ |

        `gfdlib.scaling` computes each of these from raw physical inputs
        ($U,L,H,N,\nu$, latitude) — the calculator below.
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
    from gfdlib import scaling
    return np, plt, scaling


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Purely diagnostic — every number below is a closed-form expression
        (`gfdlib.scaling`), evaluated instantly from the sliders. No
        discretization, no time-stepping.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    preset_ui = mo.ui.dropdown(
        options={
            "Custom": "custom",
            "Midlatitude synoptic weather": "synoptic",
            "Ocean mesoscale eddy": "eddy",
            "Convective thunderstorm": "storm",
            "Bathtub vortex": "bathtub",
        },
        value="Midlatitude synoptic weather",
        label="Preset",
    )
    preset_ui
    return (preset_ui,)


@app.cell(hide_code=True)
def _(mo, preset_ui):
    _presets = {
        "custom":   dict(U=10.0,  L=1.0e6, H=1.0e4, N=1.0e-2,  nu=1.0e1,  lat=45.0),
        "synoptic": dict(U=10.0,  L=1.0e6, H=1.0e4, N=1.0e-2,  nu=1.0e1,  lat=45.0),
        "eddy":     dict(U=0.1,   L=1.0e5, H=1.0e3, N=5.0e-3,  nu=1.0e1,  lat=30.0),
        "storm":    dict(U=15.0,  L=1.0e4, H=1.0e4, N=1.0e-2,  nu=1.0e1,  lat=35.0),
        "bathtub":  dict(U=0.1,   L=0.2,   H=0.1,   N=0.0,     nu=1.0e-6, lat=45.0),
    }
    _p = _presets[preset_ui.value]

    U_ui = mo.ui.slider(0.001, 50.0, step=0.001, value=_p["U"],
                        label="$U$ velocity scale (m/s)", show_value=True)
    L_ui = mo.ui.slider(1.0e2, 1.0e7, step=1.0e2, value=_p["L"],
                        label="$L$ horizontal length scale (m)", show_value=True)
    H_ui = mo.ui.slider(1.0e-2, 2.0e4, step=1.0, value=_p["H"],
                        label="$H$ vertical length scale (m)", show_value=True)
    N_ui = mo.ui.slider(0.0, 2.0e-2, step=1.0e-4, value=_p["N"],
                        label="$N$ buoyancy frequency (1/s, 0 = unstratified)",
                        show_value=True)
    nu_ui = mo.ui.slider(1.0e-6, 1.0e2, step=1.0e-6, value=_p["nu"],
                         label="$\\nu$ viscosity (m²/s — turbulent/eddy unless you set it molecular)",
                         show_value=True)
    lat_ui = mo.ui.slider(-89.0, 89.0, step=1.0, value=_p["lat"],
                          label="latitude (deg)", show_value=True)
    mo.vstack([
        mo.hstack([U_ui, L_ui], justify="start"),
        mo.hstack([H_ui, N_ui], justify="start"),
        mo.hstack([nu_ui, lat_ui], justify="start"),
    ])
    return H_ui, L_ui, N_ui, U_ui, lat_ui, nu_ui


@app.cell(hide_code=True)
def _(H_ui, L_ui, N_ui, U_ui, lat_ui, np, nu_ui, scaling):
    # --- compute every dimensionless number (reactive) ----------------------
    U, L, H, N, nu, lat = U_ui.value, L_ui.value, H_ui.value, N_ui.value, nu_ui.value, lat_ui.value
    f0 = scaling.coriolis_parameter(lat)
    f0_safe = f0 if abs(f0) > 1e-10 else np.sign(f0 + 1e-30) * 1e-10

    Ro = scaling.rossby_number(U, f0_safe, L)
    Ek = scaling.ekman_number(nu, f0_safe, H)
    Fr = scaling.froude_number(U, max(N, 1e-8), H)
    Re = scaling.reynolds_number(U, L, nu)
    Bu = scaling.burger_number(max(N, 1e-8), H, f0_safe, L)
    Lr = scaling.deformation_radius(max(N, 1e-8), H, f0_safe)
    numbers = {"Ro": Ro, "Ek": Ek, "Fr": Fr, "Re": Re, "Bu": Bu}
    return Bu, Ek, Fr, Lr, Re, Ro, f0, numbers


@app.cell(hide_code=True)
def _(Bu, Ek, Fr, Lr, Re, Ro, f0, mo, np, plt):
    # --- diagnostic 1: log-scale bar chart of the five numbers -------------
    fig1, ax1 = plt.subplots(figsize=(7, 3.2), constrained_layout=True)
    _names = ["Ro", "Ek", "Fr", "Re", "Bu"]
    _vals = [Ro, Ek, Fr, Re, Bu]
    _colors = ["#2563eb" if v < 1 else "#dc2626" for v in _vals]
    ax1.barh(_names, np.log10(np.clip(np.abs(_vals), 1e-12, 1e12)), color=_colors)
    ax1.axvline(0.0, color="k", lw=1.2, label="=1")
    ax1.set_xlabel("$\\log_{10}$(value)")
    ax1.set_title(f"$f_0={f0:.2e}\\,$s$^{{-1}}$,  $L_R={Lr/1e3:.1f}\\,$km")
    ax1.legend(fontsize=8, loc="lower right")
    ax1.grid(alpha=0.3, axis="x")
    mo.vstack([
        fig1,
        mo.md("Blue bars ($<1$, left of the line): that physics is "
              "**subdominant** at this scale — drop the term it multiplies. "
              "Red bars ($>1$): that physics **matters**, keep the term."),
    ])
    return


@app.cell(hide_code=True)
def _(Bu, Ek, Fr, Re, Ro, mo):
    # --- diagnostic 2: plain-language regime read-out -----------------------
    _lines = []
    if Ro < 0.1:
        _lines.append("**Ro ≪ 1** — acceleration is negligible next to Coriolis: "
                       "the flow is *geostrophically balanced* at leading order (Ch. 5).")
    elif Ro < 1:
        _lines.append("**Ro = O(0.1–1)** — acceleration is a correction to geostrophy, "
                       "not negligible: this is quasi-geostrophic territory (Ch. 8).")
    else:
        _lines.append("**Ro ≳ 1** — acceleration is AS LARGE AS Coriolis: there is no "
                       "steady balance to expand around. Expect inertia-gravity waves "
                       "and unbalanced motion (Ch. 6), not balanced dynamics.")

    if Ek < 0.01:
        _lines.append("**Ek ≪ 1** — friction is negligible in the interior; it only "
                       "matters in thin boundary layers (the classical Ekman layer).")
    else:
        _lines.append("**Ek = O(1) or larger** — friction competes directly with "
                       "rotation at this scale (Ch. 20's Stommel gyre is exactly this "
                       "limit, applied to an entire basin).")

    if Fr < 0.3:
        _lines.append("**Fr ≪ 1** — strongly stratified: vertical motion is strongly "
                       "suppressed, the flow is quasi-horizontal and layered (Ch. 11).")
    elif Fr > 3:
        _lines.append("**Fr ≫ 1** — stratification barely resists vertical motion at "
                       "this scale; treat the fluid as nearly unstratified.")
    else:
        _lines.append("**Fr = O(1)** — stratification and inertia are comparably "
                       "important: internal waves (Ch. 12) are the natural language.")

    if Re > 1e4:
        _lines.append("**Re ≫ 1** — inertia utterly dominates molecular viscosity: "
                       "this flow is turbulent, and any $\\nu$ used in a model of it "
                       "is really an eddy/subgrid viscosity, not the molecular value.")

    if Bu > 3:
        _lines.append("**Bu ≫ 1 ($L\\ll L_R$)** — smaller than the deformation radius: "
                       "stratification (baroclinic structure, Ch. 16) controls the "
                       "response.")
    elif Bu < 0.3:
        _lines.append("**Bu ≪ 1 ($L\\gg L_R$)** — larger than the deformation radius: "
                       "rotation dominates, flow tends toward columnar, "
                       "Taylor-Proudman-like structure (Ch. 3).")
    else:
        _lines.append("**Bu = O(1) ($L\\sim L_R$)** — rotation and stratification are "
                       "in direct competition: this is exactly the scale at which "
                       "baroclinic instability (Ch. 16) operates most efficiently.")

    mo.md("### What this scale means\n\n" + "\n\n".join("- " + l for l in _lines))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Cross from balanced to unbalanced.** Start at "Midlatitude
          synoptic weather" ($Ro\approx0.1$) and shrink $L$ toward the
          "Convective thunderstorm" preset's scale: watch $Ro$ cross 1.
          At what $L$ (roughly) does the read-out's advice flip from
          "geostrophically balanced" to "no steady balance to expand
          around"? Does that match your intuition that thunderstorms are
          NOT approximately geostrophic, while synoptic weather is?
        - **Push $Bu$ across 1.** Starting from the ocean-eddy preset,
          shrink or grow $L$ until $Bu$ crosses 1 — watch $L_R$ printed
          above the bar chart stay fixed (it only depends on $N,H,f_0$)
          while $Bu$'s regime label flips. $Bu=1$ exactly means $L=L_R$.
        - **Find where $Re$ stops mattering.** Try to make $Re<1$ using the
          slider's viscosity range. Can you, at any oceanographically or
          atmospherically plausible $U,L$? What does that tell you about
          how rare genuinely viscous (non-turbulent) large-scale flows are
          in GFD — and why every "$\nu$" used elsewhere in this book is an
          eddy viscosity standing in for unresolved turbulence, not the
          literal molecular value?
        - **The bathtub vortex.** Switch to that preset (note the *tiny*
          $L,H$ and *molecular* $\nu$) and read $Ro$: it should be
          enormous. This is the standard debunking of "does the bathtub
          drain direction depend on hemisphere" — at bathtub scale,
          rotation is utterly irrelevant to the dynamics.

        ### What you should have seen

        Five numbers, each just a ratio of two physical rates, together
        determine which terms in the full Boussinesq equations survive at
        a given scale — and therefore which reduced model in this book
        applies. Midlatitude weather sits at $Ro\ll1$ (balanced),
        thunderstorm scale does not; ocean eddies sit near $Bu\sim1$
        (baroclinically active); a bathtub vortex sits at $Ro\gg1$
        (rotation irrelevant). None of this requires solving a single
        equation — it is entirely a bookkeeping exercise on the
        coefficients that appear when you nondimensionalize, and it is
        the exercise every chapter after this one has already silently
        performed to arrive at its reduced equation set.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
