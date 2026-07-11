import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 4 — Scaling & the Dimensionless Numbers

        **Physical question.** Ch. 2 computed five numbers from raw physical
        inputs. Two of them — $Ro$ and $Bu$ — turn out to do almost all the
        work of predicting which of this book's reduced models applies at a
        given scale. After this chapter you should be able to place ANY
        flow on the $(Ro,Bu)$ plane and read off, without further
        calculation, whether it's unbalanced, baroclinically active, or
        columnar/barotropic — and see exactly where every earlier chapter's
        example flow actually sits.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        Ch. 2 nondimensionalized the full equations and found that a
        handful of coefficients decide which terms survive. Two of those
        coefficients are special: $Ro$ tells you whether a *time-mean
        balance* exists at all, and $Bu=(L_R/L)^2$ tells you, given that it
        does, whether the balance is controlled by *stratification* or by
        *rotation*. Together they define three broad regimes (symbols as in
        [NOTATION](../notation)):

        $$Ro=\frac{U}{f_0L}\ \gtrless\ 1,
          \qquad
          Bu=\left(\frac{L_R}{L}\right)^2\ \gtrless\ 1,
          \qquad L_R=\frac{NH}{f_0}.$$

        - **$Ro>1$ — unbalanced.** Acceleration is as large as Coriolis; no
          steady balance exists to expand around. Inertia-gravity waves
          (Ch. 6) dominate; concepts like PV invertibility (Ch. 7) still
          apply in principle but aren't the efficient description.
        - **$Ro<1,\ Bu>1$ ($L<L_R$) — QG, baroclinic.** The flow is smaller
          than the deformation radius: stratification actively resists
          vertical motion, and baroclinic structure (vertical shear, sloping
          density surfaces) is dynamically central. This is Ch. 16's
          territory, and it is where baroclinic instability operates most
          efficiently — not coincidentally, since $L\sim L_R$ (the $Bu\sim1$
          boundary) is exactly the Eady problem's most-unstable scale.
        - **$Ro<1,\ Bu<1$ ($L>L_R$) — QG, barotropic.** The flow is larger
          than the deformation radius: the stretching term in QGPV
          (Ch. 8) becomes negligible relative to the barotropic vorticity
          term, and the flow tends toward vertically rigid, columnar
          structure — Ch. 3's Taylor-Proudman limit, realized dynamically
          rather than by literally setting buoyancy to zero.

        This is a coarse, three-region caricature of a genuinely continuous
        transition (Vallis, *AOFD* 2nd ed., the Ro-Bu/Ro-Fr regime diagrams
        in ch. 5) — real flows cross these boundaries smoothly, and $Bu\sim1$
        is a region of enhanced activity, not a wall.
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
    from matplotlib.colors import ListedColormap
    from gfdlib import scaling
    return ListedColormap, np, plt, scaling


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        `gfdlib.scaling.classify_regime` is pure algebra (two comparisons),
        cheap enough to evaluate on a dense grid and render as a filled
        region map — fully reactive, no Run button. The marker sliders set
        $(\log_{10}Ro,\log_{10}Bu)$ directly (a proxy for "clicking" a point
        on the map); the small labeled dots are typical $(Ro,Bu)$ values for
        five earlier chapters' example flows, placed by hand from each
        chapter's stated parameters.
        """
    )
    return


@app.cell(hide_code=True)
def _(np, scaling, ListedColormap, plt):
    # --- the region map itself (fixed resolution, always drawn) ------------
    _logRo = np.linspace(-2.0, 1.3, 400)
    _logBu = np.linspace(-2.0, 2.0, 400)
    _LOGRO, _LOGBU = np.meshgrid(_logRo, _logBu, indexing="ij")
    _labels = scaling.classify_regime(10.0 ** _LOGRO, 10.0 ** _LOGBU)
    _cmap = ListedColormap(["#DD8452", "#4C72B0", "#55A868"])

    def draw_map(ax):
        im = ax.pcolormesh(_LOGRO, _LOGBU, _labels, cmap=_cmap, vmin=-0.5, vmax=2.5, shading="auto")
        ax.axvline(0.0, color="k", lw=1.2)
        ax.axhline(0.0, color="k", lw=1.2)
        ax.set_xlabel("$\\log_{10}Ro$"); ax.set_ylabel("$\\log_{10}Bu$")
        return im

    return (draw_map,)


@app.cell(hide_code=True)
def _(mo):
    logRo_ui = mo.ui.slider(-2.0, 1.3, step=0.05, value=-1.0,
                            label="$\\log_{10}Ro$", show_value=True)
    logBu_ui = mo.ui.slider(-2.0, 2.0, step=0.05, value=0.0,
                            label="$\\log_{10}Bu$", show_value=True)
    mo.vstack([mo.hstack([logRo_ui, logBu_ui], justify="start")])
    return logBu_ui, logRo_ui


@app.cell(hide_code=True)
def _(draw_map, logBu_ui, logRo_ui, mo, np, plt, scaling):
    # --- the marker point + reference chapter examples ----------------------
    _refs = [
        ("Ch.5 synoptic front",       -1.0,  0.0),
        ("Ch.16 baroclinic eddy",     -1.3,  0.05),
        ("Ch.18 2D turbulence/jets",  -1.7, -1.3),
        ("Ch.6 large-scale surge",    -1.3, -1.7),
        ("Ch.6 small-scale burst",    -0.5,  1.3),
    ]

    fig1, ax1 = plt.subplots(figsize=(7.5, 6), constrained_layout=True)
    im1 = draw_map(ax1)
    cbar1 = fig1.colorbar(im1, ax=ax1, ticks=[0, 1, 2])
    cbar1.ax.set_yticklabels(["unbalanced", "QG, baroclinic", "QG, barotropic"])

    for _name, _lr, _lb in _refs:
        ax1.plot(_lr, _lb, "o", color="white", markeredgecolor="k", ms=8, zorder=4)
        ax1.annotate(_name, (_lr, _lb), textcoords="offset points", xytext=(6, 6),
                     fontsize=7.5, color="k",
                     bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))

    _Ro, _Bu = 10.0 ** logRo_ui.value, 10.0 ** logBu_ui.value
    ax1.plot(logRo_ui.value, logBu_ui.value, "*", color="crimson", ms=20,
             markeredgecolor="k", zorder=5)
    ax1.set_title(f"you are here:  $Ro={_Ro:.2f}$,  $Bu={_Bu:.2f}$")
    mo.vstack([fig1])
    return


@app.cell(hide_code=True)
def _(logBu_ui, logRo_ui, mo, scaling):
    # --- text read-out for the marker point ---------------------------------
    _Ro, _Bu = 10.0 ** logRo_ui.value, 10.0 ** logBu_ui.value
    _label = int(scaling.classify_regime(_Ro, _Bu))
    _msgs = {
        0: ("**Unbalanced.** No steady balance to lean on — the natural "
            "description is Ch. 6's inertia-gravity-wave adjustment "
            "problem, not a balanced-flow chapter."),
        1: ("**QG, baroclinic.** Smaller than the deformation radius — "
            "stratification and vertical structure are central. This is "
            "the Ch. 16 baroclinic-instability regime, and the natural "
            "scale for developing storms and ocean mesoscale eddies."),
        2: ("**QG, barotropic.** Larger than the deformation radius — "
            "rotation dominates, the flow tends columnar and vertically "
            "rigid (Ch. 3's Taylor-Proudman limit), and Ch. 18's "
            "barotropic-turbulence machinery is the right tool."),
    }
    mo.md(_msgs[_label])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Walk the boundary.** Hold $\log_{10}Bu=0$ (exactly $L=L_R$) and
          sweep $\log_{10}Ro$ across 0: confirm the region flips from
          "QG, baroclinic" to "unbalanced" exactly at $Ro=1$, with no
          dependence at all on $Bu$ at that moment — the $Ro=1$ boundary
          is a vertical line for a reason (the classifier checks $Ro$
          first).
        - **Find yourself among the reference dots.** The two Ch. 6 points
          (small-scale burst vs. large-scale surge) sit on opposite sides
          of $Bu=1$ despite being *the same notebook* — only $\sigma/L_R$
          differs between the two presets. Does that match what you saw
          there: the small-scale burst mostly radiating away as gravity
          waves (unbalanced-leaning, high $Bu$), the large-scale surge
          surviving as a balanced remnant (Ch. 8-style dynamics, low-to-
          moderate $Bu$)?
        - **Locate Ch. 20's gyres.** The wind-driven circulation's basin
          interior is deeply in the $Ro\ll1$ regime, but its western
          boundary current is not — friction, not $Ro$, sets that layer's
          width. Where would you guess the boundary-current sliver sits on
          this map, and why might a map built purely from $Ro,Bu$ (which
          says nothing about friction) be the wrong tool for it specifically?

        ### What you should have seen

        Every notebook you've run so far, however different they looked,
        occupies a definite, predictable position on this one plane — and
        that position alone tells you which chapter's machinery is the
        natural description. The map is a caricature (three flat regions
        standing in for a smooth transition), but the caricature is honest:
        $Ro$ really does gate whether balance exists at all, and $Bu$ really
        does gate whether stratification or rotation controls the balance
        once it does. This is the last chapter of Part I; every remaining
        chapter in the book is, in this map's language, an exploration of
        one specific point (or one specific boundary) on it.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
