import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 20 — Wind-Driven Ocean Circulation

        **Physical question.** The trade winds and westerlies blow steadily
        over every subtropical ocean basin, curling clockwise (Northern
        Hemisphere) around each gyre. Sverdrup balance alone can't close
        the circulation in a bounded basin — something has to happen at
        the western edge. After this chapter you should be able to explain
        *why* western boundary currents like the Gulf Stream exist (and not
        eastern ones), and compute their width from first principles.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### How wind becomes vorticity: Ekman pumping

        The wind does not push the gyre around directly. Stress drives a
        thin surface **Ekman layer** whose transport is $90°$ to the right
        of the wind (Northern Hemisphere); where the stress *curl* is
        negative (as between the westerlies and the trades), the Ekman
        transports converge and pump water gently **downward** into the
        interior ($w_{Ek}=\nabla\times\boldsymbol\tau/(\rho f)$, at a
        stately ~30 m/*year*). The interior feels the wind only as this
        squashing of its water columns — a vorticity forcing.

        ### Sverdrup balance: the interior in one line

        Steady, linear, vertically-integrated vorticity balance on a
        $\beta$-plane (symbols as in [NOTATION](../notation)):

        $$\beta\frac{\partial\psi}{\partial x}=\frac{1}{\rho H}\nabla\times
          \boldsymbol\tau + \mathcal F,$$

        nondimensionalized on a unit square basin with idealized single-gyre
        forcing $\nabla\times\boldsymbol\tau\propto-\sin(\pi y)$. Read the
        frictionless balance ($\mathcal F=0$) physically: $\beta\psi_x$ is
        $\beta v$ — a column squashed by Ekman pumping must *lose* planetary
        vorticity, so it slides **equatorward** (Ch. 7's $Dq/Dt=0$ yet
        again). The whole subtropical interior drifts slowly south; that is
        the **Sverdrup balance** (1947), and it is beautifully verified in
        the real ocean's interior.

        But integrate it: $\psi_x=-\sin(\pi y)$ is *first-order* in $x$ and
        can satisfy $\psi=0$ at only *one* wall (conventionally the eastern):
        $\psi_{Sv}=(1-x)\sin(\pi y)$. Mass does not close; every parcel
        drifting equatorward must somehow get back poleward. Something else
        has to happen at $x=0$.

        ### Closing the basin: friction, and why the west

        Two classical choices for $\mathcal F$:

        $$\textbf{Stommel (1948):}\ \varepsilon\nabla^2\psi+\psi_x=-\sin(\pi y),
        \qquad
        \textbf{Munk (1950):}\ -\delta^3\nabla^4\psi+\psi_x=-\sin(\pi y),$$

        linear bottom drag and lateral (biharmonic) friction respectively.
        Each introduces a thin boundary layer of width $\varepsilon$
        (Stommel) or $\delta$ (Munk) in which friction competes with
        $\beta$. Why must it sit on the **western** wall? Vorticity
        bookkeeping: the wind injects negative vorticity over the whole
        gyre, steadily, so the boundary layer must be a net *source* of
        positive vorticity. A poleward return current hugging the west wall
        rubs against it in the sense that supplies exactly that; the same
        current on the east wall would *add* negative vorticity, and the
        budget could never balance. (Equivalently: relative to the westward
        Rossby-wave drift of Chs. 8–9, the western wall is "downstream" —
        energy piles up there like water at a dam.) Hence the Gulf Stream
        sits at the American coast and not the European one — an asymmetry
        set by which way the planet spins, and nothing else.
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
    from gfdlib import plotting, circulation
    return circulation, np, plotting, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        A genuinely new numerical primitive for this book: every previous
        chapter's domain was either doubly-periodic (spectral) or a channel
        with two boundaries; this is a fully **bounded rectangular basin**
        with four walls. `gfdlib.circulation` solves both problems two
        independent ways: an **exact analytic solution** (separating
        $\psi=f(x)\sin(\pi y)$ reduces the PDE to a linear ODE for $f(x)$,
        solved exactly via its characteristic polynomial's roots — quadratic
        formula for Stommel, `numpy.roots` for Munk's quartic — and a small
        `numpy.linalg.solve` for the boundary-condition constants), and a
        **general 2D finite-difference direct solve** (Kronecker-sum
        Laplacian/derivative operators, `numpy.linalg.solve` — no SciPy) that
        works for arbitrary forcing, cross-validated against the analytic
        solution with clean 2nd-order convergence. Munk's biharmonic operator
        is handled by keeping vorticity $\zeta=\nabla^2\psi$ as an
        independent field (the same pattern as ch16/17/19) rather than a
        hand-derived 4th-order stencil, closed with Thom's (1933)
        wall-vorticity formula for no-slip on the east/west walls.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"## Part A — Dial the friction, watch the boundary current")
    return


@app.cell(hide_code=True)
def _(mo):
    friction_ui = mo.ui.dropdown(
        options={"Stommel (bottom drag)": "stommel", "Munk (lateral friction)": "munk"},
        value="Stommel (bottom drag)", label="friction type",
    )
    eps_ui = mo.ui.slider(0.01, 0.15, step=0.005, value=0.05,
                          label="Stommel $\\varepsilon$", show_value=True)
    delta_ui = mo.ui.slider(0.03, 0.25, step=0.005, value=0.1,
                            label="Munk $\\delta$", show_value=True)
    n_ui = mo.ui.dropdown(options={"40": 40, "60": 60}, value="40", label="resolution $n$")
    mo.vstack([
        mo.hstack([friction_ui, n_ui], justify="start"),
        mo.hstack([eps_ui, delta_ui], justify="start"),
        mo.md(r"> Only the slider matching your chosen friction type "
              r"affects the solution below. $n=60$ is more accurate but "
              r"noticeably slower (a full basin-wide linear solve, not a "
              r"time-stepped simulation)."),
    ])
    return delta_ui, eps_ui, friction_ui, n_ui


@app.cell(hide_code=True)
def _(circulation, delta_ui, eps_ui, friction_ui, n_ui, np):
    _n = n_ui.value
    if friction_ui.value == "stommel":
        _eps = eps_ui.value
        psi_fd, x_full, y_full = circulation.solve_stommel_fd(_n, _n, _eps)
        _f_exact = circulation.stommel_f(x_full, _eps)
        param_label = f"$\\varepsilon$={_eps:.3f}"
    else:
        _delta = delta_ui.value
        psi_fd, x_full, y_full = circulation.solve_munk_fd(_n, _n, _delta)
        _f_exact = circulation.munk_f(x_full, _delta)
        param_label = f"$\\delta$={_delta:.3f}"

    f_sverdrup = 1.0 - x_full
    result = {
        "psi": psi_fd, "x": x_full, "y": y_full,
        "f_exact": _f_exact, "f_sverdrup": f_sverdrup,
        "friction": friction_ui.value, "param_label": param_label,
    }
    return (result,)


@app.cell(hide_code=True)
def _(np, plotting, plt, result):
    class _G:   # minimal shim: plotting.field only needs the domain length L
        L = 1.0

    fig1, axs1 = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True,
                              width_ratios=[1.15, 1])
    plotting.field(axs1[0], result["psi"].T, grid=_G(), signed=False,
                   title=f"$\\psi(x,y)$ — {result['friction']}, {result['param_label']}")
    axs1[0].set_xlabel("$x$ (west $\\to$ east)"); axs1[0].set_ylabel("$y$")

    _f_fd = result["psi"][len(result["y"]) // 2, :]   # psi(x, y=0.5) = f(x)
    axs1[1].plot(result["x"], result["f_sverdrup"], "k--", lw=1.2, label="Sverdrup (no friction)")
    axs1[1].plot(result["x"], result["f_exact"], lw=1, color="gray", alpha=0.6, label="exact analytic")
    axs1[1].plot(result["x"], _f_fd, lw=2, color="crimson", label="finite-difference")
    axs1[1].set_xlabel("$x$"); axs1[1].set_ylabel("$f(x)=\\psi(x,0.5)$")
    axs1[1].set_title("west$\\to$east profile through gyre center")
    axs1[1].legend(fontsize=8); axs1[1].grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Notice: **only the western third or less of the basin ever departs
        from the Sverdrup (dashed) curve.** That departure — a thin band of
        intense flow hugging $x=0$ — *is* the western boundary current. The
        eastern two-thirds of every real subtropical gyre is a broad, slow,
        almost frictionless Sverdrup flow; the Gulf Stream, Kuroshio, and
        Agulhas current are all the SAME thin sliver, concentrated at the
        western edge because that's the only place the (Sverdrup + $\beta$)
        balance can be locally overwhelmed by friction without changing sign
        everywhere else.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Both fluids

        - **Ocean (the primary case):** every subtropical ocean gyre has
          exactly this asymmetry — broad Sverdrup interior, thin intense
          western boundary current (Gulf Stream, Kuroshio, Agulhas,
          Brazil Current) — set by the same $\beta$-plane vorticity balance
          dialed above, typically closer to the Munk (lateral-friction,
          $\delta\sim$tens of km) than the Stommel limit for the real ocean.
        - **Atmosphere:** the same linear vorticity balance (forcing curl
          balanced by $\beta$ and friction) describes the depth-averaged,
          time-mean response to Ekman pumping under the subtropical highs —
          but the atmosphere's friction (turbulent boundary-layer drag) and
          much larger deformation radius put it in a different part of
          parameter space, without the ocean's dramatic western
          intensification; the *mechanism* is the same equation, the
          *appearance* is not, a useful reminder that the same balance can
          look very different depending on which term actually dominates.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Shrink the friction.** Push $\varepsilon$ or $\delta$ toward
          their smallest slider values: the boundary current gets thinner
          and the interior hugs the Sverdrup curve ever more closely — the
          textbook singular-perturbation limit.
        - **Compare boundary-layer shapes.** At similar boundary-layer
          *widths*, Stommel's profile approaches the interior
          monotonically (pure exponential decay) while Munk's can
          overshoot slightly before settling — a real structural
          difference between a 2nd-order and a 4th-order boundary layer,
          visible directly in the profile panel.
        - **Check the resolution matters less than the physics.** Switch
          $n$ from 40 to 60 at fixed friction: the finite-difference curve
          should barely move, since it's already converged onto the exact
          analytic curve — a working direct confirmation that the discrete
          solver is trustworthy, not just fast.

        ### What you should have seen

        A steady linear balance that *cannot* close itself everywhere with
        friction turned off (Sverdrup), forced to concentrate all of its
        friction into a thin sliver at exactly one wall — and which wall
        depends on the sign of $\beta$, not on which side of the basin you
        pick. That geometric asymmetry, not any special property of the
        Gulf Stream itself, is the entire reason western boundary currents
        exist and eastern ones don't.

        **Where this goes next.** This linear, steady gyre is the skeleton;
        the real Gulf Stream is its unstable, eddying descendant — run the
        boundary current through Ch. 16's baroclinic instability and
        Ch. 18's turbulence and you get the observed meanders, rings, and
        recirculations. Ch. 21 completes the circulation picture with the
        part the wind cannot do: the buoyancy-driven overturning that
        ventilates the deep ocean beneath these wind-driven gyres.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
