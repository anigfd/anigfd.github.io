import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 17 — Symmetric, Inertial & Kelvin-Helmholtz Instability

        **Physical question.** Ch. 15 and Ch. 16 covered the two mechanisms
        that dominate the large-scale storm track and ocean mesoscale. This
        chapter is a **survey** of three more local, faster mechanisms that
        matter at fronts and in stratified shear: inertial instability
        (unbalanced rotation), symmetric instability (a slantwise hybrid of
        the two), and Kelvin-Helmholtz instability (stratified shear,
        rotation-free). All three turn out to be governed by comparing a
        Richardson number to a Rossby number — one map, three mechanisms.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The two numbers being compared

        > **Richardson number**
        > $\;\mathrm{Ri}=\dfrac{N^2}{(\partial U/\partial z)^2}$ —
        > *stratification's stiffness versus the vertical shear trying to
        > overturn it.* **Rossby number** (in this chapter's local form)
        > $\;\mathrm{Ro}=-\dfrac{\partial U/\partial y}{f}$ — *relative
        > vorticity versus planetary.* Every instability in this survey is
        > a statement about where a flow sits in the (Ro, Ri) plane.

        ### Symmetric & inertial instability: two conserved labels

        For a zonal thermal-wind front $U(y,z)$, $b(y,z)$ with
        $f\,\partial U/\partial z=-\partial b/\partial y$, frictionless
        adiabatic motion in the $(y,z)$ plane conserves two parcel labels:
        **absolute momentum** $M=U-fy$ (the rotating frame's version of
        angular momentum) and buoyancy $b$. The stability question is then
        pure geometry: displace a parcel and ask whether its conserved
        labels push it back or further away.

        - Displace *vertically*: $b$-conservation resists if $N^2>0$ —
          ordinary gravitational stability.
        - Displace *horizontally*: $M$-conservation resists if the absolute
          vorticity $f(1+\mathrm{Ro})>0$ — a parcel moved poleward carries
          too little $M$ for its new surroundings and gets pulled back,
          exactly like the angular-momentum argument for orbits.
        - Displace *slantwise*, along a path between the $M$-surfaces and
          the $b$-surfaces: **if the $M$-surfaces are tilted flatter than
          the $b$-surfaces**, there is a wedge of paths along which both
          restoring mechanisms *assist* the displacement. (Compare Ch. 16's
          wedge of instability — same geometry, with $M$ playing the role
          the boundary played there.)

        The wedge exists iff the Ertel PV of the front has the opposite
        sign to $f$:

        $$q=(f-\partial U/\partial y)N^2-f(\partial U/\partial z)^2,
          \qquad qf<0 \iff (1+\mathrm{Ro})<\frac{1}{\mathrm{Ri}}$$

        (symbols as in [NOTATION](../notation)). Two limits separate the
        mechanisms: $(1+\mathrm{Ro})<0$ (absolute vorticity changes sign)
        is unstable at *every* $\mathrm{Ri}>0$ — pure **inertial**
        instability, no buoyancy needed. Finite positive $(1+\mathrm{Ro})$
        with small $\mathrm{Ri}$ is **symmetric** instability — a genuine
        hybrid, extracting energy from both the horizontal shear and the
        sloping buoyancy surfaces.

        ### Kelvin–Helmholtz instability

        For an arbitrary (not necessarily thermal-wind-balanced) stratified
        shear flow $U(z)$, $N^2(z)$, no rotation, normal modes
        $\psi'=\phi(z)e^{ik(x-ct)}$ satisfy the **Taylor–Goldstein
        equation**:

        $$(U-c)^2(\phi''-k^2\phi)-U''(U-c)\phi+N^2\phi=0$$

        — structurally Ch. 15's Rayleigh equation with a buoyancy term
        bolted on (set $N^2=0$ and it collapses back). The energy
        bookkeeping behind the famous threshold: exchanging two parcels
        across a shear layer releases kinetic energy
        $\sim\tfrac14(\Delta U)^2$ but costs potential energy
        $\sim N^2(\Delta z)^2$; the release wins when
        $\mathrm{Ri}\lesssim1/4$. The Miles (1961) / Howard (1961) theorem
        hardens the heuristic into mathematics: if
        $\mathrm{Ri}=N^2/(\partial U/\partial z)^2\geq1/4$ *everywhere*,
        the flow is stable — a genuine guarantee, not an empirical rule of
        thumb.
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
    from gfdlib import instability, symmetric
    return ListedColormap, instability, np, plt, symmetric


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        **Ri-Ro map** (Part A): pure algebra (`gfdlib.symmetric.classify`),
        no eigenvalue problem — cheap enough to be fully reactive.
        **Taylor-Goldstein growth rate** (Part B): a genuinely new numerical
        primitive, `gfdlib.instability.taylor_goldstein_growth_rate` —
        derived by keeping the buoyancy perturbation as an independent
        unknown (rather than eliminating it, which gives a harder-to-trust
        quadratic-in-$c$ pencil), producing a natively **linear** $2n\times
        2n$ generalized eigenvalue problem solved with plain
        `numpy.linalg.eigvals`. At $N^2=0$ it reduces exactly to Ch. 15's
        barotropic solver — the strongest available cross-check. It's about
        8x the matrix size (and far more of the runtime) of Ch. 15's
        solver, so Part B is Run-gated rather than fully reactive.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"## Part A — The Ri-Ro stability map")
    return


@app.cell(hide_code=True)
def _(mo):
    Ro_ui = mo.ui.slider(-2.5, 3.0, step=0.05, value=-0.3, label="Ro", show_value=True)
    Ri_ui = mo.ui.slider(-0.5, 3.0, step=0.02, value=0.6, label="Ri", show_value=True)
    mo.vstack([
        mo.hstack([Ro_ui, Ri_ui], justify="start"),
        mo.md(r"> $\mathrm{Ro}=-(\partial U/\partial y)/f$ (cross-front "
              r"relative-vorticity Rossby number); $\mathrm{Ri}=N^2/"
              r"(\partial U/\partial z)^2$ using the thermal-wind shear."),
    ])
    return Ri_ui, Ro_ui


@app.cell(hide_code=True)
def _(ListedColormap, Ri_ui, Ro_ui, np, plt, symmetric):
    _Ro = np.linspace(-2.5, 3.0, 300)
    _Ri = np.linspace(-0.5, 3.0, 300)
    _RO, _RI = np.meshgrid(_Ro, _Ri)
    _labels = symmetric.classify(_RO, _RI)

    _cmap = ListedColormap(["#4C72B0", "#DD8452", "#C44E52", "#8172B2"])
    fig0, ax0 = plt.subplots(figsize=(6.5, 5), constrained_layout=True)
    im = ax0.pcolormesh(_RO, _RI, _labels, cmap=_cmap, vmin=-0.5, vmax=3.5, shading="auto")
    cbar = fig0.colorbar(im, ax=ax0, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["stable", "symmetric", "inertial", "gravitational"])

    _Ro_curve = np.linspace(-0.99, 3.0, 200)
    ax0.plot(_Ro_curve, symmetric.critical_ri(_Ro_curve), color="k", lw=1.5)
    ax0.axvline(-1.0, color="k", lw=1.5)
    ax0.axhline(0.0, color="k", lw=1.5)

    _label = symmetric.classify(Ro_ui.value, Ri_ui.value)
    _names = {0: "stable", 1: "symmetric", 2: "inertial", 3: "gravitational"}
    ax0.plot(Ro_ui.value, Ri_ui.value, "o", color="white", markeredgecolor="black",
             markersize=10, zorder=5)
    ax0.set_xlabel("Ro"); ax0.set_ylabel("Ri")
    ax0.set_title(f"({Ro_ui.value:.2f}, {Ri_ui.value:.2f}) → {_names[int(_label)]}")
    ax0.set_ylim(-0.5, 3.0)
    fig0
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The marginal curve $\mathrm{Ri}_c=1/(1+\mathrm{Ro})$ passes through
        the textbook point $(\mathrm{Ro},\mathrm{Ri})=(0,1)$ — a front with
        zero relative vorticity is symmetrically unstable for
        $\mathrm{Ri}<1$. For $\mathrm{Ro}<-1$ the front is inertially
        unstable **regardless of stratification** — the vertical black line
        marks this unconditional boundary. Below $\mathrm{Ri}=0$ (the
        horizontal black line) the front is simply convectively unstable,
        independent of shear.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Part B — The Kelvin-Helmholtz growth-rate calculator

        The Hazel (1972) profile $U=\tanh z$, $N^2=J\,\mathrm{sech}^2z$ has
        local Richardson number $\mathrm{Ri}(z)=J\cosh^2z$, minimized (at
        $z=0$) at exactly $J$ — so $J$ is simultaneously the bulk **and**
        minimum Richardson number for this profile, and the Miles-Howard
        necessary condition for instability reduces to the single number
        $J<1/4$.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    J_ui = mo.ui.slider(0.0, 1.0, step=0.02, value=0.15, label="J (bulk Ri)", show_value=True)
    J_ui
    return (J_ui,)


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Compute growth-rate curve")
    run_btn
    return (run_btn,)


@app.cell(hide_code=True)
def _(mo):
    get_curve, set_curve = mo.state(None)
    return get_curve, set_curve


@app.cell
def _(J_ui, instability, mo, np, run_btn, set_curve):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Compute growth-rate curve** above to start."), kind="info"
    ))

    _n = 250
    _z = np.linspace(-10, 10, _n + 2)[1:-1]
    _U, _N2 = instability.hazel_profile(_z, J_ui.value)
    _k_values = np.linspace(0.05, 1.0, 40)
    _growth = instability.taylor_goldstein_growth_rate_curve(_z, _U, _N2, _k_values)
    set_curve({"J": J_ui.value, "k": _k_values, "growth": _growth})

    mo.callout(mo.md(f"**Done.** J={J_ui.value:.2f}."), kind="success")
    return


@app.cell(hide_code=True)
def _(get_curve, mo, np, plt):
    curve = get_curve()
    mo.stop(curve is None, mo.callout(
        mo.md("The growth-rate curve will appear here once you press Compute."),
        kind="warn",
    ))

    _growth = curve["growth"]
    _no_interior_peak = _growth.max() > 1e-8 and _growth[-1] > 0.9 * _growth.max()

    fig1, ax1 = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    ax1.plot(curve["k"], _growth, lw=2)
    if _growth.max() > 1e-8 and not _no_interior_peak:
        _kmax = curve["k"][np.argmax(_growth)]
        ax1.axvline(_kmax, color="crimson", ls="--", lw=1, label=f"peak at $k$={_kmax:.2f}")
        ax1.legend(fontsize=9)
    ax1.set_xlabel("$k$"); ax1.set_ylabel("growth rate $k c_i$")
    ax1.set_title(f"Taylor-Goldstein growth rate, J={curve['J']:.2f} "
                  f"({'unstable' if curve['J'] < 0.25 else 'stable — Miles-Howard'})")
    ax1.grid(alpha=0.3)

    _warning = mo.callout(mo.md(
        "**No interior peak found** — this curve is still rising at the "
        "largest $k$ tested, which is the signature of the finite-difference "
        "truncation error described above (worst near the delicate "
        "$\\mathrm{Ri}=1/4$ point), not genuine short-wave destabilization. "
        "Try a $J$ further from $0.25$ for a cleanly-resolved curve."
    ), kind="warn") if _no_interior_peak else mo.md("")

    mo.vstack([fig1, _warning])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Sweep $J$ from $0$ up past $0.25$: growth shrinks smoothly and
        vanishes by $J\approx0.4$. This simple finite-difference solver
        doesn't cleanly resolve the exact knife-edge at $\mathrm{Ri}=1/4$ —
        the equation's critical-layer indicial roots coalesce exactly there,
        a classically delicate point (professional codes use spectral
        methods for this reason) — but the qualitative conclusion, a hard
        cutoff near $\mathrm{Ri}=1/4$, is the rigorous content of the
        Miles-Howard theorem, not a fitted curve.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Compare the two Richardson-number thresholds.** At
          $\mathrm{Ro}=0$, Part A's symmetric-instability threshold is
          $\mathrm{Ri}=1$; Part B's Kelvin-Helmholtz threshold is
          $\mathrm{Ri}=1/4$. Both instabilities pit stratification against
          shear — why should a **rotating**, geostrophically-balanced parcel
          tolerate four times less shear before going unstable than a
          **non-rotating** one? (Hint: think about what a slantwise
          displacement in a rotating front can exploit that a purely
          vertical displacement cannot.)
        - **Find the inertial boundary.** In Part A, fix $\mathrm{Ri}=2$ and
          slide $\mathrm{Ro}$ down through $-1$. Confirm the front becomes
          unstable exactly at $\mathrm{Ro}=-1$, independent of where you set
          $\mathrm{Ri}$.
        - **Push J past 0.25 in Part B** and confirm the peak growth rate
          you measure decreases roughly monotonically, consistent with (but
          not a razor-sharp numerical proof of) the Miles-Howard theorem.

        ### What you should have seen

        Three instabilities, one underlying comparison: is the restoring
        force (stratification, or planetary rotation via absolute momentum
        conservation) strong enough to overcome the destabilizing shear?
        Inertial instability needs no stratification at all (a purely
        rotational effect); Kelvin-Helmholtz needs no rotation at all (a
        purely buoyancy-vs-shear effect); symmetric instability lives in
        between, and is why the Richardson-number threshold for a rotating,
        balanced front ($\mathrm{Ri}<1$) is markedly more permissive than
        for an unbalanced, non-rotating shear layer ($\mathrm{Ri}<1/4$) —
        slantwise motion has an extra degree of freedom (the horizontal
        shear) to draw on that pure vertical overturning does not.

        **Where this goes next (and where these show up).** Symmetric
        instability organizes the banded "slantwise convection" of winter
        storms' frontal zones and drains low-PV water from the wintertime
        Gulf Stream; inertial instability matters most near the equator,
        where $f\to0$ makes $\mathrm{Ro}<-1$ easy to reach;
        Kelvin–Helmholtz billows are *the* mechanism of clear-air turbulence
        and of interior ocean mixing — the microscale end of the mixing
        chain whose large-scale end is Ch. 19's eddy diffusivity, and whose
        integrated effect is the abyssal $\kappa$ of Ch. 21. This closes
        Part V's survey: every instability in the book has been an unstable
        arrangement of the same conserved quantities (PV, $M$, $b$), and
        Part VI now turns to the turbulence they collectively produce.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
