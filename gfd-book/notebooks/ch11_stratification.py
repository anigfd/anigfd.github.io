import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 11 — Stratification & the Vertical-Mode Solver

        **Physical question.** Ch. 10 split a two-LAYER fluid exactly into
        barotropic and baroclinic parts. A *continuously* stratified fluid
        admits the same trick, generalized: any $N^2(z)$ profile has its
        own complete set of orthogonal **vertical normal modes**, each
        behaving as an independent single-layer QG system with its own
        deformation radius. After this chapter you should be able to solve
        for those modes given any $N^2(z)$, and see the two-layer model of
        Chs. 10 and 16 for what it really is: the crudest possible
        truncation of this richer theory, keeping only mode 0 and mode 1.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### From layers to a continuum

        Ch. 10's two-layer split, $\psi_{bt}=(\psi_1+\psi_2)/2$,
        $\psi_{bc}=(\psi_1-\psi_2)/2$, decoupled the layer PV equations
        exactly. The continuous analog does the same for a stratified fluid
        that varies smoothly with $z$, and the derivation is worth seeing
        once in full. Take the linear, hydrostatic, Boussinesq equations
        (Ch. 2) and look for solutions in which the *horizontal* structure
        is common to $u,v,\phi$ (they are tied together by geostrophy and
        hydrostatics) while $w$ and $b$ carry their own vertical shapes:

        $$ (u,v,\phi)\propto\Phi(z),\qquad
           w\propto W(z),\qquad b\propto B(z),$$

        each multiplying the same horizontal wave. Three of the equations
        then relate the three vertical structures:

        - **Continuity** $u_x+v_y+w_z=0$: the horizontal divergence has
          $\Phi$'s shape, so $W'\propto\Phi$ — say $W'=\gamma\,\Phi$ for a
          separation constant $\gamma$.
        - **Hydrostatic balance** $\phi_z=b$: $B=\Phi'$ (up to the same
          constants).
        - **Buoyancy conservation** $b_t+wN^2=0$: $B\propto N^2W$, i.e.
          $\Phi'\propto N^2W$, so $W\propto\Phi'/N^2$.

        Substituting the third relation into the first eliminates $W$
        entirely, and the separation constant — call it $1/c^2$, with
        dimensions of (speed)$^{-2}$ — lands in exactly one place
        (symbols as in [NOTATION](../notation)):

        $$\frac{d}{dz}\!\left[\frac{1}{N^2(z)}\frac{d\Phi_n}{dz}\right]
          +\frac{1}{c_n^2}\Phi_n=0,
          \qquad \Phi_n'(0)=\Phi_n'(H)=0.$$

        The boundary conditions are just $w=0$ at the rigid bottom and lid:
        since $W\propto\Phi'/N^2$, zero vertical velocity means zero
        $\Phi'$. Notice what does **not** appear: $f_0$. Rotation enters
        only afterward, when each mode's gravity-wave speed $c_n$ is
        converted to a deformation radius $L_n=c_n/f_0$ — the eigenvalue
        problem itself is a property of the stratification alone.

        This is a genuine Sturm-Liouville eigenvalue problem: it has a
        countable, orthogonal family of solutions $\Phi_0,\Phi_1,\Phi_2,
        \ldots$, and every field in the linear problem — pressure, $u$,
        $v$ — shares this SAME set of vertical shapes. Each mode's speed
        $c_n$ converts to a deformation radius exactly as in shallow water:
        $L_n=c_n/f_0$.

        ### Mode 0 is barotropic, mode 1 is what Chs. 10/16 call "baroclinic"

        $\Phi_0=$ const solves the eigenvalue problem trivially (its
        derivative is zero everywhere, satisfying the equation with
        $c_0=\infty$) — no deformation-radius constraint at all, exactly
        Ch. 3's Taylor-Proudman limit and Ch. 10's barotropic mode. Modes
        $n\geq1$ are genuinely "baroclinic": they change sign at least once
        in $z$ (mode $n$ has $n$ internal zero-crossings — a general
        Sturm-Liouville fact, not special to this problem), each with its
        own finite $c_n$ and $L_n$. **The two-layer model is exactly the
        $n=0,1$ truncation of this theory** — it cannot represent mode 2 or
        higher at all, by construction, since two layers give only two
        degrees of freedom in $z$.

        ### A closed form to check against

        For **constant** $N=N_0$: $\Phi_n(z)=\cos(n\pi z/H)$,
        $c_n=N_0H/(n\pi)$ — verified by direct substitution (see
        `gfdlib.stratification`'s docstring and tests). This is the
        benchmark the numerical solver below is checked against before
        trusting it on a realistic, depth-varying $N^2(z)$.

        ### The WKB estimate: one integral (almost) beats the eigensolver

        For variable $N(z)$ there is a classical approximation that
        requires no matrix at all. Seek a rapidly-oscillating solution
        $\Phi\sim e^{i\theta(z)}$ with slowly-varying local wavenumber
        $\theta'(z)$: substituting into the eigenvalue equation and keeping
        the dominant ($\theta'^2$) terms gives $\theta'(z)=N(z)/c$ — the
        mode oscillates fastest in $z$ exactly where the stratification is
        strongest. Requiring $n$ half-wavelengths to fit between the two
        Neumann boundaries ($\theta(H)-\theta(0)=n\pi$) then quantizes $c$:

        $$c_n\approx\frac{1}{n\pi}\int_0^H N(z)\,dz.$$

        For constant $N$ this reproduces $N_0H/(n\pi)$ *exactly*, and it is
        the formula behind the standard global atlas of deformation radii
        built from hydrographic data (Chelton et al. 1998). But test it
        against this notebook's own solver before trusting it — the result
        is instructive. For the default pycnocline profile, WKB gets modes
        2–5 to within a few percent, **but is off by $\sim$40% for
        mode 1**. That is not a bug; it is the approximation's own premise
        failing exactly where you'd most like to use it: WKB assumed a
        *rapidly oscillating* $\Phi$, and mode 1 — with its single, slow
        sign change — oscillates least of all. The sharper the pycnocline
        (the further from constant $N$), the worse the mode-1 estimate;
        Chelton et al. apply an empirical correction for precisely this
        reason. The honest summary: $c_n\propto1/n$ with the coefficient
        $\int N\,dz/\pi$ is asymptotically exact for high modes and a
        useful first guess — never a substitute — for the gravest one.
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
    from gfdlib import stratification, plotting
    return np, plotting, plt, stratification


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        `gfdlib.stratification.vertical_modes` discretizes the
        Sturm-Liouville problem by the finite-VOLUME method (unequal
        control volumes at the two Neumann boundaries), then symmetrizes
        the resulting generalized eigenproblem into an ordinary one via a
        diagonal congruence transform, solved with plain `numpy.linalg.
        eigh` — no SciPy, and the right tool for this genuinely
        self-adjoint problem (an early, naive ghost-point attempt produced
        a non-symmetric matrix with badly wrong mode SHAPES despite
        nearly-correct eigenvalues; caught only by plotting against the
        constant-$N$ exact solution, not by the eigenvalues alone — see
        the module's docstring). Fully reactive: no time-stepping anywhere
        in this chapter.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    profile_ui = mo.ui.dropdown(
        options={
            "Constant N (validates against cos(n*pi*z/H))": "constant",
            "Idealized pycnocline (realistic ocean-like)": "pycnocline",
        },
        value="Idealized pycnocline (realistic ocean-like)",
        label="$N^2(z)$ profile",
    )
    N0_ui = mo.ui.slider(0.002, 0.02, step=0.001, value=0.01,
                        label="$N_0$ (constant-profile value)", show_value=True)
    Nmin_ui = mo.ui.slider(0.0005, 0.005, step=0.0005, value=0.001,
                           label="$N_{min}$ (pycnocline: deep value)", show_value=True)
    Nmax_ui = mo.ui.slider(0.005, 0.03, step=0.001, value=0.015,
                           label="$N_{max}$ (pycnocline: peak value)", show_value=True)
    zc_ui = mo.ui.slider(0.1, 0.9, step=0.05, value=0.25,
                         label="pycnocline depth (fraction of $H$)", show_value=True)
    thick_ui = mo.ui.slider(0.02, 0.3, step=0.02, value=0.08,
                            label="pycnocline thickness (fraction of $H$)", show_value=True)
    H_ui = mo.ui.slider(200.0, 5000.0, step=100.0, value=1000.0,
                        label="$H$ (m)", show_value=True)
    f0_ui = mo.ui.slider(3.0e-5, 1.5e-4, step=1.0e-5, value=1.0e-4,
                         label="$f_0$ (1/s)", show_value=True)
    n_modes_ui = mo.ui.slider(2, 6, step=1, value=4, label="number of modes to show",
                              show_value=True)
    mo.vstack([
        mo.hstack([profile_ui, n_modes_ui], justify="start"),
        mo.hstack([N0_ui], justify="start"),
        mo.hstack([Nmin_ui, Nmax_ui], justify="start"),
        mo.hstack([zc_ui, thick_ui], justify="start"),
        mo.hstack([H_ui, f0_ui], justify="start"),
    ])
    return H_ui, N0_ui, Nmax_ui, Nmin_ui, f0_ui, n_modes_ui, profile_ui, thick_ui, zc_ui


@app.cell(hide_code=True)
def _(
    H_ui, N0_ui, Nmax_ui, Nmin_ui, f0_ui, n_modes_ui, np, profile_ui,
    stratification, thick_ui, zc_ui,
):
    # --- solve the eigenvalue problem (cheap, reactive) ---------------------
    _H = H_ui.value
    z = np.linspace(0.0, _H, 300)

    if profile_ui.value == "constant":
        N2 = np.full_like(z, N0_ui.value ** 2)
    else:
        N2 = stratification.pycnocline_N2(
            z, N2_min=Nmin_ui.value ** 2, N2_max=Nmax_ui.value ** 2,
            z_center=zc_ui.value * _H, thickness=thick_ui.value * _H,
        )

    n_modes = n_modes_ui.value
    c, Phi = stratification.vertical_modes(z, N2, _H, n_modes)
    Lr = c / f0_ui.value
    return H_ui, Lr, N2, Phi, c, n_modes, z


@app.cell(hide_code=True)
def _(H_ui, Lr, N2, Phi, c, mo, n_modes, np, plt, z):
    fig1, axs1 = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    axs1[0].plot(np.sqrt(N2) * 1e3, z, lw=2, color="k")
    axs1[0].invert_yaxis()
    axs1[0].set_xlabel("$N(z)$  ($\\times10^{-3}$ s$^{-1}$)"); axs1[0].set_ylabel("$z$ (m, 0=bottom)")
    axs1[0].set_title("stratification profile"); axs1[0].grid(alpha=0.3)

    _colors = plt.cm.viridis(np.linspace(0.1, 0.85, n_modes))
    for _k in range(n_modes):
        axs1[1].plot(Phi[_k], z, lw=2, color=_colors[_k],
                    label=f"mode {_k+1}: $c$={c[_k]:.2f} m/s, $L_R$={Lr[_k]/1e3:.0f} km")
    axs1[1].invert_yaxis()
    axs1[1].axvline(0, color="gray", lw=0.6)
    axs1[1].set_xlabel("$\\Phi_n(z)$"); axs1[1].set_title(f"first {n_modes} baroclinic modes")
    axs1[1].legend(fontsize=7.5, loc="lower left")
    axs1[1].grid(alpha=0.3)
    mo.vstack([
        fig1,
        mo.md(f"**Mode 0 (barotropic):** $c_0=\\infty$, $L_0=\\infty$ — no "
              f"deformation-radius constraint, exactly Ch. 3's "
              f"Taylor-Proudman limit."),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read the mode shapes:** mode 1 changes sign exactly once,
        mode 2 exactly twice, mode 3 exactly three times — the standard
        Sturm-Liouville fact that mode $n$ has $n$ internal zero crossings.
        For a realistic pycnocline (strong, shallow $N^2$ peak), mode 1's
        zero crossing and its largest gradients concentrate near the
        pycnocline depth, NOT at mid-depth — the deformation radius "feels"
        wherever the stratification is actually strong, not the geometric
        center of the water column.
        """
    )
    return


@app.cell(hide_code=True)
def _(H_ui, N0_ui, c, mo, n_modes, np):
    # --- validation against the exact constant-N formula --------------------
    _H = H_ui.value
    _N0 = N0_ui.value
    _c_exact = np.array([_N0 * _H / ((k + 1) * np.pi) for k in range(n_modes)])
    _rel_err = np.abs(c - _c_exact) / _c_exact
    _lines = "\n".join(
        f"- mode {k+1}: computed $c$={c[k]:.4f}, exact $c$={_c_exact[k]:.4f}, "
        f"relative error {_rel_err[k]:.2e}"
        for k in range(n_modes)
    )
    mo.md(
        "### Validation against the exact constant-$N$ formula "
        "$c_n=N_0H/(n\\pi)$\n\n(Switch the profile dropdown to \"Constant "
        "N\" above to see this — with a pycnocline profile these numbers "
        "compare the WRONG two things on purpose, to make clear the exact "
        "formula only applies to the constant-$N$ case.)\n\n" + _lines
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Validate the solver.** Switch to "Constant N" and confirm every
          relative error above is $\lesssim10^{-3}$ — this is the same
          cross-check `tests/test_gfdlib.py` runs automatically, but seeing
          it live, at your own chosen $N_0$, is worth doing once.
        - **Concentrate the stratification.** With the pycnocline profile,
          shrink the thickness slider toward its minimum: watch mode 1's
          deformation radius shrink even though $H$ hasn't changed — a
          thinner, sharper pycnocline traps more of the vertical structure
          in a smaller region, effectively acting like a shallower
          effective $H$.
        - **Match the two-layer model.** Ch. 16's 2-layer $F=f_0^2/(g'H)$
          plays the same role as this chapter's mode-1 $1/L_1^2$. Using a
          pycnocline profile tuned to look like a sharp two-layer interface
          (very thin, very strong), compare $L_1$ here to $1/\sqrt{2F}$ from
          a Ch. 16 parameter choice with matching $g',H$ — they should be
          in the same ballpark (not exact, since two-layer and continuous
          stratification aren't identical models, but the right order of
          magnitude).
        - **Count the crossings.** For $n\_modes=6$, confirm mode 5 crosses
          zero exactly 5 times by eye in the mode-shape panel.

        ### What you should have seen

        A single eigenvalue problem, solved once, generates an entire
        hierarchy of independent vertical structures — mode 0 (no
        constraint, Taylor-Proudman-like), mode 1 (what every earlier
        two-layer calculation in this book called "the" baroclinic mode),
        and higher modes the two-layer model cannot represent at all. Real
        ocean and atmosphere stratification profiles are not constant, and
        their mode 1 deformation radius is set by wherever $N^2(z)$ is
        actually large — usually a thin pycnocline or tropopause region,
        not the water column's geometric middle. This chapter is the
        capstone of the deformation-radius idea threaded through Chs. 6, 8,
        10, and 16: every one of those chapters used SOME value of $L_R$;
        this is where that value actually comes from.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
