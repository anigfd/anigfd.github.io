import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 22 — Hamiltonian & Variational GFD

        **Physical question.** Integrate the same vortex system with a
        humble 2nd-order method and with classical 4th-order RK4, and watch
        the energy: RK4's error is far smaller — and grows *forever*, while
        the cruder method's error, if the method is chosen with the right
        geometry, stays bounded until the end of time. Why? Because the
        dynamics is **Hamiltonian**, and one integrator respects that
        structure while the other does not. After this chapter you should
        be able to say what "the structure" is, why GFD's conservation laws
        (energy, momentum, and above all PV) are *symmetries* rather than
        accidents, and when a structure-preserving method is worth its
        price.

        *This chapter opens the book's optional Part VII — the Salmon
        capstone. Nothing later in your career depends on it; everything
        earlier in this book quietly did.*
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### Least action, Noether, and the fluid's secret symmetry

        Ideal fluid mechanics can be derived from a variational principle,
        exactly like particle mechanics: parcels follow trajectories that
        make the action $\int(KE-PE)\,dt$ stationary. The payoff of saying
        it this way is **Noether's theorem**: every continuous symmetry of
        the action yields a conservation law —

        | Symmetry | Conserved quantity |
        |---|---|
        | time translation | energy |
        | space translation | momentum |
        | rotation | angular momentum |
        | **particle relabeling** | **potential vorticity / Kelvin circulation** |

        The last row is the fluid-specific entry, and it answers a question
        this book has leaned on since Ch. 7 without ever answering: *why*
        is PV materially conserved? Because a fluid, unlike a particle
        system, carries a continuum of identical parcels — and swapping the
        *labels* of parcels (without moving any actual fluid) changes
        nothing physical. That relabeling freedom is a continuous symmetry
        of the action, and it is an *infinite-dimensional* one: Noether
        then delivers not one conserved number but a conserved quantity
        **per parcel** — precisely a materially conserved scalar,
        $Dq/Dt=0$. Energy conservation is ordinary; PV conservation is the
        signature of being a fluid. (Salmon 1988 is the classic exposition;
        every conservation law used in Chs. 6–21 sits in this table.)

        ### Point vortices: a Hamiltonian system you can see

        The cleanest laboratory for these ideas is one this book already
        owns: $N$ point vortices (Ch. 7's blobs, taken to zero size) are an
        **exact** finite-dimensional reduction of 2D Euler flow, and
        Kirchhoff (1876) showed they form a Hamiltonian system with a
        twist (symbols as in [NOTATION](../notation)):

        $$\Gamma_i\frac{dx_i}{dt}=\frac{\partial H}{\partial y_i},\qquad
          \Gamma_i\frac{dy_i}{dt}=-\frac{\partial H}{\partial x_i},\qquad
          H=-\frac{1}{4\pi}\sum_{i<j}\Gamma_i\Gamma_j\ln r_{ij}^2.$$

        The twist: $x_i$ and $y_i$ are *each other's* conjugate variables —
        there is no separate momentum. **Phase space is physical space.** A
        snapshot of vortex positions is a complete dynamical state (no
        velocities needed — compare Ch. 7's invertibility, of which this is
        the point-vortex shadow), and $\Gamma$-weighted area in the plane
        is the conserved symplectic measure.

        Noether's table now reads off directly from what $H$ ignores. $H$
        depends only on separations, so:

        $$\text{translation}\Rightarrow Q=\sum_i\Gamma_ix_i,\quad
          P=\sum_i\Gamma_iy_i;\qquad
          \text{rotation}\Rightarrow L=\sum_i\Gamma_i(x_i^2+y_i^2);$$

        plus $H$ itself. Four conserved quantities — enough to make the
        3-vortex problem integrable (regular, quasi-periodic forever),
        while $N\ge4$ is generically **chaotic** (Ch. 14's sensitivity, in
        a system with four exact invariants!).

        ### Symplectic integration: solving a nearby problem exactly

        A numerical one-step method is a map of phase space. Hamiltonian
        flow is a very special map — it preserves the symplectic
        structure (in 2D: signed, $\Gamma$-weighted area). Most
        integrators, including RK4, do not; the **implicit midpoint rule**

        $$s_{n+1}=s_n+\Delta t\,f\!\left(\frac{s_n+s_{n+1}}{2}\right)$$

        does, for *any* Hamiltonian (even this non-separable one, which
        rules out the more famous leapfrog). The consequence is explained
        by **backward error analysis**: a symplectic method's numerical
        trajectory is the *exact* trajectory of a slightly-perturbed
        Hamiltonian $\tilde H=H+O(\Delta t^p)$. Energy error can therefore
        only oscillate within an $O(\Delta t^p)$ band — forever. A
        non-symplectic method solves no nearby Hamiltonian problem, and
        its energy error performs a slow random walk: small per step,
        **secular** in sum. Accuracy and structure are different axes:
        RK4 is more accurate per step; midpoint is more faithful per
        epoch. (Bonus fact, verified in `tests/`: as a Gauss method,
        implicit midpoint also conserves all *quadratic* invariants — $L$
        here — exactly, while RK4 lets $L$ drift. $Q$ and $P$, being
        linear, are conserved by every Runge-Kutta method.)
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
    from gfdlib import vortex, timestep
    return np, plt, timestep, vortex


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        `gfdlib.vortex` supplies the mutual Biot-Savart RHS, the
        Hamiltonian and its three companion invariants, and two
        integrators: `step_midpoint` (implicit midpoint, symplectic,
        2nd order, solved by fixed-point iteration) and `step_heun`
        (explicit RK2 — the fair same-order, same-cost, non-symplectic
        control). RK4 comes from `gfdlib.timestep` as everywhere else in
        the book. Known-solution checks (co-rotating pair period
        $T=2\pi^2d^2/\Gamma$, dipole speed $\Gamma/2\pi d$) live in the
        test suite.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    config_ui = mo.ui.dropdown(
        options={
            "4 vortices — chaotic": "chaos4",
            "3 vortices — integrable": "int3",
            "Co-rotating pair — exact solution known": "pair",
        },
        value="4 vortices — chaotic",
        label="configuration",
    )
    dt_ui = mo.ui.slider(0.01, 0.2, step=0.01, value=0.05,
                         label="$\\Delta t$", show_value=True)
    T_ui = mo.ui.slider(100.0, 1000.0, step=100.0, value=400.0,
                        label="run time $T$", show_value=True)
    mo.vstack([
        mo.hstack([config_ui, dt_ui, T_ui], justify="start"),
        mo.md(r"> All three integrators (Heun, implicit midpoint, RK4) run "
              r"on the identical configuration and $\Delta t$. The energy "
              r"plot is the whole story; the trajectory plot is there to "
              r"show the dynamics being integrated."),
    ])
    return T_ui, config_ui, dt_ui


@app.cell(hide_code=True)
def _(mo):
    get_sim, set_sim = mo.state(None)   # results survive the run button resetting
    return get_sim, set_sim


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Integrate (3 methods)")
    run_btn
    return (run_btn,)


@app.cell
def _(T_ui, config_ui, dt_ui, mo, np, run_btn, set_sim, timestep, vortex):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Integrate (3 methods)** above to start."), kind="info"
    ))

    _cfg = config_ui.value
    if _cfg == "chaos4":
        _Gam = np.array([1.0, 1.0, 1.0, -0.5])
        _s0 = np.array([[0.0, 1.0, 0.3, -0.8], [0.0, 0.0, 0.9, 0.4]])
    elif _cfg == "int3":
        _Gam = np.array([1.0, 1.0, 1.0])
        _s0 = np.array([[0.0, 1.1, 0.4], [0.0, 0.0, 0.8]])
    else:
        _Gam = np.array([1.0, 1.0])
        _s0 = np.array([[-0.5, 0.5], [0.0, 0.0]])

    _dt, _T = dt_ui.value, T_ui.value
    _nsteps = int(_T / _dt)
    _nsub = max(1, _nsteps // 400)

    _H0 = vortex.hamiltonian(_s0, _Gam)
    _, _, _, _L0 = vortex.invariants(_s0, _Gam)

    _sh, _sm, _sr = _s0.copy(), _s0.copy(), _s0.copy()
    _t_hist, _eH, _eM, _eR = [], [], [], []
    _LM, _LR = [], []
    _traj = []
    _t = 0.0
    for _s in range(_nsteps + 1):
        if _s % _nsub == 0:
            _t_hist.append(_t)
            _eH.append(vortex.hamiltonian(_sh, _Gam) - _H0)
            _eM.append(vortex.hamiltonian(_sm, _Gam) - _H0)
            _eR.append(vortex.hamiltonian(_sr, _Gam) - _H0)
            _LM.append(vortex.invariants(_sm, _Gam)[3] - _L0)
            _LR.append(vortex.invariants(_sr, _Gam)[3] - _L0)
            _traj.append(_sm.copy())
        if _s < _nsteps:
            _sh = vortex.step_heun(_sh, _Gam, _dt)
            _sm = vortex.step_midpoint(_sm, _Gam, _dt)
            _sr = timestep.rk4(lambda t, st: vortex.vortex_rhs(st, _Gam), _sr, _dt)
            _t += _dt

    set_sim({
        "t": np.array(_t_hist),
        "eH": np.array(_eH), "eM": np.array(_eM), "eR": np.array(_eR),
        "LM": np.array(_LM), "LR": np.array(_LR),
        "traj": np.array(_traj), "Gamma": _Gam, "H0": _H0,
    })
    mo.callout(mo.md(f"**Done.** {_nsteps} steps × 3 methods."), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the integration has run."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(np, plt, sim):
    fig1, axs1 = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)

    _traj = sim["traj"]                     # (frames, 2, N)
    _N = _traj.shape[2]
    _colors = plt.cm.plasma(np.linspace(0.15, 0.85, _N))
    for _i in range(_N):
        axs1[0].plot(_traj[:, 0, _i], _traj[:, 1, _i], lw=0.6, color=_colors[_i],
                     alpha=0.8)
        axs1[0].plot(_traj[0, 0, _i], _traj[0, 1, _i], "o", color=_colors[_i],
                     ms=7, markeredgecolor="k")
    axs1[0].set_xlabel("$x$"); axs1[0].set_ylabel("$y$")
    axs1[0].set_title("trajectories (midpoint run; dots = start)")
    axs1[0].set_aspect("equal"); axs1[0].grid(alpha=0.3)

    for _e, _lbl, _c in [(sim["eH"], "Heun (2nd, not symplectic)", "crimson"),
                         (sim["eR"], "RK4 (4th, not symplectic)", "#d97706"),
                         (sim["eM"], "midpoint (2nd, symplectic)", "#2563eb")]:
        axs1[1].plot(sim["t"], np.abs(_e) / abs(sim["H0"]), lw=1.2, color=_c, label=_lbl)
    axs1[1].set_yscale("log")
    axs1[1].set_xlabel("$t$"); axs1[1].set_ylabel("$|H(t)-H_0|/|H_0|$")
    axs1[1].set_title("energy error (log scale)")
    axs1[1].legend(fontsize=8); axs1[1].grid(alpha=0.3, which="both")
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read the energy panel:** Heun climbs steadily — a secular
        drift, the signature of *not* solving any nearby Hamiltonian
        problem. Midpoint, the same order and roughly the same cost, sits
        in a flat, bounded band: it is solving a nearby Hamiltonian
        problem *exactly*, so its energy can only oscillate. RK4 starts
        far below both (4th-order accuracy is real) but look at its
        *slope* on the log plot — it grows without bound, and given enough
        time must cross midpoint's flat band. Accuracy buys you a lower
        starting point; structure buys you the ceiling.
        """
    )
    return


@app.cell(hide_code=True)
def _(np, plt, sim):
    fig2, ax2 = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    ax2.plot(sim["t"], np.abs(sim["LR"]) + 1e-18, lw=1.2, color="#d97706",
             label="RK4")
    ax2.plot(sim["t"], np.abs(sim["LM"]) + 1e-18, lw=1.2, color="#2563eb",
             label="midpoint")
    ax2.set_yscale("log")
    ax2.set_xlabel("$t$"); ax2.set_ylabel("$|L(t)-L_0|$")
    ax2.set_title("angular impulse error — quadratic invariant")
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3, which="both")
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **The quadratic-invariant panel** is a finer diagnostic than
        energy: $L=\sum\Gamma_i(x_i^2+y_i^2)$ is *quadratic* in the state,
        and Gauss methods (implicit midpoint is the simplest one) conserve
        every quadratic invariant **exactly** — midpoint's curve is pure
        fixed-point-iteration residual, near machine precision. RK4
        genuinely drifts. (The linear invariants $Q,P$ are conserved to
        machine precision by *all* Runge-Kutta methods and are therefore
        not worth plotting — conservation-law bookkeeping is a hierarchy,
        not a single yes/no.)

        ### Try this

        - **Make RK4 lose.** Increase $\Delta t$ toward 0.2 and $T$ toward
          1000: RK4's secular error ($\propto\Delta t^4\,t$) climbs while
          midpoint's band ($\propto\Delta t^2$, flat) merely widens. Find
          a setting where RK4's energy error actually crosses above
          midpoint's — the "worse method wins" moment, live.
        - **Integrable vs. chaotic.** Switch to the 3-vortex
          configuration: the trajectories become visibly regular
          (quasi-periodic ribbons instead of tangles). Four invariants in
          a 6-dimensional phase space is enough to close the orbits;
          adding one vortex breaks integrability — Ch. 14's route to
          chaos, without ever leaving a Hamiltonian system.
        - **Watch the pair stay honest.** The co-rotating pair has an
          exact solution (period $2\pi^2d^2/\Gamma$); with a large
          $\Delta t$, watch which integrator's trajectory ring stays a
          circle longest.
        - **Break the iteration.** Edit the run cell to call
          `step_midpoint(..., n_iter=1)`: with the implicit equation
          barely solved, is the method still effectively symplectic? (The
          $L$ panel will answer before the energy panel does.)

        ### What you should have seen

        A Hamiltonian system's conservation laws are geometry, not
        bookkeeping — and a numerical method either carries that geometry
        or it doesn't, *independently of its order of accuracy*. The same
        lesson scales up: climate models run for centuries of model time,
        far beyond any hope of trajectory accuracy, and what keeps their
        statistics physical is respecting conservation structure —
        mimetic discretizations, conservative advection schemes — not
        raw order. The deepest structure of all is the fluid's
        particle-relabeling symmetry, whose Noether charge is the PV this
        entire book has been built on; Ch. 23 turns that symmetry loose
        on waves.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
