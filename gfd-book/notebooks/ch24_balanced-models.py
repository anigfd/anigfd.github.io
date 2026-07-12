import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter 24 — Balanced Models & the Slow Manifold

        **Physical question.** This book has spent 23 chapters solving
        *reduced* equations — geostrophy, QG, two-layer, planetary-scale
        balances — while the fluid itself obeys none of them exactly. Why
        did that ever work? Because rotating, stratified flow segregates
        into **fast** motions (inertia-gravity waves, Chs. 6 and 12) and
        **slow** ones (everything PV controls), and a state prepared with
        *zero fast content* stays — very nearly — fast-free: it lives on
        the **slow manifold**, and the balanced models are simply the
        dynamics restricted to it. After this final chapter you should be
        able to say what the slow manifold is, demonstrate with two runs
        of the same model why initialization onto it matters, and place
        every reduced model in this book on the same map.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations

        ### The split, one more time — now as geometry

        Ch. 6 diagonalized linear rotating shallow water at each
        wavenumber into three modes: one **slow** ($\omega=0$, the
        geostrophic mode) and two **fast**
        ($\omega=\pm\sqrt{1+k^2+l^2}\ge1$: inertia-gravity waves). Read
        that as geometry rather than as a wave census: state space splits
        into a *slow subspace* (all fast amplitudes zero — exactly the
        geostrophically balanced states $u=-\eta_y$, $v=\eta_x$) and its
        fast complement. In the **linear** system the split is exact and
        invariant: a state that starts on the slow subspace stays on it
        *forever* — it is an exact steady state on the $f$-plane (the
        book's test suite literally asserts this) — while any fast
        content rings at $\omega\ge1$, forever, without touching the slow
        part. And crucially, the frequency **gap** ($0$ versus $\ge f$)
        means the two families cannot resonate at small amplitude.

        ### Balanced initialization: projecting onto the manifold

        The practical consequence is the oldest problem in numerical
        weather prediction. Observations give you $\eta$ (pressure) far
        better than winds; initialize a model with observed $\eta$ and
        naive winds and the error is not merely a wrong forecast — it is
        a **fast one**: the projection of your error onto the fast modes
        rings through the model at gravity-wave frequency, swamping the
        meteorology (Richardson's famous 1922 hand-computed forecast
        failed exactly this way, predicting a 145 hPa pressure change in
        6 hours). The cure is to *project the initial state onto the slow
        manifold*: set the winds to their balanced values. In the linear
        model that projection is exact; the notebook below performs it
        with one slider.

        ### From subspace to manifold — and how fuzzy it really is

        Nonlinearity bends the slow *subspace* into a slow *manifold*:
        the balanced winds acquire $O(Ro)$ corrections (gradient-wind,
        Ch. 5; the QG ageostrophic circulation, Ch. 8), computable order
        by order in $Ro$ — and each of this book's reduced models is
        precisely the dynamics on this manifold, truncated at some order:

        | Model | Where it lives |
        |---|---|
        | geostrophy (Ch. 5) | the manifold at $O(1)$ — diagnostic only |
        | QG (Chs. 8, 10, 16, 18) | dynamics on the manifold at $O(Ro)$ |
        | planetary/Sverdrup balances (Chs. 20, 21) | the manifold at large scale, friction added |
        | full equations (Chs. 6, 12) | the whole space, fast modes included |

        One honest caveat closes the book. Lorenz asked whether the
        exactly-invariant slow manifold survives nonlinearity, and the
        answer is: *not quite*. Balanced motion generates a whisper of
        gravity waves — **spontaneous imbalance** — but the leak is
        *exponentially* small in Rossby number, $\sim e^{-c/Ro}$ (Vanneste
        & Yavneh 2004), which is why balance is such an unreasonably good
        approximation at $Ro\sim0.1$ even though no exact slow manifold
        exists. "Balance" is not a fact about the fluid; it is an
        asymptotic property of where you chose to start it — maintained,
        at small $Ro$, to all orders.
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
    from gfdlib import spectral, timestep, shallowwater
    return np, plt, shallowwater, spectral, timestep


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Numerical scheme

        Ch. 6's exact machinery, unchanged (`gfdlib.spectral.Grid`,
        `gfdlib.shallowwater.rhs`, RK4 at the same CFL) on an $f$-plane,
        plus one new one-line diagnostic: `gfdlib.shallowwater.divergence`
        ($\delta=u_x+v_y$, spectral). Divergence is the sharpest possible
        slow/fast meter here — the slow mode has $\delta=0$ *identically*,
        while inertia-gravity waves are made of it — so RMS $\delta$
        measures the distance from the slow manifold directly. Three runs
        share the identical height field $\eta_0$ and differ only in
        initial winds: $(u,v)=\alpha\,(-\eta_{0y},\ \eta_{0x})$ with
        $\alpha=1$ (balanced), your slider's $\alpha$, and $\alpha=0$
        (Richardson's mistake).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    alpha_ui = mo.ui.slider(0.0, 1.0, step=0.05, value=0.7,
                            label="balance fraction $\\alpha$", show_value=True)
    sigma_ui = mo.ui.slider(0.5, 4.0, step=0.25, value=2.0,
                            label="bump width $\\sigma/L_R$", show_value=True)
    T_ui = mo.ui.slider(5.0, 25.0, step=2.5, value=15.0,
                        label="run time $T$ ($f_0^{-1}$)", show_value=True)
    n_ui = mo.ui.dropdown(options={"64": 64, "96": 96}, value="64",
                          label="resolution $n$")
    mo.vstack([
        mo.hstack([alpha_ui, sigma_ui], justify="start"),
        mo.hstack([T_ui, n_ui], justify="start"),
        mo.md(r"> $\alpha=1$: winds exactly geostrophic (on the slow "
              r"subspace). $\alpha=0$: same height field, no winds at all "
              r"— Ch. 6's adjustment problem, now reread as an "
              r"initialization *error*."),
    ])
    return T_ui, alpha_ui, n_ui, sigma_ui


@app.cell(hide_code=True)
def _(mo):
    get_sim, set_sim = mo.state(None)   # results survive the run button resetting
    return get_sim, set_sim


@app.cell(hide_code=True)
def _(mo):
    run_btn = mo.ui.run_button(label="▶ Run all three initializations")
    run_btn
    return (run_btn,)


@app.cell
def _(
    T_ui, alpha_ui, mo, n_ui, np, run_btn, set_sim, shallowwater, sigma_ui,
    spectral, timestep,
):
    mo.stop(not run_btn.value, mo.callout(
        mo.md("Press **▶ Run all three initializations** above to start."),
        kind="info",
    ))

    _n, _alpha, _sigma, _T = n_ui.value, alpha_ui.value, sigma_ui.value, T_ui.value
    _L = 20.0
    _grid = spectral.Grid(_n, L=_L)
    _f = shallowwater.coriolis(_grid, beta=0.0)

    _eta0 = 0.5 * np.exp(
        -((_grid.x - _L / 2) ** 2 + (_grid.y - _L / 2) ** 2) / (2 * _sigma ** 2)
    )
    _eta_hat = _grid.fft(_eta0)
    _ug = _grid.ifft(-_grid.ddy(_eta_hat))
    _vg = _grid.ifft(_grid.ddx(_eta_hat))

    _dt = 0.4 * _grid.dx
    _nsteps = max(50, int(_T / _dt))
    _nsub = max(1, _nsteps // 300)
    _ic = _grid.n // 2

    def _run(_a):
        _state = np.stack([_eta0, _a * _ug, _a * _vg])
        _q0 = shallowwater.potential_vorticity(_state, _grid)
        _divs, _etas, _ts = [], [], []
        _t = 0.0
        for _s in range(_nsteps + 1):
            if _s % _nsub == 0:
                _delta = shallowwater.divergence(_state, _grid)
                _divs.append(float(np.sqrt((_delta ** 2).mean())))
                _etas.append(float(_state[0][_ic, _ic]))
                _ts.append(_t)
            if _s < _nsteps:
                _state = timestep.rk4(lambda t, s: shallowwater.rhs(s, _grid, _f),
                                      _state, _dt)
                _t += _dt
        _dq = float(np.abs(shallowwater.potential_vorticity(_state, _grid) - _q0).max())
        return np.array(_ts), np.array(_divs), np.array(_etas), _dq

    _t1, _d1, _e1, _dq1 = _run(1.0)
    _ta, _da, _ea, _dqa = _run(_alpha)
    _t0, _d0, _e0, _dq0 = _run(0.0)

    set_sim({
        "t": _t1,
        "div1": _d1, "diva": _da, "div0": _d0,
        "eta1": _e1, "etaa": _ea, "eta0s": _e0,
        "dq": (_dq1, _dqa, _dq0), "alpha": _alpha,
    })
    mo.callout(mo.md(
        f"**Done.** 3 runs × {_nsteps} steps. Max pointwise PV drift over "
        f"the entire run, per run: {_dq1:.1e} / {_dqa:.1e} / {_dq0:.1e} — "
        f"machine precision in all three."
    ), kind="success")
    return


@app.cell
def _(get_sim, mo):
    sim = get_sim()
    mo.stop(sim is None, mo.callout(
        mo.md("Results will appear here once the runs have completed."), kind="warn"
    ))
    return (sim,)


@app.cell(hide_code=True)
def _(np, plt, sim):
    fig1, axs1 = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    axs1[0].plot(sim["t"], sim["div0"], lw=1.5, color="crimson",
                 label="$\\alpha=0$ (no winds)")
    axs1[0].plot(sim["t"], sim["diva"], lw=1.5, color="#d97706",
                 label=f"$\\alpha={sim['alpha']:.2f}$")
    axs1[0].plot(sim["t"], sim["div1"] + 1e-18, lw=1.5, color="#2563eb",
                 label="$\\alpha=1$ (balanced)")
    axs1[0].set_yscale("log")
    axs1[0].set_xlabel("$t\\ (f_0^{-1})$"); axs1[0].set_ylabel("RMS divergence")
    axs1[0].set_title("distance from the slow manifold")
    axs1[0].legend(fontsize=8); axs1[0].grid(alpha=0.3, which="both")

    axs1[1].plot(sim["t"], sim["eta0s"], lw=1.2, color="crimson")
    axs1[1].plot(sim["t"], sim["etaa"], lw=1.2, color="#d97706")
    axs1[1].plot(sim["t"], sim["eta1"], lw=1.8, color="#2563eb")
    axs1[1].set_xlabel("$t\\ (f_0^{-1})$"); axs1[1].set_ylabel("$\\eta$ at bump center")
    axs1[1].set_title("the 'forecast' at one station")
    axs1[1].grid(alpha=0.3)
    fig1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        **How to read it:** the balanced run's divergence sits at machine
        zero — it is *on* the slow subspace and stays there, an exact
        steady state. Every other $\alpha$ rings at inertia-gravity
        frequencies, and the ringing amplitude scales exactly as
        $(1-\alpha)$ — the fast-mode content of the initialization error,
        nothing else (the dynamics is linear, so the fast and slow parts
        evolve in strict parallel). The right panel is Richardson's story
        at one grid point: the same "observations" ($\eta_0$), and either
        a quiet, meteorologically meaningful record or one buried under
        gravity-wave noise, depending purely on whether the winds were
        initialized on the manifold. Meanwhile the success message above
        reports the PV drift in every run: machine precision regardless
        of $\alpha$ — the slow content was never in danger; *balance is
        entirely about not exciting the fast modes*.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this

        - **Verify the $(1-\alpha)$ law.** Read the plateau RMS divergence
          off the log plot at $\alpha=0.5$ and $\alpha=0$: the ratio
          should be exactly 2. Why exactly? (Linearity: the fast-mode
          amplitude is a linear projection of the initialization error,
          which is proportional to $1-\alpha$.)
        - **Change what "small" means.** Shrink $\sigma/L_R$ toward 0.5:
          Ch. 6 taught that small-scale height bumps are mostly fast
          content — so even the *balanced* winds now protect a state
          that holds little slow signal. Compare the $\alpha=0$ ringing
          amplitude across $\sigma$: initialization matters most exactly
          where there is slow signal worth protecting.
        - **Where the analogy stops.** This linear model's slow subspace
          is exactly invariant, so $\alpha=1$ gives *zero* ringing
          forever. In the real (nonlinear) primitive equations, an
          identically-balanced start still leaks gravity waves —
          spontaneously, at $\sim e^{-c/Ro}$. Nothing in this notebook
          can show that (it is exponentially beyond a linear model);
          knowing which of your conclusions survive nonlinearity is the
          course's parting skill.

        ### What you should have seen — and where the book ends

        One height field, three wind initializations: identical slow
        (PV) content to machine precision, and utterly different fast
        content — from silence at $\alpha=1$ to full-volume ringing at
        $\alpha=0$, in exact proportion to the projection off the slow
        manifold. That is the entire logic of this book run in reverse.
        Chs. 5–21 built models by *assuming* the fast modes away —
        geostrophy at $O(1)$, QG at $O(Ro)$, Sverdrup at basin scale —
        and they worked because rotating stratified flow, initialized
        with even minimal care, genuinely keeps its fast and slow lives
        separate: the frequency gap protects the split, PV anchors the
        slow side (Ch. 7), waves carry the fast side (Chs. 6, 12), and
        the two meet only through the slow, patient couplings of
        Chs. 13 and 23. The reduced models are not approximations that
        happen to work; they are the dynamics of the manifold the
        atmosphere and ocean actually live near — which is why a
        24-chapter book about two fluids could be, most of the time, a
        book about one idea.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
