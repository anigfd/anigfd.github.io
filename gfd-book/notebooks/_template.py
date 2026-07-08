import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Chapter NN — <Title>

        **Physical question.** <One or two sentences: what should the reader be
        able to predict after this chapter?>
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Governing equations
        Use the SAME symbols as `NOTATION.md`. For example the barotropic
        vorticity equation:

        $$\frac{\partial \zeta}{\partial t} + J(\psi,\zeta) = \nu\nabla^2\zeta,
          \qquad \nabla^2\psi = \zeta.$$
        """
    )
    return


@app.cell
async def _(mo):
    # --- numerical scheme: import vetted primitives, don't re-implement ---
    import sys

    if sys.platform == "emscripten":
        # browser (Pyodide/WASM): install the gfdlib wheel that `make notebooks`
        # ships in the exported bundle's public/ folder
        import micropip
        await micropip.install(
            str(mo.notebook_location() / "public" / "gfdlib-0.1.0-py3-none-any.whl")
        )
    else:
        sys.path.insert(0, str(mo.notebook_dir().parent))  # repo root

    import numpy as np
    from gfdlib import spectral, timestep, diagnostics, plotting
    grid = spectral.Grid(128)          # build ONCE, reuse every step
    return np, spectral, timestep, diagnostics, plotting, grid


@app.cell(hide_code=True)
def _(mo):
    # --- interactive controls: the 1-3 dimensionless numbers of this problem ---
    res   = mo.ui.dropdown({"64": 64, "128": 128, "256": 256}, value="128", label="resolution")
    param = mo.ui.slider(0.0, 1.0, value=0.2, step=0.01, label="control parameter")
    mo.hstack([res, param])
    return res, param


@app.cell
def _(grid, param, plotting):
    # --- diagnostics: fields + a budget/spectrum ---
    fig, ax = plotting.new_fig()
    # ... compute a field f and show it:
    # plotting.field(ax, f, grid=grid, signed=True, title="vorticity")
    fig
    return (fig,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Try this
        - <exploratory prompt 1>
        - <exploratory prompt 2>

        ### What you should have seen
        <the expected qualitative result, so the reader can check intuition>
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
