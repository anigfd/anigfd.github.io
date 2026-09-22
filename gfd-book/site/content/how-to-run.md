---
title: "How to run these notebooks"
weight: 3
---

Every chapter embeds a live [marimo](https://marimo.io) notebook compiled
to WebAssembly. There are three ways to use them, in increasing order of
commitment.

## 1. In the browser — nothing to install

Just use the chapter page. The embedded notebook runs entirely in your
browser tab via [Pyodide](https://pyodide.org) (CPython compiled to
WebAssembly): there is no server, and nothing you do can affect anyone
else — reload the page to reset everything.

What to expect:

- **First start is the slow part.** The notebook boots a full Python
  runtime plus NumPy and matplotlib on first use — typically 10–40
  seconds and a few tens of MB, depending on your connection. It is
  cached afterward; subsequent chapters start much faster.
- **Sliders react instantly; simulations wait for you.** Cheap panels
  (dispersion relations, analytic solutions, eigenvalue problems)
  recompute live as you drag. Time-stepping simulations sit behind a
  **▶ Run** button so a half-dragged slider doesn't launch a hundred
  runs — set the controls, then press Run.
- **The browser is 3–10× slower than native Python.** Every simulation
  chapter exposes a *resolution* control and a run-length/steps control:
  if a run feels slow, drop the resolution first — for the physics these
  notebooks demonstrate, $64^2$ almost always shows the same phenomenon
  as $128^2$.
- Phones work for reading, but the simulations want a laptop's CPU and
  screen width.

## 2. Locally — for the exercises

Each chapter's *computational* exercise asks you to modify the notebook,
and the natural place to do that is a local checkout:

```bash
git clone https://github.com/aneeshcs/gfd_textbook_aneesh
cd gfd_textbook_aneesh/gfd-book
pip install -r requirements.txt        # marimo, numpy, matplotlib, pytest
marimo edit notebooks/ch06_geostrophic-adjustment.py
```

`marimo edit` opens the notebook in your browser as an editable, reactive
document — change any cell and everything downstream recomputes. The
notebooks import the shared library `gfdlib` straight from the repo (the
import cell at the top of every notebook handles the local-vs-browser
difference automatically), so edits to `gfdlib` take effect immediately
too.

Three things worth knowing before you start editing:

- **Notebooks are plain Python files.** A marimo notebook is a `.py` file
  you can read, diff, and version-control like any other code — there is
  no JSON, and `git diff` on your modified notebook shows exactly what
  you changed.
- **Cells form a dataflow graph, not a top-to-bottom script.** marimo
  re-runs a cell when a variable it *uses* changes, regardless of cell
  order. Names prefixed with `_` are private to their cell; everything
  else is shared.
- **The numerics live in [`gfdlib`](https://github.com/aneeshcs/gfd_textbook_aneesh/tree/main/gfd-book/gfdlib), not in the notebooks.** If you want
  to change *how* something is computed (a solver, a time-stepper), edit
  `gfdlib` and run `make test` — the test suite pins every primitive to
  an analytic solution, a conservation law, or a convergence order, and
  it will tell you precisely which physical guarantee your change broke.

## 3. Rebuilding the site — for contributors

To reproduce the full pipeline (export every notebook to WASM and preview
the site):

```bash
make test        # gfdlib correctness — must pass
make notebooks   # export all notebooks -> site/static/nb/
make serve       # hugo server -D, preview at localhost:1313
```

See the contributor docs in the repo — `docs/architecture.md` for how the
marimo → WASM → Hugo pipeline works and the performance budget behind it,
and `docs/authoring.md` for the add-a-chapter checklist. Symbols and sign
conventions for all chapters are fixed on the
[notation]({{< relref "notation.md" >}}) page.
