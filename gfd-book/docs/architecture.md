# Architecture — how a notebook becomes a browser-runnable chapter

This page explains the build pipeline end to end: what happens between
editing `notebooks/chNN_slug.py` and a reader turning a slider on the live
site, and the constraints (all of them performance constraints) that shape
how chapters are written.

## The pipeline at a glance

```
gfdlib/*.py ──(pip wheel)──▶ notebooks/public/gfdlib-0.1.0-py3-none-any.whl
                                        │
notebooks/chNN_slug.py ──(marimo export html-wasm)──▶ site/static/nb/chNN_slug/
                                        │                       │
                                        │       Hugo shortcode {{< marimo src="/nb/chNN_slug/" >}}
                                        ▼                       ▼
                              reader's browser boots Pyodide, micropip-installs
                              the bundled wheel, and runs the notebook locally
```

There is no server-side compute anywhere. The exported bundle is static
files; all numerics run in the reader's browser via Pyodide (CPython
compiled to WebAssembly).

## Step by step

### 1. `gfdlib` is packaged as a wheel

`make wheel` (a dependency of `make notebooks`) runs
`pip wheel --no-deps -w notebooks/public .` whenever any `gfdlib/*.py` or
`pyproject.toml` changes. The wheel lands in `notebooks/public/` because
marimo copies a notebook's `public/` folder into every WASM export — that is
how the shared library travels with each chapter bundle.

**Only put things in `notebooks/public/` that every chapter should carry.**
Because all notebooks are siblings in one flat `notebooks/` directory, they
all share this one `public/` folder — marimo has no way to scope it to a
single chapter, so anything placed there (the 40 KB wheel; fine) gets
duplicated into all 24 exports. A chapter-specific precomputed dataset (see
the performance-budget note below) belongs in `notebooks/data/<name>.npz`
instead: `make notebooks` copies `notebooks/data/chNN_slug.npz` into only
that chapter's own `site/static/nb/chNN_slug/public/` after export, so it
isn't paid for by the other 23 chapters. (Discovered when a single 5 MB file
in the shared folder was quietly adding ~125 MB to the built site.)

### 2. Every notebook has a dual-mode import cell

Each notebook (see `notebooks/_template.py`) begins with:

```python
if sys.platform == "emscripten":          # running under Pyodide in the browser
    import micropip
    await micropip.install(
        str(mo.notebook_location() / "public" / "gfdlib-0.1.0-py3-none-any.whl")
    )
else:                                      # running locally (marimo edit / pytest)
    sys.path.insert(0, str(mo.notebook_dir().parent))   # repo root

from gfdlib import spectral, timestep, ...
```

Locally the notebook imports `gfdlib` straight from the repo; in the browser
it installs the bundled wheel. Same import line either way, so the code the
student reads never branches.

### 3. Export to WASM

`make notebooks` runs, for every `notebooks/ch*.py`:

```bash
marimo export html-wasm notebooks/chNN_slug.py -o site/static/nb/chNN_slug --mode run
```

`--mode run` gives the reader a read-and-interact app (sliders work, code is
visible but not editable). Use `--mode edit` for chapters where the reader
should modify code live. `site/static/nb/` is git-ignored — bundles are
build products, rebuilt from source.

### 4. Hugo embeds the bundle

Each chapter page in `site/content/partN/chNN_slug.md` embeds its notebook
with the `marimo` shortcode:

```
{{< marimo src="/nb/chNN_slug/" >}}
```

`make serve` previews the whole site (`hugo server -D` from `site/`).
Deploying means building the Hugo site and publishing `site/public/` — never
touch the live site without explicit approval.

## The performance budget (the real constraint)

Pyodide is single-threaded and ~3–10× slower than native CPython. Every
notebook design decision traces back to this:

- **Grids stay modest**: 64²–256² for spectral solvers. The included solver
  does 128² × 400 steps in under a second native, so a few seconds in the
  browser.
- **Precompute once, outside the time loop**: build one
  `gfdlib.spectral.Grid` at the top of the notebook and reuse it. Wavenumber
  arrays, dealias masks, and integrating-factor exponentials must never be
  rebuilt per step.
- **Expose speed/detail trade-offs to the reader**: a resolution slider and
  a steps-per-frame slider are standard controls.
- **Precompute genuinely heavy runs offline** and ship the fields as data
  the widget scrubs through (e.g. a full baroclinic life cycle), rather than
  integrating live — put the file in `notebooks/data/<name>.npz`, not
  `notebooks/public/` (see above).
- **Lazy-load**: Pyodide does not boot until the reader clicks Run.

## Testing and CI

- `make test` runs `pytest -q tests/` — numeric correctness of every
  `gfdlib` primitive (analytic solutions, conservation laws, convergence
  orders). It must pass before every commit.
- CI (GitHub Actions, on every push/PR) runs the test suite **and exports
  every notebook to WASM** as a smoke check, so an export-breaking change
  is caught before it ships.

## Reproducibility

`requirements.txt` pins exact versions (marimo, numpy, matplotlib, pytest)
rather than floors — the source of truth for authoring-side versions. This
also pins the reader-side Pyodide runtime: marimo's own version fully
determines which Pyodide build `marimo export html-wasm` bundles (see
`PYODIDE_VERSION` in `marimo/_pyodide/pyodide_constraints.py`), so pinning
marimo is what makes every reader's in-browser environment reproducible,
not just the authoring environment. The repo is tagged per "edition"; the
tag marks a commit where the pinned versions, all 24 chapters, and a clean
CI run (tests + WASM export of every notebook) coincide.

marimo notebooks are plain `.py` files, so chapter diffs review like
ordinary code.
