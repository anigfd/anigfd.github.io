# Architecture — how a notebook becomes a browser-runnable chapter

This page explains the build pipeline end to end: what happens between
editing `notebooks/chNN_slug.py` and a reader turning a slider on the live
site, and the constraints (all of them performance constraints) that shape
how chapters are written.

## The pipeline at a glance

```
gfdlib/*.py ──(pip wheel)──▶ notebooks/public/gfdlib-0.1.0-py3-none-any.whl
                                        │
notebooks/chNN_slug.py ──(marimo export html-wasm)──▶ site/static/nb/chNN_slug.html
                                        │                  (+ ONE shared assets/ and public/)
                                        │       Hugo shortcode {{< marimo src="/nb/chNN_slug.html" >}}
                                        ▼                       ▼
                              reader's browser boots Pyodide, micropip-installs
                              the shared wheel, and runs the notebook locally
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

`notebooks/public/` also carries chapter-specific precomputed datasets
(currently `ch16_baroclinic-instability.npz`, 5 MB). Because every chapter
exports into one directory (step 3), marimo emits a single shared `public/`,
so each file there is stored exactly once for the whole book. Nothing in it
is preloaded: an exported page references `public/` only from the notebook
source it embeds, so a file is fetched over HTTP when — and only when — the
chapter that wants it asks for it. A dataset used by one chapter therefore
costs the other 23 nothing.

This was not true under the older per-chapter export form, where `public/`
was copied into all 24 self-contained bundles and a single 5 MB file quietly
added ~125 MB to the built site. That is why a `notebooks/data/` directory
used to exist; the shared export form removed the need for it.

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
marimo export html-wasm notebooks/chNN_slug.py -o site/static/nb/chNN_slug.html --mode run
```

**The `-o` form is load-bearing.** Naming a *file* makes its parent the output
directory, so all 24 chapters share one `assets/` tree (~687 files, 27 MB —
the marimo frontend and the Pyodide runtime) and one `public/`. The whole
built book is then ~35 MB. Naming a *directory* instead
(`-o site/static/nb/chNN_slug/`) gives each chapter a self-contained ~27 MB
copy of that identical tree: ~650 MB across the book, and the reader
re-downloads the runtime on every chapter instead of hitting a warm cache.

Two consequences the Makefile handles, and that must survive any edit to it:

- `make clean` runs before a full export. marimo *merges* into an existing
  `assets/` (`shutil.copytree(..., dirs_exist_ok=True)`), so chunks from an
  older marimo version would otherwise accumulate forever.
- With one shared `assets/`, every export rehashes every chapter's chunk
  filenames, so a reader holding a stale cache gets "Failed to fetch
  dynamically imported module". `scripts/patch_chunk_reload.py` injects a
  handler that reloads once, guarded by a `sessionStorage` flag against
  reload loops.

marimo also drops a stray copy of `CLAUDE.md` into the output directory; the
Makefile deletes it, so the authoring instructions are not published.

`--mode run` gives the reader a read-and-interact app (sliders work, code is
visible but not editable). Use `--mode edit` for chapters where the reader
should modify code live. `site/static/nb/` is git-ignored — bundles are
build products, rebuilt from source.

### 4. Hugo embeds the bundle

Each chapter page in `site/content/partN/chNN_slug.md` embeds its notebook
with the `marimo` shortcode:

```
{{< marimo src="/nb/chNN_slug.html" >}}
```

`make serve` previews the whole site (`hugo server -D` from `site/`).
Deployment is automatic: merging to `main` runs `.github/workflows/deploy.yml`,
which re-exports every notebook, builds the Hugo site and publishes
`site/public/` to GitHub Pages at https://anigfd.github.io/. Nothing built is
ever committed.

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
  integrating live — put the file in `notebooks/public/<name>.npz` and fetch
  it with `mo.notebook_location() / "public" / ...` (see above).
- **Lazy-load**: Pyodide does not boot until the reader clicks Run.

## Testing, CI and deployment

`.github/workflows/deploy.yml` runs on every push to `main` and every pull
request, and does the whole pipeline in order:

1. `make test` — `pytest -q tests/`, numeric correctness of every `gfdlib`
   primitive (analytic solutions, conservation laws, convergence orders).
   A hard gate: nothing ships if the numerics are wrong. It must also pass
   locally before every commit.
2. `make notebooks` — exports every chapter. A fully-reactive notebook
   executes end to end during the WASM export, so a clean export is a real
   integration test, not just a parse check.
3. A shape check on the result. The per-chapter export form also "succeeds"
   while costing ~27 MB per chapter, so CI asserts what a correct build looks
   like: exactly one `assets/` tree, the wheel and ch16's `.npz` in the shared
   `public/`, an `.html` per chapter each carrying the chunk-reload handler,
   and no stray `CLAUDE.md`.
4. `hugo --minify --source site`, then `actions/upload-pages-artifact`.

A separate `deploy` job, gated on `main`, publishes that artifact with
`actions/deploy-pages`. A pull request therefore gets the tests, the export
and the shape check, but never a deploy. Nothing generated is committed — the
published site is rebuilt from source every time.

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
