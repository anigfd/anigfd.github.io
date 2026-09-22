# CLAUDE.md — An interactive Geophysical Fluid Dynamics Textbook

You are helping build an **interactive GFD textbook**: long-form text on a Hugo
site + runnable **marimo** notebooks exported to **WebAssembly** (Pyodide) so
readers execute every model in the browser, no install.

- Live site: https://anigfd.github.io/
- Audience: first-year AOS / applied-math / physics PhD students.
- Pedagogy: **Vallis-style** (systematic scaling, balance, instability;
  physically grounded). Applications **balanced & fluid-first** — teach each
  mechanism generically, then show BOTH an atmospheric and an oceanic instance.
  Salmon-flavored variational/Hamiltonian material lives ONLY in Part VII and is
  optional/skippable — never a prerequisite.

## The one rule that keeps this a *textbook* and not a notebook gallery
Every notebook imports vetted primitives from **`gfdlib/`** (pure NumPy,
Pyodide-safe — NO SciPy, NO compiled deps). Do not re-implement a spectral
inversion, a Jacobian, a time-stepper, or a spectrum four different ways. If you
need a new primitive, add it to `gfdlib` with a docstring and a test in
`tests/`, then use it. Notation and sign conventions are fixed in `NOTATION.md`
and must match between the printed equation and the stepped equation.

## Repository layout
```
gfdlib/            # shared Pyodide-safe primitives (spectral, timestep, diagnostics, plotting)
notebooks/         # one marimo .py per chapter: chNN_slug.py  (plain Python — diffs cleanly)
tests/             # pytest: numeric correctness of gfdlib (run before every commit)
site/content/      # Hugo pages: part/chapter tree, long-form text
site/static/nb/    # exported WASM bundles land here (git-ignored; produced by `make`)
Makefile           # `make notebooks` exports every changed .py -> site/static/nb/
NOTATION.md        # symbol table + sign conventions (linked from every chapter)
docs/              # contributor docs: architecture.md, gfdlib.md (API tour), authoring.md
```

## Build / preview / test commands
- Export all notebooks to WASM:  `make notebooks`
  (first builds the gfdlib wheel into `notebooks/public/`; marimo copies that
  folder into the export, and notebooks micropip-install the wheel when
  running under Pyodide — see the import cell in `notebooks/_template.py`)
- Export one:  `make nb-one NB=ch06_geostrophic-adjustment`
  (use `--mode edit` for chapters where the reader should edit code live)
- Preview site:  `make serve`  (runs `hugo server -D` from `site/`)
- Run numeric tests:  `make test`  (pytest on `tests/` — MUST pass before commit)
- Deploy: automatic. Merging to `main` runs `.github/workflows/deploy.yml`,
  which re-exports every notebook, builds the site and publishes to GitHub
  Pages. Nothing built is ever committed.

## Build facts that are easy to break
- Chapters export as `-o site/static/nb/chNN_slug.html` — a **file**, not a
  directory. That form makes the parent the output directory, which is what
  gives one shared `assets/` and one shared `public/`. Switching to
  `-o .../chNN_slug/` silently multiplies the built site by ~27 MB per
  chapter (33 MB → ~650 MB) and makes readers re-download Pyodide on every
  chapter instead of hitting a warm cache.
- `make clean` must run before a full export: marimo *merges* into an
  existing `assets/` rather than replacing it, so stale chunks accumulate.
- With one shared `assets/`, every export rehashes every chapter's chunk
  filenames, so a stale browser cache breaks the page.
  `scripts/patch_chunk_reload.py` injects a reload-once handler; the Makefile
  runs it. Do not drop it.
- `marimo export` drops a stray `CLAUDE.md` into the output directory; the
  Makefile removes it. It must not reach the published site.
- Precomputed chapter data goes in `notebooks/public/<name>.npz`, fetched
  with `mo.notebook_location() / "public" / ...`. There is one shared
  `public/`, but nothing in it is preloaded — only the chapter that asks for
  a file downloads it.
- The Hugo shortcode takes the file: `{{< marimo src="/nb/chNN_slug.html" >}}`.

## Performance budget (Pyodide is single-threaded, ~3-10x slower than native)
- Spectral grids 64^2-256^2. The included solver does 128^2 x 400 steps in <1s native.
- Precompute wavenumber arrays / FFT setup ONCE, outside the time loop (build a
  `gfdlib.spectral.Grid` at notebook top and reuse it).
- Expose a **resolution** slider and a **steps-per-frame** slider so the reader
  trades speed for detail.
- For genuinely heavy runs (e.g. a baroclinic life-cycle), precompute offline and
  ship the field as data the widget scrubs — don't integrate it live.
- Lazy-load: don't boot Pyodide until the reader clicks Run.

## Notebook template (every chapter follows it)
Copy `notebooks/_template.py`. Standard section order:
1. Title + the physical question (markdown)
2. The governing equations — SAME symbols as NOTATION.md
3. The numerical scheme in ONE cell, importing from gfdlib
4. Interactive controls (marimo UI: sliders for the 1-3 dimensionless numbers)
5. 2-4 diagnostics (fields + a budget/spectrum) using gfdlib.plotting + diagnostics
6. "Try this" exploratory prompts
7. "What you should have seen" (the expected result)

## House conventions
- matplotlib: always `fig, ax = ...` then `ax.plot(...)`; NEVER bare `plt.plot`/`plt.savefig`.
- Colormaps: diverging+symmetric for SIGNED fields (vorticity, PV, height anomaly);
  perceptually-uniform sequential for non-negative fields. Use `gfdlib.plotting`.
- Every chapter page ends with "Further reading" pointing at the matching
  Vallis / Salmon / Pedlosky sections.
- Git: work on a branch, never commit straight to main; `make test` must pass first.

## Chapter map
See `PLAN.md` (the full curriculum). Existing notebooks to integrate onto the
template + gfdlib: geostrophic adjustment (ch06), 2D turbulence (ch18),
internal gravity waves (ch12), Rayleigh-Benard->Lorenz (ch14).
