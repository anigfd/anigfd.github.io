# Interactive Geophysical Fluid Dynamics

An interactive GFD textbook: long-form text on a Hugo site plus runnable
**[marimo](https://marimo.io)** notebooks exported to **WebAssembly**
(Pyodide), so readers execute every model directly in the browser — no
install, no server. A runnable companion to Vallis, Salmon, and Pedlosky.

- **Live site:** https://www.aneeshcs.com/gfd/
- **Audience:** first-year PhD students in atmospheric & oceanic sciences,
  applied math, or physics. Assumed background: multivariable calculus,
  linear algebra, a first fluids course, comfort with Python/NumPy —
  no prior GFD.
- **Thesis:** every core GFD concept has a *minimal runnable model* — small
  enough to execute in a browser tab, rich enough to exhibit the real
  phenomenon. The text derives the model; the notebook lets the student
  falsify their intuition against it.

## Quick start

```bash
pip install -r requirements.txt   # marimo, numpy, matplotlib, pytest
make test                         # numeric correctness of gfdlib (must pass)
marimo edit notebooks/ch06_geostrophic-adjustment.py   # run a chapter locally
make notebooks                    # export ALL notebooks -> site/static/nb/
make serve                        # preview the Hugo site (hugo server -D)
```

## Chapter map

Six main parts plus an optional capstone (~24 chapters). Each chapter is a
Hugo page (`site/content/partN/`) wrapping **one primary interactive
notebook** (`notebooks/chNN_slug.py`) built on shared primitives in
`gfdlib/`. ✅ = notebook implemented; 📝 = text-only or planned.

| Part | Chapters | Notebooks | gfdlib modules |
|---|---|---|---|
| **I — Foundations** | 1 kinematics & Okubo-Weiss ✅ · 2 governing eqns/Boussinesq scaling ✅ · 3 rotation (inertial oscillations, Taylor-Proudman) ✅ · 4 Ro-Bu regime map ✅ | 4 of 4 | `kinematics`, `rotation`, `scaling` |
| **II — Balance & adjustment** | 5 geostrophic/thermal wind ✅ · 6 rotating shallow water & adjustment ✅ · 7 vorticity & PV ✅ | 3 of 3 | `balance`, `shallowwater`, `pv` |
| **III — Quasi-geostrophy** | 8 the QG approximation ✅ · 9 Rossby waves & ray tracing ✅ · 10 two-layer QG 📝 | 2 of 3 | `qg`, `rossby` |
| **IV — Stratified flow & waves** | 11 stratification 📝 · 12 internal gravity waves ✅ · 13 wave–mean interaction 📝 | 1 of 3 | `internalwaves` |
| **V — Instabilities** | 14 convection → Lorenz ✅ · 15 barotropic ✅ · 16 baroclinic (Eady + 2-layer) ✅ · 17 symmetric/KH survey ✅ | 4 of 4 | `convection`, `instability`, `baroclinic`, `symmetric` |
| **VI — Turbulence & circulation** | 18 geostrophic turbulence ✅ · 19 eddy transport & mixing ✅ · 20 wind-driven gyres (Stommel/Munk) ✅ · 21 overturning circulation ✅ | 4 of 4 | `diagnostics`, `mixing`, `circulation`, `overturning` |
| **VII — Structure** (optional capstone) | 22–24: Hamiltonian GFD, wave activity, balanced models | 📝 optional | — |

Every notebook also draws on the cross-cutting modules `spectral`
(pseudo-spectral grid), `timestep` (RK4, integrating-factor RK4, leapfrog),
`diagnostics`, and `plotting`. The full curriculum rationale lives in
[`PLAN.md`](PLAN.md).

## Repository layout

```
gfdlib/            # shared, Pyodide-safe numeric primitives (pure NumPy — no SciPy)
notebooks/         # one marimo notebook per chapter: chNN_slug.py (plain Python, diffs cleanly)
tests/             # pytest suite for gfdlib — `make test` must pass before every commit
site/content/      # Hugo pages: the part/chapter tree with the long-form text
site/static/nb/    # exported WASM bundles (git-ignored; produced by `make notebooks`)
docs/              # contributor documentation (architecture, gfdlib API, authoring guide)
NOTATION.md        # the symbol table & sign conventions every chapter links to
PLAN.md            # the full curriculum plan and design rationale
CLAUDE.md          # working conventions for AI-assisted authoring
```

## Documentation

- [`docs/architecture.md`](docs/architecture.md) — how a notebook becomes a
  browser-runnable chapter (marimo → WASM → Hugo), and the performance
  budget that shapes every design decision.
- [`docs/gfdlib.md`](docs/gfdlib.md) — a guided tour of the shared library:
  what each module provides and which chapter uses it.
- [`docs/authoring.md`](docs/authoring.md) — the step-by-step checklist for
  adding a new chapter.
- [`NOTATION.md`](NOTATION.md) — symbols and sign conventions. The equation
  printed in the text is the equation stepped in the solver, always.

## The one design rule

**Notebooks import vetted primitives from `gfdlib/`; they never re-implement
them.** One spectral inversion, one Arakawa/spectral Jacobian, one RK4 —
implemented once, tested in `tests/`, reused everywhere. This is what keeps
24 notebooks a *textbook* rather than a notebook gallery. If a chapter needs
a new primitive, it is added to `gfdlib` with a docstring and a test first.

## Contributing a chapter

Work on a branch, follow [`docs/authoring.md`](docs/authoring.md), and open
a PR. In short: copy `notebooks/_template.py`, use `NOTATION.md` symbols,
import from `gfdlib`, add tests for any new primitive, and make sure
`make test` and `make notebooks` both succeed.
