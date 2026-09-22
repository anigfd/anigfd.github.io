# Authoring guide — adding a chapter

The checklist for going from "chapter N should exist" to a merged PR.
Conventions here are load-bearing: they are what keep 20+ notebooks feeling
like one book. See [`architecture.md`](architecture.md) for how the build
works and [`gfdlib.md`](gfdlib.md) for what the library already provides.

## 0. Before writing anything

- Read the chapter's entry in [`PLAN.md`](../PLAN.md) — it specifies the
  physical question, the notebook concept, and the knob (the 1–3
  dimensionless parameters the reader varies).
- Read [`NOTATION.md`](../NOTATION.md). The symbols and sign conventions
  there are mandatory; if the chapter needs a new symbol, add it to
  NOTATION.md in the same PR.
- Check [`gfdlib.md`](gfdlib.md) for primitives you can reuse. Chapters
  routinely share machinery (ch. 21 reuses ch. 14's convection solver
  wholesale).

## 1. Branch

```bash
git checkout main && git pull
git checkout -b chNN-short-slug
```

Never commit to main directly. One chapter per branch/PR.

## 2. New primitives go in `gfdlib` first

If the chapter needs numerics that don't exist yet:

- add them to an existing module or a new `gfdlib/<topic>.py` — **pure
  NumPy only** (no SciPy, no compiled deps; Pyodide must run it);
- write a docstring stating the equation solved and its conventions;
- write tests in `tests/test_gfdlib.py` that check *correctness*, not
  execution: an analytic solution, a conservation law held to tolerance,
  or a measured convergence order;
- register a new module in `gfdlib/__init__.py` (docstring example line,
  the `from . import ...` line, and `__all__`).

Heads-up: because every chapter branch appends to the same import list and
test file, two open chapter PRs will conflict there. The conflicts are
purely additive — resolve by keeping both sides' additions.

## 3. Write the notebook

Copy `notebooks/_template.py` to `notebooks/chNN_slug.py`. Standard section
order (the template has the scaffolding):

1. **Title + the physical question** — what the reader can do/explain after
   this chapter that they couldn't before.
2. **Governing equations** — same symbols as NOTATION.md. The printed
   equation must be the stepped equation.
3. **The numerical scheme in ONE cell** — imports from `gfdlib`, all
   precomputation (`spectral.Grid`, masks, operators) outside the time loop.
4. **Interactive controls** — marimo sliders for the dimensionless numbers;
   include a resolution and/or steps-per-frame slider for stepped models.
5. **2–4 diagnostics** — fields plus at least one budget or spectrum, drawn
   with `gfdlib.plotting`.
6. **"Try this"** — exploratory prompts (find the transition, break the
   balance).
7. **"What you should have seen"** — the expected result, stated plainly.

House style:

- matplotlib: always `fig, ax = plt.subplots(...)`, never bare `plt.plot`;
- diverging + symmetric-limit colormaps for signed fields, sequential for
  non-negative ones (use `gfdlib.plotting`);
- keep the performance budget: grids 64²–256², precompute once, lazy
  heavy work behind a Run button.

Iterate locally with `marimo edit notebooks/chNN_slug.py`.

## 4. Write the chapter page

Create `site/content/partN/chNN_slug.md` (front matter: `title`, `weight`,
`part` — match the neighbors). Standard sections: **Overview** (the physical
story), **The model** (equations + what the gfdlib functions do), the
embedded notebook via `{{< marimo src="/nb/chNN_slug/" >}}`, **Both fluids**
(one atmospheric and one oceanic instance of the mechanism — mandatory),
**Exercises** (three tiers: analytic, computational, exploratory), and
**Further reading** (the matching Vallis / Pedlosky / Salmon sections).

## 5. Verify

```bash
make test        # gfdlib correctness — must pass
make notebooks   # every notebook must export to WASM without error
make serve       # eyeball the chapter page + embedded notebook
```

## 6. PR

Push the branch, open a PR, let CI run (it repeats the test suite and the
export smoke check). Merge only with CI green. If main moved while the PR
was open, expect the additive conflicts from step 2 — merge main, keep both
sides, re-run `make test`, push.
