# Interactive Geophysical Fluid Dynamics

Long-form GFD text + runnable **marimo** notebooks exported to **WebAssembly**
(Pyodide) so readers execute every model in the browser. Companion to
Vallis and Salmon; live at https://www.aneeshcs.com/gfd/

## Quick start
```bash
pip install -r requirements.txt          # marimo, numpy, matplotlib, pytest
make test                                # numeric correctness of gfdlib
marimo edit notebooks/_template.py       # author a chapter from the template
make notebooks                           # export all notebooks -> site/static/nb/
make serve                               # preview the Hugo site
```

## How it fits together
- `gfdlib/` — one vetted, Pyodide-safe implementation of each numeric primitive.
- `notebooks/chNN_slug.py` — one marimo notebook per chapter (plain Python).
- `site/` — Hugo content tree; exported WASM bundles land in `site/static/nb/`.
- `NOTATION.md` — symbols & sign conventions (the printed equation == the stepped equation).
- `PLAN.md` — the full curriculum. `CLAUDE.md` — conventions for AI-assisted authoring.

See `CLAUDE.md` before adding a chapter.
