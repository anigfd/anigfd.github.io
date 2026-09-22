# Building *An interactive Geophysical Fluid Dynamics Textbook*
### A plan for a browser-runnable GFD textbook (marimo + WebAssembly)

**Author target site:** `anigfd.github.io` — Hugo, notebooks exported to HTML, run in-browser via Pyodide/WASM.
**Audience:** first-year Ph.D. students in atmospheric & oceanic sciences, or applied-math / physics students focused on ocean–atmosphere physics.
**Assumed background:** multivariable calculus & linear algebra; a first course in fluid mechanics or classical mechanics; comfort with Python/NumPy; no prior GFD required.
**Pedagogical spine:** Vallis-style — systematic scaling, balance, and instability arguments, physically grounded, connected to runnable models. Applications *balanced and fluid-first*: teach each mechanism generically, then instantiate it in **both** an atmospheric and an oceanic example. A later "Structure" part carries the Salmon-flavored variational/Hamiltonian material as an elegant capstone rather than a prerequisite.

---

## 1. What this book is, and what makes it different

There is no shortage of excellent GFD texts (Vallis, Salmon, Pedlosky, Gill, Cushman-Roisin & Beckers, Holton & Hakim, Olbers–Willebrand–Eden). What none of them can do on paper is let the reader **turn the knob and watch the balance break**. Your four existing notebooks already prove the format works. The book's thesis:

> Every core GFD concept has a *minimal runnable model* — small enough to execute in a browser tab, rich enough to exhibit the real phenomenon. The text derives the model; the notebook lets the student falsify their intuition against it.

Three design commitments follow from that thesis:

1. **One idea, one model, one knob.** Each notebook isolates a single mechanism and exposes the 1–3 dimensionless parameters that govern it (Ra, Ro, Bu, Fr, Re, Ek, β...). The student's job is to find the transition.
2. **Derivation and code share notation.** The equation printed in the text is the equation stepped in the solver. No silent nondimensionalization gaps.
3. **The book is a graph, not a line.** A "physics-first" reader path and a "math-structure" reader path both exist; notebooks are nodes reachable from either.

### Your four existing notebooks, situated
| Existing notebook | Where it belongs in the arc | Role |
|---|---|---|
| Rayleigh–Bénard → Lorenz (LBM D2Q9) | Part V, Convection & low-order chaos | Capstone on convection + gateway to dynamical-systems thinking |
| Geostrophic adjustment & shallow water (pseudo-spectral RK4, β-plane) | Part II/III hinge, Rotating shallow water | The central "balance vs. waves" notebook — arguably the spine of the whole book |
| 2D Navier–Stokes turbulence (vorticity–streamfunction) | Part VI, Geostrophic turbulence | Dual-cascade centerpiece |
| Internal gravity waves (stratified) | Part IV, Stratified flow & waves | Anisotropic dispersion showcase |

They span four of the six major parts already — you're further along than a chapter count suggests. The gaps are the *connective* chapters (equations of motion, vorticity/PV, QG, instabilities) and the *closing* chapters (wind-driven & buoyancy-driven circulation, wave–mean interaction, structure).

---

## 2. Book architecture

Six parts, ~24 chapters. Each chapter = short text (web pages) + **one primary interactive notebook** + optional "mini" widgets. Difficulty rises left-to-right; the four existing notebooks are marked ★.

### Part I — Foundations (the rotating, stratified fluid)
1. **Kinematics & the material derivative** — fields, trajectories vs. streamlines, Reynolds transport. *Notebook:* particle advection in a prescribed 2D flow; watch tracer filamentation. *Knob:* strain vs. rotation of the velocity-gradient tensor (Okubo–Weiss).
2. **Governing equations & Boussinesq approximation** — mass, momentum, thermodynamics; anelastic vs. Boussinesq; buoyancy. *Notebook:* interactive scaling calculator — enter a length/velocity/latitude, get Ro, Ek, Fr, Re, Bu and which terms survive.
3. **Effects of rotation** — rotating frame, Coriolis & centrifugal, the tangent plane, f- and β-planes. *Notebook:* inertial oscillations & the Taylor–Proudman column; slider from Ro≫1 to Ro≪1.
4. **Scaling & the dimensionless numbers** — formal nondimensionalization, distinguished limits. *Notebook:* the "regime map" — a Ro–Bu plane you click to see which reduced model applies.

### Part II — Balance and its adjustment
5. **Geostrophic & hydrostatic balance** — thermal wind, the balance hierarchy. *Notebook:* thermal-wind reconstruction — draw a temperature field, get the geostrophic shear.
6. **Rotating shallow water** — the RSW equations, PV, energy. *Notebook:* ★ **geostrophic adjustment** (your existing one, extended with a PV-conservation diagnostic and a Rossby-radius slider).
7. **Vorticity and potential vorticity** — vorticity equation, Ertel & shallow-water PV, invertibility. *Notebook:* PV staircases / PV inversion — place PV anomalies, invert for the flow.

### Part III — Quasi-geostrophy (the workhorse)
8. **The QG approximation** — asymptotics from RSW, the QGPV equation. *Notebook:* single-layer QG on a β-plane — inject a vortex, watch Rossby-wave radiation.
9. **Rossby waves** — dispersion, group velocity, β-refraction, WKB rays. *Notebook:* Rossby-wave ray tracer on a background flow; reproduce the "great circle" ray paths.
10. **Two-layer QG & available potential energy** — baroclinic/barotropic modes. *Notebook:* two-layer QG spin-up.

### Part IV — Stratified flow and waves
11. **Stratification & buoyancy frequency** — N², internal modes. *Notebook:* vertical-mode solver — enter N²(z), get the eigenmodes.
12. **Internal gravity waves** — anisotropic dispersion, beams, reflection, critical layers. *Notebook:* ★ **internal gravity waves** (existing; add a variable-N² and a critical-layer demo).
13. **Wave–mean-flow interaction (intro)** — Stokes drift, radiation stress, the transformed Eulerian mean at a first pass. *Notebook:* wave packet on a shear; momentum deposition at the critical layer.

### Part V — Instabilities and the routes to disorder
14. **Convection & Rayleigh–Bénard** — linear onset, the Lorenz truncation, chaos. *Notebook:* ★ **Rayleigh–Bénard → Lorenz** (existing; add the bifurcation diagram as an interactive sweep).
15. **Barotropic instability** — Rayleigh/Kuo criteria, the shear layer. *Notebook:* growth-rate solver + nonlinear roll-up of a jet.
16. **Baroclinic instability** — the Eady & Charney problems, energetics. *Notebook:* Eady growth-rate calculator + a 2-layer life-cycle run — the single most important instability in the book, shown in **both** storm-track and ocean-eddy settings.
17. **Symmetric & other instabilities (survey)** — inertial, symmetric, Kelvin–Helmholtz. *Notebook:* the Ri–Ro stability map.

### Part VI — Turbulence and the general circulation
18. **Geostrophic (2D) turbulence** — dual cascade, spectra, jets. *Notebook:* ★ **2D Navier–Stokes turbulence** (existing; add a β-plane switch to grow zonal jets and a spectrum overlay).
19. **Eddy transport & mixing** — diffusivities, PV mixing, the eddy-mean decomposition. *Notebook:* passive-tracer stirring + effective diffusivity diagnostic.
20. **Wind-driven ocean circulation** — Sverdrup, Stommel, Munk; western boundary currents. *Notebook:* barotropic gyre solver — dial wind stress & friction, get the Gulf-Stream-like WBC.
21. **Buoyancy-driven / overturning circulation** — Hadley cell, MOC, abyssal recipes. *Notebook:* a minimal overturning box/2D model.

### Part VII — Structure (the Salmon capstone, optional path)
22. **Hamiltonian & variational GFD** — least action, particle relabeling, why PV is conserved. *Notebook:* symplectic vs. non-symplectic integrators on a point-vortex system — watch energy drift.
23. **Conservation laws & wave activity** — Noether's theorem in GFD, pseudomomentum, non-acceleration. *Notebook:* wave-activity budget in a shear.
24. **Balanced models & the slow manifold** — a unifying retrospective; where each reduced model lives.

> Part VII is deliberately terminal and skippable. It rewards the reader who wants the *why behind the why*, in Salmon's spirit, without gating the physically-motivated main path.

---

## 3. Shared infrastructure (build this once, reuse everywhere)

The difference between "four nice notebooks" and "a textbook" is consistency. Four things to standardize up front:

**3.1 A common notation & symbol table.** One page, linked from every chapter. Fix conventions early: sign of N², direction of β, nondimensionalization scheme, layer indexing. Nothing erodes a multi-author-feel faster than f₀ meaning three things.

**3.2 A tiny shared library, `gfdlib` (pure NumPy/Pyodide-safe).** Not a framework — a handful of vetted primitives every notebook imports so solvers don't get re-implemented four different ways:
- spectral helpers: `rfft2`/`irfft2` wrappers, dealiasing (2/3 rule), wavenumber grids, Poisson/`∇⁻²` inversion, Jacobian `J(ψ,ζ)` (Arakawa or spectral);
- time steppers: RK4, IMEX for stiff diffusion, leapfrog+Robert filter;
- diagnostics: energy/enstrophy/PV budgets, isotropic spectra, Nusselt/Sverdrup numbers;
- plotting: a consistent colormap policy (diverging for signed fields, perceptually-uniform sequential otherwise) and a standard figure frame.

Ship it as a single `.py` file the notebooks fetch, or inline-vendored per notebook for zero-dependency robustness (WASM sandboxes can't always pip-install).

**3.3 A notebook template.** Every notebook opens the same way: *Physical question → the equations (same symbols as the text) → the numerical scheme in one cell → interactive controls → 2–4 diagnostics → "Try this" prompts → "What you should have seen."* marimo's reactive cells make the controls-to-diagnostics wiring clean; standardize the layout so students build muscle memory.

**3.4 An exercises & self-check pattern.** Three tiers per chapter: (a) *analytic* (pencil derivations continuing the text), (b) *computational* (modify the notebook — change a boundary condition, add a term), (c) *exploratory* (open-ended: "find the Ro at which the jet destabilizes"). Optionally, lightweight in-notebook auto-checks (assert a conservation law holds to tolerance).

---

## 4. Technical build & publishing pipeline

Your stack already works; this is about making it repeatable and fast.

**4.1 Authoring.** One marimo notebook per chapter (`.py`, git-tracked — marimo notebooks are plain Python, so they diff and review cleanly, a real advantage over `.ipynb`). Keep prose in the notebook light; the *long-form* exposition lives in Hugo pages, with the notebook embedded/linked. This separation lets you edit text without touching code.

**4.2 Export to WASM.** `marimo export html-wasm notebook.py -o site/ --mode run` (or `--mode edit` to let students edit code live). This produces a self-contained bundle that boots Pyodide in the browser. Automate it: a `Makefile`/CI job that re-exports every changed notebook and drops the HTML into the Hugo static tree.

**4.3 Performance budget (the real constraint).** Pyodide is ~3–10× slower than native and single-threaded by default. Design rules:
- keep grids modest (64²–256² for spectral; the LBM one already shows this is fine);
- precompute wavenumber arrays and FFT plans once, outside the time loop;
- expose a "resolution" and "steps-per-frame" control so the student trades speed for detail;
- for anything genuinely heavy, precompute offline and ship the result as data the notebook visualizes (e.g., a baroclinic life-cycle can be a precomputed field the widget scrubs through);
- lazy-load: don't boot Pyodide until the reader clicks "Run".

**4.4 Site structure (Hugo Blox).** A `Part → Chapter` content tree; each chapter page carries an abstract, the long-form text, the embedded notebook (iframe or marimo island), exercises, and "further reading" pointing at the corresponding Vallis/Salmon/Pedlosky sections. Add: a **dependency graph** landing page (which chapter needs which), a **notation** page, and a **"how to run these"** page.

**4.5 Reproducibility & maintenance.** Pin the marimo and Pyodide versions per release; tag the repo per "edition." A CI check that every notebook exports without error catches breakage before it ships. `requirements`-in-notebook comments so a reader can also run locally.

---

## 5. Suggested build order (phased)

The order optimizes for **shippable value early** and **reusable infrastructure first**.

**Phase 0 — Foundation (infrastructure, ~1 sprint).** Build `gfdlib`, the notebook template, the notation page, the export Makefile/CI, and retrofit your four existing notebooks onto the template + shared library. Deliverable: four consistent notebooks + a build you can run with one command. *This makes every later chapter cheaper.*

**Phase 1 — The spine (Part II–III).** Chapters 5–9. Geostrophic balance → RSW (you have adjustment) → vorticity/PV → QG → Rossby waves. This is the conceptual core of GFD and connects two of your existing notebooks. After this phase the book is already a coherent short course.

**Phase 2 — Instabilities (Part V).** Chapters 14–16, anchored by baroclinic instability (Ch. 16) — the highest-value new notebook in the whole book, and the natural bridge to your turbulence notebook. Reuse the convection→Lorenz notebook here.

**Phase 3 — Waves & stratification (Part IV).** Chapters 11–13 around your internal-wave notebook.

**Phase 4 — Circulation & turbulence (Part VI).** Chapters 18–21; your turbulence notebook plus the wind-driven gyre.

**Phase 5 — Foundations & structure (Parts I, VII).** Chapters 1–4 and 22–24. Foundations are written last on purpose: once the later notebooks exist, the intro can point forward concretely ("you'll use this in Ch. 16").

Rationale: you already own notebooks in Parts II, IV, V, VI. Phasing the *connective* chapters first turns a scattered set into a spine as fast as possible.

---

## 6. Risks & mitigations

| Risk | Mitigation |
|---|---|
| WASM too slow for spectral solvers at useful resolution | Resolution slider + precomputed heavy runs; the existing pseudo-spectral notebook already proves 128²-class works |
| Notebooks drift apart in style/notation | `gfdlib` + template + notation page enforced from Phase 0 |
| Scope creep (24 chapters is a lot) | Parts I & VII are explicitly optional; a viable "minimum book" is Phases 0–2 (spine + instabilities) |
| Pyodide/marimo version breakage over time | Pin versions per edition; CI export check |
| Text competing with Vallis instead of complementing | Position as *the runnable companion*: every chapter cites the paper text it pairs with, adds the model they can't |

---

## 7. Minimum viable book vs. full book

- **MVB (recommended first public milestone):** Phase 0 + Phase 1 + Ch. 16 (baroclinic instability). Nine notebooks, a coherent spine from balance through the central instability, all four existing notebooks integrated. This is already a citable, useful resource.
- **Full book:** all six/seven parts, ~24 notebooks.

---

## 8. Immediate next actions (what I can help you build now)

1. **Scaffold `gfdlib`** — I can draft the shared spectral/time-stepping/diagnostics module as a Pyodide-safe `.py` and a smoke-test notebook.
2. **Draft the notebook template** — a marimo template with the standard section layout and the controls→diagnostics wiring.
3. **Prototype the highest-value new notebook** — the baroclinic-instability (Eady + 2-layer life cycle) notebook, since it anchors Phase 2 and links your existing turbulence and adjustment notebooks.
4. **Write a chapter** — pick any chapter above and I'll draft the long-form text (Vallis-grounded, with the both-fluids examples) plus its notebook spec.
5. **Set up the export pipeline** — a `Makefile` + CI config that turns `notebooks/*.py` into WASM HTML in your Hugo tree.

Tell me which of these to start on and I'll build it.
