---
title: "Ch 15. Barotropic instability"
weight: 515
part: 5
---

## Overview

Whether a shear flow spontaneously rolls up into vortices or just sits
there, however hard you perturb it, is decided entirely by the *curvature*
of its own velocity profile. Rayleigh's inflection-point criterion — later
extended to a rotating planet by Kuo — says instability is impossible
unless the meridional gradient of absolute vorticity changes sign
somewhere in the flow. This chapter turns that criterion into a genuine
growth-rate calculator (a new kind of numerics for this book: a matrix
eigenvalue problem, not a time-stepping PDE), then hands the fastest-growing
mode to the same nonlinear machinery from Ch. 18 to watch it roll up into a
chain of Kelvin-Helmholtz vortices.

## The model

Notebook: `notebooks/ch15_barotropic-instability.py` → exported to `/nb/ch15_barotropic-instability/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$(U-c)(\phi''-k^2\phi)+(\beta-U'')\phi=0.$$

Part A solves this as a generalized eigenvalue problem, rewritten
$B^{-1}A$ and diagonalized with plain `numpy.linalg.eigvals` — no SciPy.
Part B reuses `gfdlib.spectral.Grid` and `gfdlib.timestep.ifrk4_step`
unchanged from ch18, with a new initial condition,
`gfdlib.instability.double_shear_layer`. The solver is independently
verified against two facts that don't depend on this code being right in
advance: a profile with no inflection point gives *exactly* zero growth at
every wavenumber tested, and the classical $\tanh(y/\delta)$ shear layer's
growth-rate peak lands at $k\delta\approx0.445$, matching the published
value (Michalke 1964, $k\delta\approx0.4446$).

{{< marimo src="/nb/ch15_barotropic-instability/" >}}

## Both fluids

- **Atmosphere:** the polar-front jet and subtropical jet both meander and
  shed cutoff eddies via barotropic instability acting on their horizontal
  shear — the same inflection-point mechanism as this notebook's shear
  layer, just embedded in a jet rather than isolated.
- **Ocean:** western boundary currents (the Gulf Stream, Kuroshio) are
  strongly barotropically unstable, shedding rings and meanders at a
  wavelength set by the same growth-rate-curve logic computed here.

## Exercises

1. *(analytic)* Derive the Rayleigh-Kuo equation from the linearized
   barotropic vorticity equation, starting from $\zeta_t+u\zeta_x+v\zeta_y+
   \beta v=0$ linearized about $U(y)$. Confirm your equation matches the
   classical $\beta=0$ (Rayleigh) form in Drazin & Reid.
2. *(computational)* Modify `ch15_barotropic-instability.py`'s Part A to
   add a fourth profile, a *piecewise-linear* jet (a well-known
   textbook case with an exact analytic growth-rate formula). Compare your
   numerical curve to the analytic one.
3. *(exploratory)* In Part A, find the value of $\beta$ (for the jet
   profile) at which the growth-rate curve's maximum first drops to zero —
   the flow's stabilization threshold.

## Further reading

Drazin, P. G. & Reid, W. H., *Hydrodynamic Stability*, ch. 4; Michalke, A.
(1964), *J. Fluid Mech.* **19**, 543–556; Kuo, H.-L. (1949), *J. Meteorol.*
**6**, 105–122; Vallis, *AOFD* 2nd ed., §9.2–9.3.
