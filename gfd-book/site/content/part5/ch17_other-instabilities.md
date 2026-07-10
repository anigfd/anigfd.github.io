---
title: "Ch 17. Symmetric & other instabilities"
weight: 517
part: 5
---

## Overview

Chs. 15 and 16 covered the two mechanisms that dominate the large-scale
storm track and ocean mesoscale. This chapter is a **survey** of three
faster, more local mechanisms that matter at fronts and in stratified
shear layers: inertial instability (unbalanced rotation), symmetric
instability (a slantwise hybrid of inertial and buoyant), and
Kelvin-Helmholtz instability (stratified shear, no rotation needed at all).
All three turn out to reduce to comparing a Richardson number against a
Rossby number — one map, three mechanisms, and a genuinely rigorous
theorem (Miles-Howard) governing the third.

## The model

Notebook: `notebooks/ch17_other-instabilities.py` → exported to
`/nb/ch17_other-instabilities/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}):

$$q=(f-\partial U/\partial y)N^2-f(\partial U/\partial z)^2,\qquad
  (U-c)^2(\phi''-k^2\phi)-U''(U-c)\phi+N^2\phi=0.$$

`gfdlib.symmetric` derives the Ertel-PV sign criterion for symmetric and
inertial instability directly (not from a memorized formula); at zero
relative vorticity it reduces to the textbook value $\mathrm{Ri}=1$, and
for $(1+\mathrm{Ro})<0$ the front is unconditionally unstable at every
$\mathrm{Ri}$ — pure inertial instability. `gfdlib.instability.
taylor_goldstein_growth_rate` solves the Taylor-Goldstein equation as a
natively linear $2n\times2n$ generalized eigenvalue problem (keeping the
buoyancy perturbation as an independent unknown, rather than eliminating it
into a harder-to-trust quadratic pencil); at zero stratification it reduces
exactly to ch. 15's barotropic solver. Validating it against the rigorous
Miles-Howard theorem ($\mathrm{Ri}\geq1/4$ everywhere $\Rightarrow$ stable)
caught a real numerical issue: initial growth-rate curves didn't vanish
cleanly above the threshold, traced (not dismissed) to genuine
finite-difference truncation error that is worst exactly at the delicate
$\mathrm{Ri}=1/4$ point, where the equation's critical-layer indicial roots
coalesce — a classically hard point for simple discretizations. The
notebook flags this honestly rather than hiding it.

{{< marimo src="/nb/ch17_other-instabilities/" >}}

## Both fluids

- **Atmosphere:** symmetric instability is the mechanism behind
  banded precipitation structures in fronts and rapidly-intensifying
  cyclones — slantwise convection along sloping absolute-momentum surfaces,
  distinct from ordinary upright thunderstorm convection. Kelvin-Helmholtz
  billows form in strongly-sheared, weakly-stratified layers such as jet
  streaks and clear-air-turbulence zones.
- **Ocean:** symmetric instability is a major mechanism restratifying
  mixed-layer fronts after wintertime convection, competing with
  submesoscale baroclinic instability (ch. 16) to set the depth and
  structure of the mixed layer. Kelvin-Helmholtz billows are directly
  observed at the base of the equatorial undercurrent and in energetic
  tidal-mixing fronts, where they are a leading mechanism for diapycnal
  mixing.

## Exercises

1. *(analytic)* Starting from the Ertel PV $q=(f-\partial U/\partial
   y)N^2-f(\partial U/\partial z)^2$ and the thermal-wind relation
   $f\partial U/\partial z=-\partial b/\partial y$, derive the nondimensional
   form $q/(fN^2)=(1+\mathrm{Ro})-1/\mathrm{Ri}$ used in
   `gfdlib.symmetric`.
2. *(computational)* Modify `ch17_other-instabilities.py`'s Part B to plot
   the peak growth rate (not the whole curve) as a function of $J$ swept
   from $0$ to $0.5$. Confirm it decreases monotonically and reaches zero
   by around $J\approx0.4$, consistent with — but, given the numerical
   caveats documented in the notebook, not a razor-sharp proof of — the
   Miles-Howard $\mathrm{Ri}=1/4$ threshold.
3. *(exploratory)* Compare the two Richardson-number thresholds head to
   head: symmetric instability at $\mathrm{Ro}=0$ requires $\mathrm{Ri}<1$;
   Kelvin-Helmholtz requires $\mathrm{Ri}<1/4$. Why does a rotating,
   geostrophically-balanced parcel tolerate four times less shear before
   destabilizing than a non-rotating one? What extra degree of freedom does
   slantwise motion have that pure vertical overturning does not?

## Further reading

Miles, J. W. (1961), *J. Fluid Mech.* **10**, 496–508 ("On the stability of
heterogeneous shear flows"); Howard, L. N. (1961), *J. Fluid Mech.* **10**,
509–512 ("Note on a paper of John W. Miles"); Hazel, P. (1972), *J. Fluid
Mech.* **51**, 39–61 ("Numerical studies of the stability of inviscid
stratified shear flows"); Hoskins, B. J. (1974), *Q. J. R. Meteorol. Soc.*
**100**, 480–482 ("The role of potential vorticity in symmetric stability
and instability"); Vallis, *AOFD* 2nd ed., §6.9-6.10 (symmetric & inertial
instability), §9.4 (Kelvin-Helmholtz); Pedlosky, *Geophysical Fluid
Dynamics* 2nd ed., §7.13.
