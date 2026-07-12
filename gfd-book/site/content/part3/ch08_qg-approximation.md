---
title: "Ch 8. The quasi-geostrophic approximation"
weight: 308
part: 3
---

## Overview

Quasi-geostrophy is what you get when you take rotating shallow water
seriously in the limit every large-scale flow actually lives in: small
Rossby number. The leading-order flow is exactly geostrophic; the whole
content of the dynamics is in the small ageostrophic correction, and QG
theory shows that correction is completely captured by conserving a single
scalar, the quasi-geostrophic potential vorticity. This chapter derives
that equation from Ch. 6's RSW, and uses an isolated vortex to show what's
new relative to Ch. 7's barotropic PV: a deformation-radius stretching term
that screens vortices, caps the Rossby-wave dispersion relation, and sends
an isolated eddy drifting westward while it sheds a trailing wake — the
*beta-gyre*.

## The model

Notebook: `notebooks/ch08_qg-approximation.py` → exported to `/nb/ch08_qg-approximation/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$q_t+J(\psi,q)=0,\qquad q=\nabla^2\psi-\frac{\psi}{L_R^2}+\beta y.$$

Reuses Ch. 18's pseudo-spectral machinery and Ch. 6's Helmholtz PV
inversion directly — nondimensionalized by $L_R$, the QGPV balance operator
$(\nabla^2-1)\psi=q-\beta y$ is *the same operator* as the shallow-water PV
inversion, just now driving a fully nonlinear evolution. No new inversion
code exists or is needed; `gfdlib.qg` supplies only the modified
Rossby-wave dispersion relation.

## What to look for

1. **Screening.** Set $\beta=0$ and look at the initial streamfunction: a
   vortex of radius $\sigma\lesssim L_R$ has its far-field velocity cut off
   exponentially beyond $L_R$ (the $-\psi/L_R^2$ term turns the inversion
   from a long-range logarithm into a short-range Yukawa-like kernel).
   Compare against Ch. 7, where the same blob's influence extended across
   the whole domain.
2. **Westward drift from nothing.** Turn $\beta$ up: an isolated,
   perfectly symmetric vortex begins to translate westward without being
   pushed. Watch the first frames closely — advection of planetary
   vorticity builds an antisymmetric secondary dipole (the *beta-gyre*)
   across the vortex, and it is that dipole which carries the parent
   along at a speed of order $\beta L_R^2$.
3. **The wake.** As the vortex drifts it sheds a trailing Rossby-wave
   wake. Check the dispersion-relation panel: the wake's dominant
   wavelength sits near the wavenumber whose westward phase speed matches
   the vortex's own drift — the wave the vortex can resonate with.
4. **The capped dispersion curve.** In the dispersion panel, note that
   $|\omega|$ has a maximum near $k\sim1/L_R$: unlike Ch. 7's barotropic
   waves, QG Rossby waves cannot propagate faster than $\beta L_R^2$, and
   long waves all bunch at that speed — why oceanic eddies of many sizes
   drift westward at nearly the same rate.

{{< marimo src="/nb/ch08_qg-approximation/" >}}

## Both fluids

- **Atmosphere:** synoptic-scale weather systems are quasi-geostrophic to
  good approximation — QG theory (via the omega equation, not covered here)
  is still the standard first-pass tool for diagnosing where a system will
  intensify.
- **Ocean:** mesoscale eddies (the ocean's weather, $L\sim50$–150 km) are
  observed to drift predominantly westward at a speed set by $\beta L_R^2$
  — precisely the beta-gyre mechanism this notebook reproduces from a single
  isolated vortex.

## Exercises

1. *(analytic)* Starting from $q=\nabla^2\psi-\psi/L_R^2+\beta y$, show that
   in the limit $L_R\to\infty$ (or equivalently rescaling lengths $\gg L_R$)
   QG reduces exactly to Ch. 7's barotropic PV. What physically has been
   assumed away in that limit?
2. *(computational)* Modify `ch08_qg-approximation.py` to place two vortices
   of opposite sign, each individually beta-drifting. Do they still form a
   self-propelling dipole (Ch. 7), or does the beta-drift dominate?
3. *(exploratory)* Using the dispersion-relation panel, predict which
   vortex radius $\sigma$ (relative to $L_R$) should shed the most visible
   wake, then confirm it in the snapshot browser.

## Further reading

Vallis, *AOFD* 2nd ed., §5.1–5.5 (QG scaling and the QGPV equation);
Pedlosky, *Geophysical Fluid Dynamics*, ch. 6; Chelton et al. (2011),
*Progress in Oceanography* **91**, 167–216 (observed westward propagation
of mesoscale eddies); Flierl (1987), *J. Phys. Oceanogr.* (beta-drift of
isolated vortices).
