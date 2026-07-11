---
title: "Ch 11. Stratification & the buoyancy frequency"
weight: 411
part: 4
---

## Overview

Ch. 10 split a two-LAYER fluid exactly into barotropic and baroclinic
parts. A continuously stratified fluid admits the same trick, generalized:
any $N^2(z)$ profile has its own complete, orthogonal set of **vertical
normal modes**, each behaving as an independent single-layer QG system with
its own deformation radius. This chapter solves for those modes directly,
and shows the two-layer model used everywhere earlier in the book is
exactly the crudest possible truncation of this richer theory — mode 0 and
mode 1 only.

## The model

Notebook: `notebooks/ch11_stratification.py` → exported to
`/nb/ch11_stratification/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}).

$$\frac{d}{dz}\!\left[\frac{1}{N^2(z)}\frac{d\Phi_n}{dz}\right]
  +\frac{1}{c_n^2}\Phi_n=0,\qquad \Phi_n'(0)=\Phi_n'(H)=0,\qquad L_n=\frac{c_n}{f_0}.$$

`gfdlib.stratification.vertical_modes` discretizes this Sturm-Liouville
problem by the finite-volume method and solves it as a genuinely symmetric
eigenproblem (`numpy.linalg.eigh`), validated against the exact
constant-$N$ solution $\Phi_n=\cos(n\pi z/H)$, $c_n=N_0H/(n\pi)$.

{{< marimo src="/nb/ch11_stratification/" >}}

## Both fluids

- **Ocean:** the first baroclinic deformation radius, computed from a
  realistic pycnocline profile, sets the observed $\sim30$–$50\,$km scale
  of mesoscale eddies at midlatitudes — smaller than the constant-$N$
  formula would suggest, because real stratification concentrates in a
  thin near-surface layer rather than spreading uniformly through the
  water column.
- **Atmosphere:** the tropopause plays the ocean pycnocline's role — a
  region of sharply enhanced $N^2$ near which the first baroclinic mode's
  structure concentrates, setting the $\sim1000\,$km scale of synoptic
  weather systems.

## Exercises

1. *(analytic)* Verify by direct substitution that $\Phi_n=\cos(n\pi z/H)$
   solves $d/dz[(1/N_0^2)d\Phi/dz]+(1/c_n^2)\Phi=0$ with
   $c_n=N_0H/(n\pi)$, and that $\Phi_n'(0)=\Phi_n'(H)=0$ for every integer
   $n$ (matching `gfdlib.stratification`'s docstring and the notebook's
   built-in validation panel).
2. *(computational)* Modify `ch11_stratification.py` to add a two-pycnocline
   profile (a shallow seasonal thermocline stacked above a deeper permanent
   one) and see whether a genuinely new SHAPE of mode 2 emerges — does it
   develop an extra local extremum near the shallow feature?
3. *(exploratory)* Using the pycnocline profile, sweep the thickness slider
   from thick to thin at fixed $N_{min},N_{max}$, and record $L_1$ at each
   setting. Does $L_1$ approach a two-layer-like limit as the pycnocline
   sharpens? Compare to a matching Ch. 16 $F=f_0^2/(g'H)$ calculation with
   $g'$ estimated from your $N_{max},$ thickness.

## Further reading

Gill, A. E. (1982), *Atmosphere-Ocean Dynamics*, §6.11-6.12 (vertical
normal modes); Chelton, D. B. et al. (1998), *J. Phys. Oceanogr.* **28**,
433–460 (global deformation-radius climatology); Vallis, *AOFD* 2nd ed.,
§5.5-5.6; Pedlosky, *Geophysical Fluid Dynamics* 2nd ed., §6.9-6.11.
