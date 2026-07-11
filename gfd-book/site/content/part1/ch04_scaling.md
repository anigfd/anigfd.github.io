---
title: "Ch 4. Scaling & the dimensionless numbers"
weight: 104
part: 1
---

## Overview

Ch. 2 computed five dimensionless numbers from raw physical inputs. Two of
them, $Ro$ and $Bu$, do almost all the work of predicting which of this
book's reduced models applies at a given scale. This chapter turns that
observation into a map: place any flow on the $(Ro,Bu)$ plane and read off
directly whether it's unbalanced, baroclinically active, or
columnar/barotropic — and see exactly where every earlier chapter's example
flow actually lands.

## The model

Notebook: `notebooks/ch04_scaling.py` → exported to `/nb/ch04_scaling/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}).

$$Ro=\frac{U}{f_0L}\ \gtrless\ 1,\qquad
  Bu=\left(\frac{L_R}{L}\right)^2\ \gtrless\ 1,\qquad L_R=\frac{NH}{f_0}.$$

`gfdlib.scaling.classify_regime` labels every $(Ro,Bu)$ point unbalanced /
QG-baroclinic / QG-barotropic; the notebook renders it as a filled region
map with five reference points from earlier chapters' example flows placed
directly on it.

{{< marimo src="/nb/ch04_scaling/" >}}

## Both fluids

- **Atmosphere:** synoptic weather sits near $Bu\sim1$ — not a coincidence,
  since that is exactly the scale at which baroclinic instability (Ch. 16)
  most efficiently converts available potential energy into storms.
- **Ocean:** the wind-driven gyre interior (Ch. 20) sits deep in $Ro\ll1$,
  but its western boundary current is controlled by friction, a physics
  this $(Ro,Bu)$-only map cannot see — a useful reminder that even a good
  map has a domain of validity.

## Exercises

1. *(analytic)* Show algebraically that $Bu=(Fr/Ro)^2$ (both are ratios
   built from the same $U,L,H,N,f_0$), so the $(Ro,Bu)$ map is equivalent
   to an $(Ro,Fr)$ map — confirm this against the definitions in
   `gfdlib.scaling`.
2. *(computational)* Modify `ch04_scaling.py` to add a sixth reference
   point for Ch. 12's internal-wave beams, using that chapter's typical
   $\hat N,\hat f,\hat\omega$ to estimate an effective $Ro,Bu$ (this
   requires some judgment — internal waves aren't really a "balanced
   flow" scale in the same sense as the other five points, which is itself
   worth noting in your answer).
3. *(exploratory)* Sweep the marker along the $Bu=1$ line from very small
   to very large $Ro$. At what $Ro$ does the marker cross into "unbalanced"?
   Does that crossing point depend on where along the $Bu=1$ line you
   started? (It shouldn't — re-read `classify_regime`'s logic to see why
   the $Ro$ boundary is checked independently of $Bu$.)

## Further reading

Vallis, *AOFD* 2nd ed., §2.7-2.8 (regime diagrams); Pedlosky, *Geophysical
Fluid Dynamics* 2nd ed., §6.1-6.3; Cushman-Roisin & Beckers, *Introduction
to Geophysical Fluid Dynamics*, ch. 4-5.
