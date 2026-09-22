---
title: "Ch 21. Buoyancy-driven / overturning circulation"
weight: 621
part: 6
---

## Overview

The Hadley cell and the ocean's meridional overturning circulation (MOC)
look nothing alike on a map, but they share a single mechanism: heat a
fluid unevenly along one boundary and it organizes into one overturning
cell — rising where it's heated, sinking where it's cooled — not many
small convective rolls. This chapter closes Part VI (and the book's
physically-motivated main path) by building that mechanism directly, and
by deriving Munk's (1966) classic "abyssal recipe" — the simple balance
that first let oceanographers infer the deep ocean's overturning rate from
its observed stratification.

## The model

Notebook: `notebooks/ch21_overturning.py` → exported to
`/nb/ch21_overturning.html`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}):

$$w\frac{\partial T}{\partial z}=\kappa\frac{\partial^2T}{\partial z^2},
  \qquad
  \frac{\partial\theta}{\partial t}+J(\psi,\theta)=\nabla^2\theta+Q(x,z).$$

`gfdlib.overturning.abyssal_profile` solves Munk's balance exactly (a
linear 2nd-order ODE, general solution $A+Be^{wz/\kappa}$, boundary
constants from a $2\times2$ `numpy.linalg.solve`) — large $w/\kappa$
compresses the transition from bottom to surface values into a thin layer
near the top, the real ocean's near-uniform abyssal temperature capped by
a thin thermocline. `gfdlib.overturning.rhs_overturning` reuses ch. 14's
Boussinesq vorticity-streamfunction machinery (`convection.ChannelGrid`)
unchanged, replacing uniform bottom heating with a differential,
surface-concentrated heating field — Rossby's (1965) "horizontal
convection" mechanism — and was verified to produce exactly the expected
circulation sense: rising motion under the heated column, sinking under
the cooled one, in every configuration tested.

## What to look for

1. **One ratio controls the profile.** In Part A, only $w/\kappa$
   matters (verify: doubling both leaves the curve unchanged). Crank the
   ratio up and the exponential's e-folding depth $\kappa/w$ shrinks:
   a nearly uniform abyss capped by a thin thermocline — the observed
   shape of essentially every mid-ocean temperature profile. Munk's
   famous move was to run this logic backwards: from the *observed*
   thermocline depth and an estimate of $w$, infer the ocean's interior
   mixing rate $\kappa\sim10^{-4}\,$m$^2\,$s$^{-1}$.
2. **One cell, not many.** Run Part B: differential heating along one
   boundary produces a single domain-filling overturning cell — rising
   over the heated end, sinking over the cooled end — rather than Ch. 14's
   array of counter-rotating rolls. The forcing's own asymmetry organizes
   the flow at the largest available scale.
3. **No threshold this time.** Sweep $Ra$ downward: unlike Ch. 14, the
   circulation never switches off — it just weakens smoothly. Horizontal
   convection has no conduction-only rest state to bifurcate from,
   because a horizontal temperature gradient along a boundary *cannot*
   be balanced by pure conduction (exercise 3 makes you say why).
4. **Hunt the asymmetry.** Horizontal-convection theory predicts the two
   branches of the cell are not mirror images: the surface-*cooled* end is
   convectively destabilized and should concentrate into a narrower,
   faster sinking branch as you raise $Ra$, while the heated end stays
   broad and gentle. Test it across the $Ra$ slider's range. The real
   counterparts are extreme versions of this fingerprint — the MOC's few
   localized deep-water-formation sites versus basin-wide diffuse
   upwelling.

{{< marimo src="/nb/ch21_overturning.html" >}}

## Both fluids

- **Atmosphere:** the Hadley cell is the thermally-direct overturning
  circulation between the heated tropics and the cooler subtropics —
  rising near the equator, sinking near 30° latitude, with the trade
  winds and subtropical highs as its surface signature. Held & Hou's
  (1980) angular-momentum-conserving theory extends this same "single
  overturning cell" logic with rotation added.
- **Ocean:** the meridional overturning circulation sinks dense water at
  high latitudes (North Atlantic Deep Water, Antarctic Bottom Water),
  spreads it through the deep ocean, and slowly upwells it back toward
  the surface elsewhere — set by exactly the diffusion/upwelling balance
  in Part A, and driven by exactly the differential (here: buoyancy-,
  not just heat-) forcing mechanism in Part B.

## Exercises

1. *(analytic)* Solve $w\,dT/dz=\kappa\,d^2T/dz^2$ from scratch (general
   solution $T=A+Be^{wz/\kappa}$, then the two boundary conditions) and
   confirm it matches `gfdlib.overturning.abyssal_profile`.
2. *(computational)* Modify `ch21_overturning.py`'s Part B to use a
   heating pattern with TWO full cycles across the domain
   ($\cos(4\pi x/L_x)$ instead of $\cos(2\pi x/L_x)$). Does the flow
   organize into four cells, or does it prefer to merge into fewer,
   larger ones — and if so, why might that connect to the inverse
   cascade of Ch. 18?
3. *(exploratory)* In Part B, hold $Q_0$ fixed and sweep $Ra$ from low to
   high. Is there a sharp threshold before any circulation appears (like
   Ch. 14's $Ra_c$), or does the overturning strength grow smoothly from
   zero? What does that tell you about how "horizontal convection" differs
   from Rayleigh-Bénard convection at a basic level?

## Further reading

Munk, W. H. (1966), *Deep-Sea Res.* **13**, 707–730 ("Abyssal recipes");
Rossby, H. T. (1965), *Deep-Sea Res.* **12**, 9–16 ("On thermal
convection driven by non-uniform heating from below: an experimental
study"); Held, I. M. & Hou, A. Y. (1980), *J. Atmos. Sci.* **37**,
515–533 ("Nonlinear axially symmetric circulations in a nearly inviscid
atmosphere"); Vallis, *AOFD* 2nd ed., §14.1-14.2 (Hadley cell), §15.1-15.3
(thermohaline circulation); Pedlosky, *Geophysical Fluid Dynamics* 2nd
ed., §6.1-6.5.
