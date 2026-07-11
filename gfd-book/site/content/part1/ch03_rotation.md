---
title: "Ch 3. Effects of rotation"
weight: 103
part: 1
---

## Overview

Two questions about what rotation alone does to a fluid. Give a parcel a
push and remove every other force — what does it do? Not travel in a
straight line: it coasts around a circle at the **inertial period**
$T=2\pi/f$, the signature seen in nearly every filtered current-meter
record from the upper ocean. Force a rapidly-rotating, homogeneous fluid at
one level only — does the rest of the column notice? It does, rigidly: the
**Taylor-Proudman theorem** says a steady, balanced, homogeneous rotating
flow cannot vary along the rotation axis. This chapter derives both exactly,
and shows the second is not a separate calculation from Ch. 5's thermal
wind — it is that relation's zero-buoyancy limit.

## The model

Notebook: `notebooks/ch03_rotation.py` → exported to `/nb/ch03_rotation/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}).

$$\frac{du}{dt}=fv,\qquad\frac{dv}{dt}=-fu
  \qquad\Longrightarrow\qquad
  b\equiv\text{const}\ \Rightarrow\ \frac{\partial u_g}{\partial z}=\frac{\partial v_g}{\partial z}=0.$$

`gfdlib.rotation.inertial_trajectory` is the closed-form circular solution
of the first equation, cross-checked against an independent
`gfdlib.timestep.rk4` integration; Part B reuses `gfdlib.balance`'s exact
thermal-wind machinery from Ch. 5 unchanged, sweeping the buoyancy contrast
toward zero.

{{< marimo src="/nb/ch03_rotation/" >}}

## Both fluids

- **Ocean:** inertial oscillations are the dominant signal in
  storm-generated near-surface currents once tides and the mean flow are
  filtered out; Taylor columns are directly observed steering deep,
  weakly-stratified abyssal currents around seamounts and ridges far more
  than the overlying wind or buoyancy forcing would predict.
- **Atmosphere:** inertial oscillations appear as the low-level jet's
  nocturnal acceleration once daytime turbulent friction relaxes at
  sunset — a textbook case of a "sudden push" (removal of drag) exciting
  the same circular response derived here; strongly rotating, weakly
  stratified planetary atmospheres (e.g. the gas giants' deep interiors)
  show Taylor-column-organized convection.

## Exercises

1. *(analytic)* Verify by direct substitution that $w(t)=w_0e^{-ift}$
   solves $dw/dt=-ifw$ with $w=u+iv$, and that integrating once more gives
   a circle of radius $|w_0|/f$ centered at
   $z_0+w_0/(if)$ (matching `gfdlib.rotation.inertial_trajectory` and its
   test in `tests/test_gfdlib.py`).
2. *(computational)* Modify `ch03_rotation.py`'s Part B to plot the shear
   range as a function of $L_y$ instead of $\Delta b$ (holding $\Delta b$
   fixed) — does narrowing the front also drive the flow toward
   Taylor-Proudman columnarity, or only $\Delta b\to0$ does?
3. *(exploratory)* Using Part A's sliders, estimate how many inertial
   periods it would take a storm-excited ocean current (typical
   $f\sim10^{-4}\,$s$^{-1}$ at midlatitudes) to complete even a single
   orbit — convert the nondimensional period on the hodograph panel to
   hours. Does this match the day-or-so timescale commonly cited for
   inertial oscillations in oceanographic textbooks?

## Further reading

Vallis, *AOFD* 2nd ed., §2.5-2.6 (inertial oscillations),
§4.1-4.2 (Taylor-Proudman); Pedlosky, *Geophysical Fluid Dynamics* 2nd ed.,
§2.3-2.5; Cushman-Roisin & Beckers, *Introduction to Geophysical Fluid
Dynamics*, §7.1-7.3.
