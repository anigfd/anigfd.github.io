---
title: "Ch 7. Vorticity and potential vorticity"
weight: 207
part: 2
---

## Overview

Potential vorticity is the single most useful diagnostic quantity in this
subject because of one property: **invertibility**. Knowing $q(x,y)$
everywhere — plus a balance condition and boundary conditions — fixes the
entire velocity field, uniquely, with no other information needed. And
because $q$ is also materially conserved, that same fact tells you how the
flow evolves: invert for the velocity, advect $q$ with it, invert again. This
chapter builds that two-step algorithm from the vorticity equation and uses
it to explain two classic results: why two vortices merge or orbit forever
depending on their separation, and why regions where eddies have mixed PV
into a **staircase** of flat, well-mixed bands develop sharp **jets** exactly
at the risers between bands — not within the quiet plateaus.

## The model

Notebook: `notebooks/ch07_vorticity-pv.py` → exported to `/nb/ch07_vorticity-pv.html`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$q_t+J(\psi,q)=0,\qquad \zeta=\nabla^2\psi=q-\beta y.$$

Reuses the exact pseudo-spectral machinery from Ch. 18
(`gfdlib.spectral.Grid`, `gfdlib.timestep.ifrk4_step`) — this chapter's new
contribution is `gfdlib.pv`, which builds initial PV fields: a periodic-safe
Gaussian blob for placing vortices, and an idealized alternating-band
staircase. A fully reactive "invert, don't integrate" panel shows the balance
step alone (no time-stepping); a second, Run-gated panel evolves the same
field forward to show the resulting dynamics.

## What to look for

1. **Inversion is instant — use it first.** The upper panel needs no Run
   button: it takes the PV field you built and solves for $\psi$ and the
   winds. Drag the separation slider and watch the velocity field between
   two like-signed vortices reorganize — everything you see is *diagnosed*
   from $q$ alone, which is the invertibility principle doing its work.
2. **Find the merger threshold.** In the *merger* preset, run the
   evolution at a few separations: well-separated vortices orbit each
   other essentially forever, while below a critical separation (a few
   core radii) they wrap around each other and merge into one core,
   throwing off filament arms. The transition is surprisingly sharp — this
   is the same vortex-merger physics that drives Ch. 18's inverse cascade.
3. **The dipole travels.** The *opposite-signed pair* preset produces the
   one configuration that self-propels: each vortex advects the other in
   the same direction. Confirm the speed falls as you increase the
   separation.
4. **Jets live at the risers.** In the *staircase* preset, look at the
   zonal-mean wind: it peaks exactly at the sharp PV jumps between the
   flat, well-mixed bands, not inside the plateaus. Then use the reactive
   panel (no Run needed) to vary the number of bands and watch exercise 3's
   question take shape.

{{< marimo src="/nb/ch07_vorticity-pv.html" >}}

## Both fluids

- **Atmosphere:** the tropopause's PV field is routinely mapped and
  inverted operationally — "PV thinking" lets forecasters read cyclogenesis
  directly off a PV map, and jet streams are themselves regions of sharp PV
  gradient exactly as in this chapter's staircase.
- **Ocean:** mesoscale eddies are, to a good approximation, isolated PV
  anomalies that interact by mutual advection just like this notebook's
  vortex pairs; multiple oceanic jets (e.g. in the Antarctic Circumpolar
  Current) are believed to be maintained by the same PV-staircase mechanism
  demonstrated here.

## Exercises

1. *(analytic)* Show that $Dq/Dt=0$ combined with $\zeta=\nabla^2\psi$ is a
   *closed* system — that is, that no other equation or unknown is needed to
   advance $q$ forward in time. This closure is what "invertibility" buys
   you.
2. *(computational)* Modify `ch07_vorticity-pv.py`'s merger preset to use
   three like-signed vortices arranged in a triangle instead of two. Do they
   merge pairwise, or all at once?
3. *(exploratory)* In the staircase preset, use the "invert, don't
   integrate" panel (no need to press Run) to find how the jet strength
   scales with the number of bands at fixed amplitude — does adding more,
   narrower bands make each jet stronger or weaker?

## Further reading

Vallis, *AOFD* 2nd ed., §4.1–4.6 (vorticity and PV) and §4.8 (invertibility);
Hoskins, McIntyre & Robertson (1985), *Quart. J. Roy. Meteor. Soc.* **111**,
877–946 ("On the use and significance of isentropic potential vorticity
maps"); Dritschel & McIntyre (2008), *J. Atmos. Sci.* **65**, 855–874
(PV staircases and multiple jets).
