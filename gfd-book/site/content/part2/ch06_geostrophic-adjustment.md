---
title: "Ch 6. Rotating shallow water & geostrophic adjustment"
weight: 206
part: 2
---

## Overview

Perturb a rotating fluid away from rest and two things happen at once: fast
**inertia-gravity waves** radiate the imbalance outward, and rotation traps a
piece of the disturbance behind as a standing, balanced flow. This is
**geostrophic adjustment**, and it is the mechanism by which almost every
localized disturbance in the atmosphere or ocean — a convective burst, a
wind-stress event, an eddy shed from an unstable current — settles into the
slowly evolving balanced motion that the rest of this book is about. The
**deformation radius** $L_R=c/f_0$ is the scale that decides the outcome:
disturbances much smaller than $L_R$ disperse almost entirely as gravity
waves, while disturbances much larger than $L_R$ survive almost unchanged as
geostrophic flow. What is conserved throughout — and what therefore *predicts*
the balanced end state before you integrate a single time step — is the
linear **potential vorticity** $q=\zeta-\hat\eta$.

## The model

Notebook: `notebooks/ch06_geostrophic-adjustment.py` → exported to `/nb/ch06_geostrophic-adjustment.html`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$\hat\eta_t+\hat u_x+\hat v_y=0,\qquad
  \hat u_t-\hat f\hat v=-\hat\eta_x,\qquad
  \hat v_t+\hat f\hat u=-\hat\eta_y,
  \qquad \hat f=1+\hat\beta\,(y-y_0).$$

Pseudo-spectral on a doubly-periodic domain (`gfdlib.spectral.Grid`),
classical RK4 (`gfdlib.timestep.rk4` — the system is exactly linear, so no
integrating factor is needed), with the new `gfdlib.shallowwater` module
supplying the RHS and two potential-vorticity primitives:
`potential_vorticity` (diagnose $q=\zeta-\hat\eta$) and `invert_pv` (solve
$(\nabla^2-1)\eta_{bal}=q$ for the balanced height field carrying a given
PV). Three presets: a small-scale burst ($\sigma\ll L_R$), a large-scale surge
($\sigma\gg L_R$), and a $\beta$-plane run that shows the balanced remnant
drift westward as a Rossby wave.

## What to look for

1. **Waves win at small scales.** Run the *small-scale burst* preset and
   scrub the frame slider from the start: rings of inertia-gravity waves
   carry essentially the whole height anomaly away (at the nondimensional
   wave speed $c=1$ — check the ring radius against elapsed time), and the
   final $\eta$ is nearly flat. The overlay of $\eta_{bal}=$
   `invert_pv`$(q_0)$ predicted that end state *before the run started* —
   compare them at the last frame.
2. **Balance wins at large scales.** Switch to the *large-scale surge*
   preset: now most of the height anomaly survives, and what changes
   instead is the wind field, which spins up into a rim current around the
   anomaly. Same equations, same initial shape — only $\sigma/L_R$
   differs.
3. **Energetics of the remnant.** Watch the energy panel in the
   large-scale case: even when balance "wins," at most one third of the
   released potential energy ends up in the balanced flow — the rest
   radiates. Rossby's classic result, visible as the gap between the
   curves.
4. **The remnant is not stuck.** In the *β-plane* preset the balanced
   remnant itself drifts westward — adjustment hands the disturbance off
   to Rossby-wave dynamics (Ch. 8–9), which is where its story continues.

{{< marimo src="/nb/ch06_geostrophic-adjustment.html" >}}

## Both fluids

- **Atmosphere:** a convective outflow or a jet-stream meander is released
  unbalanced; within a day it sheds inertia-gravity waves and leaves behind a
  geostrophically balanced anomaly — the everyday mechanism that keeps
  synoptic weather systems in near-balance despite continuous forcing.
- **Ocean:** wind-driven upwelling or a detrainment event from an unstable
  boundary current perturbs the thermocline; geostrophic adjustment converts
  part of that perturbation directly into a mesoscale eddy, while the rest
  escapes as fast barotropic and internal gravity waves.

## Exercises

1. *(analytic)* Starting from $q=\zeta-\hat\eta$ and the momentum/continuity
   equations, derive $q_t=-\hat\beta\hat v$. Confirm it reduces to exact
   pointwise conservation on the $f$-plane.
2. *(computational)* Modify `ch06_geostrophic-adjustment.py` to add linear
   (Rayleigh) drag $-r\hat u,-r\hat v$ to the momentum equations. Does $q$
   still evolve as $q_t=-\hat\beta\hat v$, or does drag act as an additional
   PV source? Check against the modified notebook's own diagnostic.
3. *(exploratory)* Sweep $\sigma/L_R$ finely between the "small" and "large"
   presets and plot the fraction of $\eta_0$'s peak amplitude that survives in
   $\eta_{bal}=$ `invert_pv(q0)`. At what $\sigma/L_R$ does the survival
   fraction cross 50%?

## Further reading

Vallis, *AOFD* 2nd ed., §3.7–3.9 (geostrophic adjustment) and §5.2 (linear PV);
Gill, *Atmosphere-Ocean Dynamics*, ch. 7; Salmon, *Lectures on GFD*, §2.4;
Rossby (1938), *J. Mar. Res.* **1**, 239–263.
