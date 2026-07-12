---
title: "Ch 18. Geostrophic (2D) turbulence"
weight: 618
part: 6
---

## Overview

Rotation and stratification make large-scale atmospheric and oceanic flow
approximately two-dimensional — and 2D turbulence is a *different animal* from
the 3D turbulence of engineering. Because vortex stretching is absent, the flow
conserves enstrophy $Z=\tfrac12\langle\zeta^2\rangle$ alongside energy
$E=\tfrac12\langle|\nabla\psi|^2\rangle$, and that second invariant forces
energy **upscale**: small vortices merge into ever-larger ones (the *inverse
cascade*), while enstrophy drains to small scales through filaments. On a
rotating planet the upscale march does not continue forever — when eddies grow
big enough to feel $\beta$, Rossby-wave dynamics arrest the cascade near the
Rhines scale $k_\beta\simeq\sqrt{\beta/2U}$ and reorganize the flow into
**zonal jets**. This chapter's notebook lets you watch both acts: the dual
cascade, and its arrest.

## The model

Notebook: `notebooks/ch18_geostrophic-turbulence.py` → exported to `/nb/ch18_geostrophic-turbulence/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$\frac{\partial \zeta}{\partial t} + J(\psi,\zeta) + \beta v
  = -\nu(-\nabla^2)^{n_\nu}\zeta, \qquad \nabla^2\psi = \zeta .$$

Pseudo-spectral, 2/3-rule dealiased, with an integrating factor that treats
dissipation *and* the Rossby-wave propagator exactly
(`gfdlib.timestep.ifrk4_step`). Presets: decaying McWilliams (1984)
turbulence, a hyperviscous variant, and a $\beta$-plane jet run.

## What to look for

1. **Merger all the way up.** Run the *decaying (McWilliams 1984)* preset
   and scrub: a fine-grained random field organizes itself into a
   handful of large, long-lived vortices through repeated mergers —
   Ch. 7's two-vortex merger applied recursively. This is the inverse
   cascade in physical space.
2. **The two invariants part ways.** Watch the energy and enstrophy
   curves together: energy stays nearly flat while enstrophy collapses.
   That asymmetry is the entire logic of 2D turbulence — dissipation
   (which acts at small scales) can reach the enstrophy but barely
   touches the energy, so the energy has nowhere to go but up-scale.
   Fjørtoft's theorem (exercise 1) is this plot, proved.
3. **The spectrum migrates.** In the spectral panel, follow the energy
   peak: it marches steadily toward lower $k$ while a power-law tail
   feeds enstrophy toward the dissipation range. Compare the tail's slope
   with the $k^{-3}$ enstrophy-cascade prediction — and note where (and
   why) coherent vortices make it steeper.
4. **Jets from a knob.** Switch to the *β-plane* preset: the isotropic
   inverse cascade proceeds until eddies reach the Rhines scale, then
   anisotropizes into east-west bands — watch the Hovmöller panel of
   $\bar u(y,t)$ develop persistent stripes. Sweep $\beta$ and check the
   jet spacing against $\pi\sqrt{2U/\beta}$ (exercise 3 makes this
   quantitative). Zero knobs were labeled "make jets"; $\beta$ plus an
   arrested cascade is sufficient.

{{< marimo src="/nb/ch18_geostrophic-turbulence/" >}}

## Both fluids

- **Atmosphere:** the banded winds of Jupiter and Saturn are the textbook
  Rhines-arrested inverse cascade; on Earth, the eddy-driven midlatitude jet is
  maintained by exactly this upscale momentum transfer from baroclinic eddies
  (Ch. 16 supplies the eddies).
- **Ocean:** the Southern Ocean and the subtropical gyres carry multiple
  quasi-zonal jets ("striations") with spacing near the local Rhines scale, and
  the mesoscale eddy field ($L\sim100$ km) inverse-cascades energy toward the
  basin scale until $\beta$ — and bottom drag — intervene.

## Exercises

1. *(analytic)* From conservation of $E$ and $Z$, prove Fjørtoft's theorem: if
   spectral energy spreads from a middle wavenumber $k_1$ to $k_0<k_1<k_2$ with
   $k_2/k_1 = k_1/k_0 = 2$, more energy must go to $k_0$ than to $k_2$. Where
   does the enstrophy go?
2. *(computational)* Modify `ch18_geostrophic-turbulence.py` to add linear
   (Ekman) drag $-r\zeta$ to the linear operator. How does the final vortex
   size — and the jet amplitude in the $\beta$ preset — depend on $r$?
3. *(exploratory)* Sweep $\beta$ at fixed energy: find the smallest $\beta$ at
   which the Hovmöller diagram of $\bar u(y,t)$ shows persistent stripes, and
   check the jet spacing against $\pi/k_\beta$ with $U=\sqrt{2E}$.

## Further reading

Vallis, *AOFD* 2nd ed., §§11.1–11.4 (two-dimensional turbulence) and ch. 12
(geostrophic turbulence and jets); Salmon, *Lectures on GFD*, ch. 4;
Kraichnan (1967); Rhines (1975); McWilliams (1984).
