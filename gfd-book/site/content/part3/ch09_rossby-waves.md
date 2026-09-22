---
title: "Ch 9. Rossby waves"
weight: 309
part: 3
---

## Overview

A Rossby wave's group velocity depends on its own wavenumber, so a wave
packet moving through a slowly varying background doesn't travel in a
straight line — its wavenumber refracts as it goes, and it follows a curved
*ray*, exactly as light refracts through a lens. This chapter builds that
ray theory from the dispersion relation and group velocity you've already
used (ch. 18, ch. 8), then takes it somewhere the rest of this book hasn't
gone: onto the sphere itself. With the right background flow, the rays this
notebook traces are **exact great circles** — a genuine theorem (Hoskins &
Karoly 1981), not a plotting trick, and the mechanism behind how real
atmospheric teleconnections (like the Pacific–North America pattern excited
by El Niño) arc across the globe.

## The model

Notebook: `notebooks/ch09_rossby-waves.py` → exported to `/nb/ch09_rossby-waves.html`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$\frac{d\mathbf x}{dt}=\frac{\partial\omega}{\partial\mathbf k},\qquad
  \frac{d\mathbf k}{dt}=-\frac{\partial\omega}{\partial\mathbf x}.$$

`gfdlib.rossby` deliberately avoids a hand-derived symbolic ray equation:
`ray_rhs` differentiates the single scalar `dispersion_omega` by central
differences, so correctness rests on one formula rather than a chain of
algebra, and `great_circle_deviation` checks the result against an
independent, coordinate-free geometric criterion — no PDE, no time-stepping
in the usual sense, just `gfdlib.timestep.rk4` on a 3-variable ray. Cheap
enough that the whole notebook is reactive, with no Run button.

## What to look for

1. **The arc, not the straight line.** Launch a fan of rays from
   $\phi_0=20°$: each one climbs poleward, turns at its own **turning
   latitude** (where its meridional wavenumber passes through zero), and
   arcs back toward the equator. No ray crosses its turning latitude —
   refraction by the planetary vorticity gradient bends it back, exactly
   as a lens bends light back toward the dense medium.
2. **Great circles, verified.** The `great_circle_deviation` readout stays
   at machine precision for solid-body rotation — the Hoskins–Karoly
   theorem holding *exactly*, not approximately, in front of you. (Then
   see exercise 2 for how quickly it breaks when the background flow isn't
   solid-body.)
3. **Rotation strength sets the reach.** Increase $\Omega_s$: the
   stationary wavenumber $K_s$ falls everywhere (stronger westerlies can
   hold longer waves stationary), and since a ray turns where $K_s(\phi)$
   has dropped to its own zonal wavenumber, lowering the whole $K_s$
   profile slides each turning latitude *equatorward*. Confirm it in the
   fan: the same rays that grazed the pole at small $\Omega_s$ turn low
   at large $\Omega_s$.
4. **Where rays bunch, anomalies live.** Notice the rays launched at
   slightly different angles converge near their turning latitudes (the
   ray density is highest there). Ray theory says wave amplitude
   accumulates where rays crowd — one reason observed teleconnection
   centers of action sit at preferred latitudes rather than smearing along
   the whole path.

{{< marimo src="/nb/ch09_rossby-waves.html" >}}

## Both fluids

- **Atmosphere:** stationary Rossby wave trains excited by tropical
  convection anomalies propagate poleward along great-circle-like arcs
  before curving back equatorward — the observed teleconnection patterns
  linking ENSO to midlatitude weather are, to good approximation, segments
  of exactly the rays this notebook traces.
- **Ocean:** oceanic Rossby waves (much slower, non-dispersive westward
  phase propagation dominates over the ray-refraction effects shown here)
  carry the memory of wind-stress anomalies across entire ocean basins,
  setting the timescale for basin-scale adjustment.

## Exercises

1. *(analytic)* From $\omega=-\beta k/(k^2+l^2)$, derive
   $\mathbf c_g=\beta(k^2-l^2,2kl)/(k^2+l^2)^2$ and show
   $\mathbf c_g\cdot\mathbf k>0$ always — group velocity always has an
   eastward component relative to the phase, even though phase itself
   always moves westward.
2. *(computational)* Modify `ch09_rossby-waves.py` to use
   $U(\phi)=U_0\cos^2\phi$ instead of solid-body rotation. Confirm with the
   great-circle deviation diagnostic that the rays are no longer exact great
   circles, and describe qualitatively how they differ.
3. *(exploratory)* Find the launch angle (at fixed $\phi_0,\Omega_s$) that
   sends a ray closest to the pole without crossing it. How does that
   critical angle change as you increase $\Omega_s$?

## Further reading

Hoskins, B. J. & Karoly, D. J. (1981), *J. Atmos. Sci.* **38**, 1179–1196
("The steady linear response of a spherical atmosphere to thermal and
orographic forcing"); Vallis, *AOFD* 2nd ed., §7.1–7.4 (Rossby waves) and
§13.4 (ray tracing); Karoly, D. J. (1983), *Tellus* **35A**, 190–200.
