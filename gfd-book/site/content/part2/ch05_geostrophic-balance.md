---
title: "Ch 5. Geostrophic & hydrostatic balance"
weight: 205
part: 2
---

## Overview

Two balances — geostrophic (Coriolis against pressure gradient) and
hydrostatic (gravity against vertical pressure gradient) — together with the
chain rule produce one of the most useful diagnostic relations in the
subject: **thermal wind**. It says the vertical shear of the geostrophic
wind is fixed entirely by the horizontal temperature gradient, with no need
to know the wind itself. This is how the jet stream can be *predicted* from
a weather map showing only temperature, and it's why the jet sits where it
does: exactly where the temperature gradient reverses sign, at the
tropopause.

## The model

Notebook: `notebooks/ch05_geostrophic-balance.py` → exported to `/nb/ch05_geostrophic-balance/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$u_g=-\frac{\phi_y}{f},\qquad v_g=\frac{\phi_x}{f},\qquad \phi_z=b
  \quad\Longrightarrow\quad fu_{g,z}=-b_y.$$

Unlike every other notebook in this book, there is no PDE and no time
integration here — `gfdlib.balance` provides closed-form expressions for an
idealized frontal buoyancy field and the geostrophic wind it implies, exact
to floating-point precision. Move a slider and the reconstruction updates
live.

{{< marimo src="/nb/ch05_geostrophic-balance/" >}}

## Both fluids

- **Atmosphere:** this is literally how forecasters read a thickness chart —
  tightly packed isotherms on a temperature map mean strong vertical wind
  shear is guaranteed nearby, jet stream included, before a single wind
  observation is consulted.
- **Ocean:** the same relation underlies the **dynamic method** used for a
  century to estimate ocean currents from hydrographic (temperature/salinity)
  sections alone — geostrophic currents relative to a reference level are
  computed directly from the density field, exactly as this notebook
  reconstructs $u_g$ from $b$.

## Exercises

1. *(analytic)* Derive $fu_{g,z}=-b_y$ from geostrophic and hydrostatic
   balance, starting from the horizontal momentum and hydrostatic equations.
   Show that at low Rossby number the acceleration terms are the ones that
   drop out, not the Coriolis or pressure-gradient terms.
2. *(computational)* Modify `ch05_geostrophic-balance.py` to use a
   *double*-front buoyancy field (two `tanh` fronts of opposite sign,
   representing a polar front and a subtropical front). How many jet cores
   appear, and where?
3. *(exploratory)* At what value of $Ro=U_{jet}/(fL_y)$ (read from the
   jet-core plot's title) would you say geostrophic balance has stopped
   being a good approximation? Push $f$ down until you're not sure anymore.

## Further reading

Vallis, *AOFD* 2nd ed., §2.7–2.8 (geostrophic and hydrostatic balance) and
§2.9 (thermal wind); Holton & Hakim, *An Introduction to Dynamic
Meteorology*, §3.4; Gill, *Atmosphere-Ocean Dynamics*, §7.10 (the dynamic
method).
