---
title: "Ch 1. Kinematics & the material derivative"
weight: 101
part: 1
---

## Overview

Drop a patch of dye into a fluid near a vortex embedded in a straining
current: some of it survives as a coherent blob orbiting the vortex core,
the rest is stretched into thin filaments and swept away. This opening
chapter builds the two ideas every later chapter leans on without saying
so — the **material derivative** ($D/Dt=\partial_t+\mathbf u\cdot\nabla$,
the rate of change *following a parcel*) and the local **strain-vorticity**
decomposition of a velocity field — using nothing but a prescribed 2D flow
and particle trajectories.

## The model

Notebook: `notebooks/ch01_kinematics.py` → exported to `/nb/ch01_kinematics/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}).

$$\frac{d\mathbf x}{dt}=\mathbf u(\mathbf x,t),\qquad
  W=S_n^2+S_s^2-\zeta^2,\qquad
  S_n=u_x-v_y,\ \ S_s=v_x+u_y,\ \ \zeta=v_x-u_y.$$

`gfdlib.kinematics.velocity_field` evaluates a steady superposition of
uniform strain and a smooth (Lamb-Oseen-like) vortex in closed form;
`gfdlib.kinematics.okubo_weiss` differentiates it with central differences
to get the Okubo-Weiss parameter $W$ — the standard diagnostic (Okubo 1970,
Weiss 1991) for separating vorticity-dominated cores ($W<0$) from
strain-dominated, filament-producing regions ($W>0$). A ring of tracer
particles is released around the vortex and advected with
`gfdlib.timestep.rk4`.

{{< marimo src="/nb/ch01_kinematics/" >}}

## Both fluids

- **Atmosphere:** the Okubo-Weiss decomposition is the standard method for
  automatically detecting and tracking coherent vortices (cyclones,
  anticyclones, polar vortex fragments) in reanalysis wind fields — a
  storm's core sits in $W<0$ territory, its surrounding frontal shear zones
  in $W>0$.
- **Ocean:** mesoscale eddy-detection algorithms applied to satellite
  altimetry are built almost entirely on this same $W<0$ criterion; the
  filamentation you watched here is exactly how an eddy's edge sheds
  material into the surrounding strain field, the first step in eddy decay.

## Exercises

1. *(analytic)* For pure strain alone ($\Gamma=0$), show by direct
   differentiation of $u=-\alpha x,\ v=\alpha y$ that $S_n=-2\alpha$,
   $S_s=0$, $\zeta=0$, hence $W=4\alpha^2$ everywhere — confirm this matches
   `gfdlib.kinematics.okubo_weiss` at $\Gamma=0$ (see `tests/test_gfdlib.py`).
2. *(computational)* Modify `ch01_kinematics.py` to release the tracer ring
   at several different radii at once (e.g. $0.5\sigma$, $1.3\sigma$,
   $2.5\sigma$) and color each ring differently. Does every particle in the
   innermost ring stay coherent regardless of $\alpha$, for the range of
   $\alpha$ the sliders allow?
3. *(exploratory)* Find, by trial and error with the sliders, the
   approximate ratio $\Gamma/(\alpha\sigma^2)$ at which the $W=0$ contour
   just barely encloses a nonzero region. (This is the vortex-survival
   threshold; compare your empirical value to the literature on vortex
   stripping by ambient strain, e.g. Legras & Dritschel 1993.)

## Further reading

Okubo, A. (1970), *Deep-Sea Res.* **17**, 445–454; Weiss, J. (1991),
*Physica D* **48**, 273–294; Legras, B. & Dritschel, D. G. (1993), *Fluid
Dyn. Res.* **10**, 159–176 ("Vortex stripping and the generation of high
vorticity gradients in two-dimensional flows"); Vallis, *AOFD* 2nd ed.,
§2.1-2.3; Salmon, *Lectures on Geophysical Fluid Dynamics*, ch. 1.
