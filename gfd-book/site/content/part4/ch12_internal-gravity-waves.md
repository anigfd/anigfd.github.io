---
title: "Ch 12. Internal gravity waves"
weight: 412
part: 4
---

## Overview

Stratified fluids support a wave family with a genuinely strange dispersion
relation: frequency depends only on the *angle* of the wavevector, not its
magnitude. Shake a stratified tank at one point and the response isn't a
circular ripple — it's four narrow beams radiating out at a fixed angle,
carrying energy in a direction **perpendicular** to their own phase
propagation. This chapter builds that dispersion relation from the linear
equations of motion, then uses it to explain two things every internal wave
does: it picks its own propagation angle from $\omega$, $N$, and $f$ alone,
and it **reflects** wherever the local stratification $N(z)$ drops below its
own frequency — a mechanism identical in spirit to a light ray hitting total
internal reflection, except here the "refractive index" is the local
buoyancy frequency.

## The model

Notebook: `notebooks/ch12_internal-gravity-waves.py` → exported to `/nb/ch12_internal-gravity-waves.html`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$\frac{\partial^2 q}{\partial t^2} = -\Big(N^2(z)\,\psi_{xx}+f^2\,\psi_{zz}\Big)+S,
  \qquad q=\nabla^2\psi,\qquad u=\psi_z,\ w=-\psi_x.$$

The $(x,z)$ plane reuses `gfdlib.spectral.Grid` exactly as the horizontal
chapters do, with the grid's $y$-axis standing in for depth. A leapfrog
scheme (`gfdlib.internalwaves.step_leapfrog`) evaluates $\psi_{xx},\psi_{zz}$
spectrally but multiplies by $N^2(z)$ in physical space, so — unlike the
textbook two-layer or constant-$N$ treatments — the stratification is free to
vary with depth without any special-case machinery. Two presets: uniform $N$
(the classic St. Andrew's Cross) and a cosine-shaped thermocline that creates
a **wave duct**, bounded above and below by genuine turning-level reflections
the solver produces on its own.

## What to look for

1. **Four beams, one angle.** Run the *St. Andrew's Cross* preset: the
   point forcing radiates four beams whose inclination is set by
   $\omega/N$ alone — steepen them toward vertical by raising
   $\hat\omega$ toward $\hat N$, flatten them by lowering it. No property
   of the *forcing* (its size, its strength) moves that angle; only the
   frequency does.
2. **Phase across, energy along.** Zoom your attention onto one beam and
   scrub consecutive frames: the crests march *across* the beam while the
   beam itself extends *along* its own axis. That is
   $\mathbf c_g\perp\mathbf k$ seen directly — energy propagating
   perpendicular to phase, the signature strangeness of internal waves.
3. **The waveguide.** Switch to the *thermocline duct* preset: where the
   local $N(z)$ falls below $\omega$, the beam bends back — a turning
   level. The solver was never told about reflection; it emerges from
   $N^2(z)\psi_{xx}$ alone. Watch the reflected beam interfere with the
   incident one inside the duct.
4. **Rotation closes the window from below.** Turn $\hat f$ up from zero:
   nothing dramatic happens until $\hat\omega$ approaches $\hat f$, and
   then the beams flatten toward horizontal and the response stalls —
   the other edge of the internal-wave frequency window $f<\omega<N$.
   Waves exist only inside it; both presets live or die by where
   $\hat\omega$ sits in that window.

{{< marimo src="/nb/ch12_internal-gravity-waves.html" >}}

## Both fluids

- **Atmosphere:** mountain waves are internal gravity waves forced by flow
  over topography; they propagate obliquely upward at the angle this
  chapter's dispersion relation predicts, and can reflect at a
  critical/turning level where a wind reversal or a stratification change
  intervenes — the same geometry as this chapter's duct, with a different
  cause.
- **Ocean:** the barotropic tide flowing over ocean-ridge topography
  generates internal tides — low-mode internal waves radiating away at a
  fixed angle set by the tidal frequency and the local $N(z)$ profile; their
  breaking supplies roughly 90% of the energy for deep-ocean diapycnal
  mixing.

## Exercises

1. *(analytic)* From $\omega^2=N^2\cos^2\theta+f^2\sin^2\theta$, show that
   $\omega$ is bounded between $f$ and $N$ for any real $\theta$, and that
   group velocity is perpendicular to phase velocity (differentiate
   $\omega(k,m)$ and dot the result into $\mathbf k=(k,m)$).
2. *(computational)* Modify `ch12_internal-gravity-waves.py` to force with
   two point sources at different depths instead of one. Do the resulting
   beam patterns superpose linearly, as the governing equation's linearity
   predicts?
3. *(exploratory)* In the thermocline-duct preset, sweep $\hat\omega$ upward
   at fixed $\hat N_{mean},\hat N_{amp}$: find the value at which the two
   turning levels merge with the domain boundary and the wave stops being
   trapped.

## Further reading

Vallis, *AOFD* 2nd ed., §6.1–6.4 (internal waves) and §6.11 (internal tides);
Gill, *Atmosphere-Ocean Dynamics*, ch. 6; Sutherland, *Internal Gravity
Waves* (2010); Garrett & Munk (1972, 1979).
