---
title: "Ch 13. Wave–mean-flow interaction"
weight: 413
part: 4
---

## Overview

Two ways a linear wave leaves a permanent mark on the mean state it rides
on. A floating buoy does not end up exactly where it started once a wave
train passes, despite the wave's Eulerian-mean velocity being exactly zero
at every fixed point — it **Stokes-drifts**, a direct consequence of the
wave field being evaluated at the parcel's true, displaced position. And a
Rossby wave riding on a sheared mean flow cannot reach the latitude where
its Doppler-shifted phase speed matches the flow: its wavenumber diverges
as it approaches that **critical layer**, setting up exactly the structure
that (once dissipation or finite amplitude enters) lets it deposit
momentum into the mean flow there.

## The model

Notebook: `notebooks/ch13_wave-mean.py` → exported to `/nb/ch13_wave-mean/`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}).

$$u_S=\frac{A^2k}{2\omega},\qquad
  \omega(y,k,l)=U(y)k-\frac{\beta k}{k^2+l^2},\qquad
  \frac{dy}{dt}=\frac{2\beta kl}{(k^2+l^2)^2},\quad
  \frac{dl}{dt}=-k\frac{dU}{dy}.$$

`gfdlib.wavemean.wave_velocity`/`stokes_drift` derive and verify the Stokes
formula against direct particle advection (`gfdlib.timestep.rk4`, reusing
Ch. 1's machinery); `ray_rhs_shear` extends Ch. 9's ray tracing to a
Doppler-shifted background shear, verified to conserve $\omega$ exactly
along the ray.

## What to look for

1. **The drift is real, and it is quadratic.** In Part A, watch a particle
   ride the wave: it oscillates, but each cycle ends slightly downstream
   of where it began. Read off the measured drift and check it against
   $u_S=A^2k/2\omega$; then double the amplitude $A$ and confirm the drift
   *quadruples*. A linear wave produces a mean effect only at second
   order — this factor of four is that statement made tangible.
2. **The Eulerian mean is exactly zero.** The wave is
   $u'=A\cos(kx-\omega t)$, which averages to nothing at every fixed
   point over a period — yet the particle undeniably moves. The entire drift lives in the
   correlation between where the parcel *is* and what the wave is doing
   *there*. Lagrangian and Eulerian means are different animals.
3. **The ray runs out of road.** In Part B, launch a ray into the shear:
   $l$ grows linearly (at the rate $-k\,dU/dy$ — check the slope), the
   group velocity collapses as $1/l^3$, and the ray asymptotes to its
   critical latitude $y_c$ without ever arriving. Move $\Lambda$ and watch
   $y_c$ move with it: the wall sits exactly where the intrinsic frequency
   $\omega-U(y)k$ is Doppler-shifted to zero.
4. **Same wall, both directions.** Launch with $l_0>0$ and $l_0<0$ from
   the same point (exercise 3 continues this): rays can leave a shear
   region freely but approach their critical layer only asymptotically —
   which is why, once dissipation acts, the momentum they carry is
   deposited *there* and nowhere else. That one-way delivery is the QBO
   mechanism in the "Both fluids" section below.

{{< marimo src="/nb/ch13_wave-mean/" >}}

## Both fluids

- **Ocean:** Stokes drift is directly measured by surface drifters and is
  a leading-order term in how floating debris and larvae are transported —
  distinct from, and typically larger than, the mean Eulerian current at
  the surface.
- **Atmosphere:** critical-layer absorption of upward-propagating waves is
  the accepted driving mechanism of the stratospheric quasi-biennial
  oscillation (QBO) — waves deposit momentum preferentially where their
  phase speed matches the slowly-descending mean wind, reversing that wind
  and shifting the critical layer downward in a self-sustaining cycle.

## Exercises

1. *(analytic)* Redo the Stokes-drift derivation for
   $u'(x,t)=A\sin(kx-\omega t)$ instead of $\cos$ — confirm you get the
   identical result $A^2k/(2\omega)$ (the phase offset cannot matter, since
   $\langle\sin^2\rangle=\langle\cos^2\rangle=1/2$).
2. *(computational)* Modify `ch13_wave-mean.py`'s Part B to use a
   $\tanh$-shaped shear $U(y)=U_0\tanh(y/L)$ instead of linear shear, and
   locate the critical layer numerically (where does $l$ start diverging?)
   rather than from a closed-form $y_c$.
3. *(exploratory)* In Part B, launch several rays with different $l_0$
   (some positive, some negative) from the same $y_0$. Do they all
   approach the same $y_c$, or does $y_c$ depend on $l_0$ through
   $\omega_0$? Explain why launching with the "wrong sign" of $l_0$ might
   send a ray toward a completely different (or no) critical layer.

## Further reading

Bühler, O., *Waves and Mean Flows* (2nd ed.), ch. 1-2 (Stokes drift), ch. 9
(critical layers); Andrews, D. G. & McIntyre, M. E. (1978), *J. Fluid Mech.*
**89**, 609–646; Vallis, *AOFD* 2nd ed., §10.1-10.3; Baldwin, M. P. et al.
(2001), *Rev. Geophys.* **39**, 179–229 (the QBO).
