---
title: "Ch 14. Convection, Rayleigh–Bénard & Lorenz chaos"
weight: 514
part: 5
---

## Overview

Heat a fluid layer from below and nothing happens — until it does. Below a
critical temperature contrast, conduction quietly carries the heat and the
fluid stays at rest; cross the threshold and convection rolls appear,
carrying heat far more efficiently. This chapter derives that threshold from
the linearized Boussinesq equations, then takes the same equations somewhere
unexpected: keeping only the three most energetic modes turns the full PDE
into three coupled ODEs — the Lorenz (1963) system — whose behavior as you
turn up the heating traces the entire route from steady rolls to deterministic
chaos. This is the calculation that founded chaos theory, and it comes
directly from a fluid-dynamics problem, not an abstract map.

## The model

Notebook: `notebooks/ch14_convection-lorenz.py` → exported to `/nb/ch14_convection-lorenz.html`.
Governing equations use the symbols in [notation]({{< relref "notation.md" >}}):

$$\zeta_t+J(\psi,\zeta)=Pr\nabla^2\zeta+Pr\,Ra\,\theta_x,\qquad
  \theta_t+J(\psi,\theta)=\nabla^2\theta+\psi_x,\qquad \nabla^2\psi=\zeta.$$

Two connected models, both built on `gfdlib.convection`: the full 2D PDE on
a periodic-$x$/Dirichlet-$z$ channel (`ChannelGrid`, a pure-NumPy tridiagonal
Poisson solve — no SciPy), and its three-variable Lorenz truncation
$(X,Y,Z)$, integrated with the same `gfdlib.timestep.rk4` used everywhere
else in this book. A new interactive **bifurcation diagram** sweeps $r$
across hundreds of values simultaneously (vectorized, not looped) and plots
every local maximum of $Z(t)$, making the transition to chaos at
$r_H\approx24.74$ visible in one figure.

## What to look for

1. **A genuine threshold.** Run the 2D simulation just below and just
   above $Ra_c=27\pi^4/4\approx657.5$ (the free-slip critical value —
   the slider crosses it): below, the seeded perturbation dies and pure
   conduction remains; above, rolls emerge with horizontal wavelength
   near $2\sqrt2\,H$, the $k_c=\pi/\sqrt2$ mode that goes unstable first.
   Instability as a bifurcation, not a matter of degree.
2. **Watch the truncation be good, then fail.** At modest supercriticality
   the Lorenz variables track the 2D simulation's story faithfully —
   steady rolls are the fixed points $C_\pm$. Push $r$ (the Lorenz stand-in
   for $Ra/Ra_c$) past $\approx24.74$ and the trajectory abandons the
   fixed points for the butterfly. The three-mode model is no longer
   quantitatively the PDE — but the qualitative lesson it discovered is
   the one that matters.
3. **Sensitive dependence, measured.** Run several trajectories with the
   `trajectories` slider: initial conditions differing by $10^{-8}$ stay
   visually identical for a while, then separate completely — and the
   separation panel shows the gap growing along a straight line on the
   log axis (the Lyapunov exponent) until it saturates at attractor size.
   That exponential is why weather forecasts gain so little lead time from
   enormous gains in observation accuracy.
4. **Order inside chaos.** In the bifurcation diagram, look past the onset
   at $r_H\approx24.74$ for narrow vertical windows where the smear of
   $Z$-maxima briefly collapses toward a few clean branches — periodic
   windows, the fine structure exercise 3 sends you hunting for with a
   narrower sweep and more $r$ values.

{{< marimo src="/nb/ch14_convection-lorenz.html" >}}

## Both fluids

- **Atmosphere:** Rayleigh-Bénard convection is the textbook idealization of
  boundary-layer thermals and shallow cumulus convection — warm air rising,
  cool air sinking, in cells set by the same onset physics derived here.
- **Ocean:** the same instability drives open-ocean deep convection (e.g., in
  the Labrador and Weddell Seas) where surface cooling destabilizes the water
  column and triggers overturning plumes that ventilate the deep ocean.

## Exercises

1. *(analytic)* Derive the growth rate $\sigma(k)$ from the linearized
   equations and show $\sigma=0$ gives $Ra_c(k)=(k^2+\pi^2)^3/k^2$; confirm
   $Ra_c$ is minimized at $k_c=\pi/\sqrt2$.
2. *(computational)* Modify `ch14_convection-lorenz.py`'s 2D simulation to
   use rigid (no-slip) rather than free-slip boundaries by changing the
   vorticity boundary condition. How does $Ra_c$ shift?
3. *(exploratory)* In the bifurcation diagram, hunt for periodic windows at
   $r>24.74$ (narrow ranges where the smear briefly collapses to a simple
   curve). Zoom in by narrowing `bif_rmax` and increasing the number of $r$
   values around a candidate window.

## Further reading

Lorenz, E. N. (1963). *Deterministic Nonperiodic Flow*, J. Atmos. Sci. **20**,
130–141; Saltzman, B. (1962), J. Atmos. Sci. **19**, 329–341; Chandrasekhar,
S. (1961). *Hydrodynamic and Hydromagnetic Stability*; Strogatz, S. H. (1994).
*Nonlinear Dynamics and Chaos*; Vallis, *AOFD* 2nd ed., §9.1–9.3.
