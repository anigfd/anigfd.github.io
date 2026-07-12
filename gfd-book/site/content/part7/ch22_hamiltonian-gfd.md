---
title: "Ch 22. Hamiltonian & variational GFD"
weight: 722
part: 7
---

## Overview

Integrate the same vortex system with a humble 2nd-order method and with
classical RK4, and watch the energy: RK4's error is far smaller — and grows
forever — while the cruder method's error, if the method carries the right
geometry, stays bounded until the end of time. The geometry is Hamiltonian
structure, and this chapter uses it to answer a question the whole book has
leaned on without answering: GFD's conservation laws are Noether charges of
symmetries, and PV — the book's central quantity since Ch. 7 — is the
charge of the fluid's own peculiar symmetry, **particle relabeling**.

## The model

Notebook: `notebooks/ch22_hamiltonian-gfd.py` → exported to
`/nb/ch22_hamiltonian-gfd/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}).

$$\Gamma_i\frac{dx_i}{dt}=\frac{\partial H}{\partial y_i},\qquad
  \Gamma_i\frac{dy_i}{dt}=-\frac{\partial H}{\partial x_i},\qquad
  H=-\frac{1}{4\pi}\sum_{i<j}\Gamma_i\Gamma_j\ln r_{ij}^2.$$

`gfdlib.vortex` implements the Kirchhoff point-vortex system — an exact
finite-dimensional Hamiltonian reduction of 2D Euler in which positions
themselves are the conjugate pair — with its four Noether invariants
(H, Q, P, L), a symplectic implicit-midpoint integrator, and a same-order
non-symplectic control (Heun). Known solutions (co-rotating pair period,
dipole translation speed) are verified in the test suite.

{{< marimo src="/nb/ch22_hamiltonian-gfd/" >}}

## Both fluids

- **Atmosphere & ocean together:** this chapter's payoff is not one
  phenomenon but a design principle for the models of both fluids. Climate
  and ocean models run for centuries of model time — far beyond trajectory
  accuracy — and what keeps their statistics physical is conservation
  structure (mimetic grids, conservative advection, careful energy/enstrophy
  bookkeeping descended from Arakawa's schemes), exactly the
  structure-vs-accuracy distinction the energy-drift experiment makes
  visible in miniature.
- **Ocean, concretely:** interacting mesoscale eddies are routinely modeled
  as point-vortex systems; the 3-vortex integrable / 4-vortex chaotic
  transition sets genuine limits on predicting eddy encounters (e.g.
  Gulf Stream ring interactions).

## Exercises

1. *(analytic)* Verify by direct differentiation that Kirchhoff's
   $H=-\frac{1}{4\pi}\sum_{i<j}\Gamma_i\Gamma_j\ln r_{ij}^2$ generates
   exactly `gfdlib.vortex.vortex_rhs` through
   $\Gamma_i\dot x_i=\partial H/\partial y_i$,
   $\Gamma_i\dot y_i=-\partial H/\partial x_i$, and show
   $\frac{d}{dt}\sum_i\Gamma_ix_i=0$ directly from the antisymmetry of the
   double sum.
2. *(computational)* Modify `ch22_hamiltonian-gfd.py` to add leapfrog
   applied naively to $(x,y)$ as a fourth method. It fails (leapfrog's
   symplecticity requires a separable Hamiltonian, and Kirchhoff's is
   maximally non-separable) — document *how* it fails compared to Heun.
3. *(exploratory)* For the chaotic 4-vortex configuration, find the
   $(\Delta t, T)$ combination at which RK4's secular energy error first
   crosses above midpoint's bounded band. How does that crossing time scale
   as you halve $\Delta t$? (Prediction: RK4's error $\sim\Delta t^4t$,
   midpoint's band $\sim\Delta t^2$, so $t_\times\sim\Delta t^{-2}$.)

## Further reading

Salmon, R. (1988), *Ann. Rev. Fluid Mech.* **20**, 225–256 ("Hamiltonian
fluid mechanics"); Salmon, *Lectures on Geophysical Fluid Dynamics*, ch. 7;
Shepherd, T. G. (1990), *Adv. Geophys.* **32**, 287–338; Aref, H. (1983),
*Ann. Rev. Fluid Mech.* **15**, 345–389 (point-vortex motion); Hairer,
Lubich & Wanner, *Geometric Numerical Integration* 2nd ed., ch. VI
(symplectic methods, backward error analysis).
