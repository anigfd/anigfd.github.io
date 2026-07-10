---
title: "Ch 16. Baroclinic instability"
weight: 516
part: 5
---

## Overview

Ch. 15's barotropic instability draws energy from horizontal shear in the
mean flow's kinetic energy. Mid-latitude storms and ocean mesoscale eddies
draw energy from somewhere else entirely: the potential energy locked up in
a sloping temperature surface — the same thermal-wind shear from ch. 5. This
chapter builds the two classic models of that mechanism: the Eady problem,
a boundary-value eigenvalue problem giving a closed growth-rate curve from
stratification and shear alone, and the 2-layer (Phillips) model, whose
nonlinear evolution shows what the instability actually *does* — grow,
peak, break, and equilibrate. It is arguably the single most consequential
instability in geophysical fluid dynamics.

## The model

Notebook: `notebooks/ch16_baroclinic-instability.py` → exported to
`/nb/ch16_baroclinic-instability/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}):

$$\sigma=\mu c_i,\qquad \mu=\frac{kNH}{f_0},\qquad
  q_i=\nabla^2\psi_i+F(\psi_j-\psi_i)+\beta y.$$

`gfdlib.baroclinic` solves both growth-rate problems as small ($2\times2$)
eigenvalue problems, derived directly rather than taken from a memorized
closed form: `eady_growth_rate` matches the Eady interior solution to the
linearized thermal boundary condition at $z=0,H$, verified against the
published benchmark ($\mu_{max}\approx1.61$, $\sigma_{max}\approx0.31$,
cutoff $\mu_c\approx2.399$). `invert_2layer` is a closed-form $2\times2$
Helmholtz solve; `rhs_2layer` is the nonlinear tendency for
`gfdlib.timestep.ifrk4_step`, cross-checked against an independently derived
linear eigenvalue to $<0.1\%$ once transients settle. A fixed mean shear
with no other sink turned out to be an unlimited energy source — verified,
not assumed: without bottom drag, nonlinear growth continued past any
physically reasonable amplitude regardless of hyperviscosity, until
overflow. `rhs_2layer` therefore includes the standard Phillips (1954)
bottom Ekman drag term, which produces a genuine growth-peak-decay-
equilibrate life cycle (confirmed by direct integration to $t=500$ at the
notebook's full resolution). That full run is ~$10^5$ time steps — too
heavy to integrate live in the browser — so it's precomputed offline and
shipped as data for Part B's widget to scrub; a second, smaller live panel
lets you test parameter changes yourself over the growth phase.

{{< marimo src="/nb/ch16_baroclinic-instability/" >}}

## Both fluids

- **Atmosphere:** the storm track — synoptic-scale (~1000 km) cyclones and
  anticyclones extracting available potential energy from the pole-to-equator
  temperature gradient via exactly this mechanism, with an e-folding time of
  order a day.
- **Ocean:** mesoscale eddies (~50-100 km, an order of magnitude smaller,
  set by the much smaller oceanic deformation radius) extracting potential
  energy from sloping density surfaces at ocean fronts and along the edges
  of currents like the Gulf Stream, with e-folding times of weeks rather than
  days — the same instability, the same equations, wildly different scales.

## Exercises

1. *(analytic)* Show that the Eady problem's interior solution
   $\phi(z)=A\cosh(\mu z/H)+B\sinh(\mu z/H)$ follows from $q'=0$ with rigid
   lids, and that matching it to the linearized thermal boundary condition
   at $z=0,H$ gives a $2\times2$ generalized eigenvalue problem for $c$
   rather than a scalar equation.
2. *(computational)* Modify `ch16_baroclinic-instability.py`'s B2 panel to
   sweep $\beta$ at fixed $F,\Delta U$ and locate the marginal-stability
   threshold numerically (where Phillips' necessary condition,
   $\beta_2=\beta-F\Delta U$ changing sign, stops holding). Compare it to
   the analytic threshold $\beta=F\Delta U$.
3. *(exploratory)* In B1's precomputed life cycle, compare the wavenumber
   content of $q_1'$ early in the run (still near-sinusoidal, dominated by
   the seeded mode) to late in the run (post-breaking). What does the
   change tell you about how baroclinic eddies cascade energy once they
   reach finite amplitude?

## Further reading

Eady, E. T. (1949), *Tellus* **1**, 33–52 ("Long waves and cyclone waves");
Phillips, N. A. (1954), *Tellus* **6**, 273–286 ("Energy transformations and
meridional circulations associated with simple baroclinic waves in a
two-level, quasi-geostrophic model"); Vallis, *AOFD* 2nd ed., §6.1-6.8 (the
Eady and Phillips problems); Pedlosky, *Geophysical Fluid Dynamics* 2nd ed.,
§7.1-7.13.
