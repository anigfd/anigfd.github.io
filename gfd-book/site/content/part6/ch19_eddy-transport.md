---
title: "Ch 19. Eddy transport & mixing"
weight: 619
part: 6
---

## Overview

Ch. 18's turbulent eddies don't just cascade energy and enstrophy — they
stir anything riding along with the flow: heat, salt, chemical tracers,
potential vorticity itself. This chapter builds a passive tracer directly
into ch. 18's turbulence machinery and asks the practical question every
coarse ocean or atmosphere model has to answer: if you can't afford to
resolve the eddies, what single number stands in for everything they do to
a tracer? The answer is an *effective diffusivity* — and this notebook
computes one directly from a turbulent simulation, rather than taking it as
a free parameter.

## The model

Notebook: `notebooks/ch19_eddy-transport.py` → exported to
`/nb/ch19_eddy-transport/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}):

$$\frac{\partial c'}{\partial t}+J(\psi,c')=-\Gamma v+\kappa\nabla^2c',
  \qquad K_{eff}=-\frac{\overline{v'c'}}{\Gamma}.$$

`gfdlib.mixing.rhs_coupled` stacks vorticity and tracer perturbation into
one state (the same two-field pattern as ch. 16's two-layer model), reusing
ch. 18's exact turbulence machinery for $\zeta$ unchanged, with the tracer's
mean gradient split off as an ordinary source term — the identical "mean +
periodic perturbation" trick used for $\beta y$ (ch. 18) and mean shear
(ch. 16). Because the perturbation equation is *linear* in $c'$, $K_{eff}$
must be exactly independent of the imposed gradient $\Gamma$ for a fixed
flow — checked numerically (two runs differing only in $\Gamma$ agree to
machine precision) rather than assumed, alongside a control run with the
velocity field set identically to zero, which gives $K_{eff}=0$ exactly (no
possible eddy flux without a flow).

{{< marimo src="/nb/ch19_eddy-transport/" >}}

## Both fluids

- **Atmosphere:** eddy mixing of potential vorticity and chemical tracers
  by synoptic-scale storms sets the effective diffusivity used in every
  coarse-resolution climate and chemistry-transport model — the same
  flux-gradient closure computed directly here.
- **Ocean:** mesoscale eddies (ch. 16's baroclinic instability, at ocean
  scales) are the dominant mechanism transporting heat, carbon, and
  nutrients across ocean fronts and gyre boundaries; every non-eddy-
  resolving ocean climate model parameterizes this transport with an
  eddy diffusivity of exactly the kind measured in this notebook (the
  Gent-McWilliams scheme, in practice).

## Exercises

1. *(analytic)* Substitute $C_{total}=\Gamma y+c'$ into the advection-
   diffusion equation and show the mean-gradient term becomes the source
   $-\Gamma v$ in the perturbation equation for $c'$, with no other change.
2. *(computational)* Modify `ch19_eddy-transport.py` to report the
   time-mean $K_{eff}$ over the last quarter of the run (instead of just
   plotting the full time series) and print it alongside $u_{rms}\sqrt{2E}$
   for a direct mixing-length comparison.
3. *(exploratory)* Run at two different resolutions ($n=64$ and $n=128$,
   same $k_0$, seed, and total simulated time) and compare the late-time
   $K_{eff}$. Does resolving finer-scale filamentation change the measured
   diffusivity much, or does it converge quickly?

## Further reading

Green, J. S. A. (1970), *Q. J. R. Meteorol. Soc.* **96**, 157–185
("Transfer properties of the large-scale eddies and the general circulation
of the atmosphere"); Gent, P. R. & McWilliams, J. C. (1990), *J. Phys.
Oceanogr.* **20**, 150–155 ("Isopycnal mixing in ocean circulation
models"); Vallis, *AOFD* 2nd ed., §10.1-10.3 (eddy diffusivity & mixing
length); Pedlosky, *Geophysical Fluid Dynamics* 2nd ed., §7.14-7.16.
