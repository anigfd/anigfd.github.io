---
title: "Ch 24. Balanced models & the slow manifold"
weight: 724
part: 7
---

## Overview

The book's closing retrospective, with a runnable demonstration at its
center. Twenty-three chapters solved *reduced* equations while the fluid
obeys none of them exactly — and it worked because rotating, stratified
flow keeps its fast (inertia-gravity) and slow (PV-controlled) lives
separate. A state initialized with zero fast content lives on the **slow
manifold**, and every balanced model in this book is the dynamics
restricted to it: geostrophy at $O(1)$, QG at $O(Ro)$, the planetary
balances at basin scale. The demonstration is Richardson's 1922 forecast
disaster in miniature: one height field, three wind initializations,
identical PV to machine precision — and divergence ringing in exact
proportion to the projection *off* the manifold.

## The model

Notebook: `notebooks/ch24_balanced-models.py` → exported to
`/nb/ch24_balanced-models/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}).

$$(u,v)_{t=0}=\alpha\,(-\eta_{0y},\ \eta_{0x}),\qquad
  \delta=u_x+v_y,\qquad
  \text{RMS}\,\delta\ \propto\ (1-\alpha).$$

Ch. 6's rotating-shallow-water machinery unchanged, plus one new one-line
diagnostic (`gfdlib.shallowwater.divergence`): the slow mode has
$\delta=0$ identically while inertia-gravity waves are made of it, so RMS
divergence directly meters the distance from the slow manifold.

{{< marimo src="/nb/ch24_balanced-models/" >}}

## Both fluids

- **Atmosphere:** balanced initialization descends directly from
  Richardson's failure — every operational forecast system since has
  projected its analysis toward the slow manifold (nonlinear normal-mode
  initialization historically; digital filters and careful data
  assimilation today) before integrating.
- **Ocean:** state estimation faces the same problem with a slower clock —
  an ocean reanalysis initialized off-balance rings with near-inertial
  waves and spurious convection; and the "spontaneous imbalance" caveat
  (the slow manifold leaks gravity waves at $\sim e^{-c/Ro}$) is an active
  research topic precisely because ocean fronts reach $Ro\sim1$, where the
  leak stops being exponential and becomes Ch. 12's internal-wave
  generation in earnest.

## Exercises

1. *(analytic)* Show that for the linear $f$-plane system, the
   geostrophic initialization $u=-\eta_y$, $v=\eta_x$ is an exact steady
   state (substitute into `gfdlib.shallowwater.rhs`'s equations), and
   that an arbitrary initial state's fast-mode amplitude — hence its
   later divergence ringing — is a linear functional of the
   initialization error, giving the notebook's exact $(1-\alpha)$
   scaling.
2. *(computational)* Modify `ch24_balanced-models.py` to add a fourth run
   initialized with the *PV-consistent* balanced state
   `shallowwater.invert_pv(q0)` (Ch. 6's adjustment end-state) instead of
   the raw $\eta_0$ — for a small-scale bump ($\sigma<L_R$), how different
   are "balancing the winds to the observed height" and "balancing
   everything to the observed PV"?
3. *(exploratory)* Turn on $\beta$ (edit the `coriolis` call): the slow
   mode is no longer exactly steady (Rossby waves live on the slow
   manifold and *move*). Can you still cleanly separate slow evolution
   from fast ringing in the divergence record — and what does that tell
   you about why divergence, rather than any time derivative, was the
   right diagnostic all along?

## Further reading

Lorenz, E. N. (1986), *J. Atmos. Sci.* **43**, 1547–1557 ("On the
existence of a slow manifold"); Leith, C. E. (1980), *J. Atmos. Sci.*
**37**, 958–968 (nonlinear normal-mode initialization); Vanneste, J. &
Yavneh, I. (2004), *J. Atmos. Sci.* **61**, 211–223 (exponentially small
imbalance); Vallis, *AOFD* 2nd ed., §5.9 (the slow manifold and balance);
Richardson, L. F. (1922), *Weather Prediction by Numerical Process*
(the origin story).
