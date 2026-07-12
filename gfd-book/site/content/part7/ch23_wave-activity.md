---
title: "Ch 23. Conservation laws & wave activity"
weight: 723
part: 7
---

## Overview

Can waves accelerate the jet they ride on? For conservative,
small-amplitude waves the answer is a theorem, and the answer is *no*:
zonal-mean flow and wave **pseudomomentum** change in exact opposition,
$\partial_t(\bar u+A)=0$ pointwise in latitude — the **non-acceleration
theorem**. Waves are a loan of momentum, repaid in full unless the wave
dies; only transience and dissipation (Ch. 13's critical layers being the
canonical executioner) convert the loan into a grant. This chapter derives
the theorem in three short steps — Noether's zonal-translation charge, the
Taylor identity, the zonal-mean momentum equation — and then verifies it
pointwise against the full nonlinear barotropic model.

## The model

Notebook: `notebooks/ch23_wave-activity.py` → exported to
`/nb/ch23_wave-activity/`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}).

$$A=\frac{\overline{q'^2}}{2\bar q_y},\qquad
  \frac{\partial A}{\partial t}=-\overline{v'q'}
  =\frac{\partial}{\partial y}\overline{u'v'}=-\frac{\partial\bar u}{\partial t}
  \;\;\Longrightarrow\;\;
  \frac{\partial}{\partial t}(\bar u+A)=0.$$

A small Rossby-wave packet on a sinusoidal shear $\bar U=U_0\cos y$ is
evolved with Ch. 18's full nonlinear pseudo-spectral solver;
`gfdlib.wavemean.wave_activity` diagnoses $\bar u(y)$, $A(y)$, and
$\bar q_y(y)$ from every snapshot. The residual $\Delta(\bar u+A)$ comes
out below 1% of either term (also enforced in the test suite), and a
hyperviscosity slider provides the theorem's off switch.

{{< marimo src="/nb/ch23_wave-activity/" >}}

## Both fluids

- **Atmosphere:** the quasi-biennial oscillation is non-acceleration's
  loophole running as an engine — waves launched by tropical convection
  break at slowly-descending critical lines, depositing momentum that
  reverses the stratospheric winds every ~14 months; sudden stratospheric
  warmings are the same ledger settled violently in midwinter. The
  Eliassen–Palm flux plotted daily in stratosphere diagnostics is exactly
  this chapter's wave-activity flux, generalized to stratified flow.
- **Ocean:** eddy parameterizations built on PV fluxes (Ch. 19's
  caution about up-gradient fluxes included) are constrained by the same
  identity — where eddies are steady and conservative, their momentum
  forcing of the mean gyres must vanish, however energetic they look.

## Exercises

1. *(analytic)* Derive the Taylor identity
   $\overline{v'q'}=-\partial_y\overline{u'v'}$ from $q'=\nabla^2\psi'$,
   $u'=-\psi'_y$, $v'=\psi'_x$: expand
   $\overline{\psi'_x(\psi'_{xx}+\psi'_{yy})}$, and show every term is
   either a pure $x$-derivative of a product (zonal average zero) or
   $-\partial_y(\overline{u'v'})$.
2. *(computational)* Modify `ch23_wave-activity.py` to also plot the
   displacement form of wave activity,
   $\tfrac12\bar q_y\overline{\eta^2}$ with $\eta=-q'/\bar q_y$, and
   confirm it is identical to $A$ (it is algebraically the same object —
   this is a consistency check on the diagnostics, and a good way to
   internalize that $A$ measures contour displacement).
3. *(exploratory)* Using the $\nu$ slider as a dial rather than a switch,
   measure the final-time residual $\max_y|\Delta(\bar u+A)|$ as a
   function of $\nu$. Is the momentum "lost" to dissipation proportional
   to $\nu$ at small $\nu$, and where in $y$ does the mean flow actually
   end up changed?

## Further reading

Andrews, D. G. & McIntyre, M. E. (1976), *J. Atmos. Sci.* **33**,
2031–2048; Charney, J. G. & Drazin, P. G. (1961), *J. Geophys. Res.*
**66**, 83–109; Bühler, *Waves and Mean Flows* 2nd ed., ch. 4–5; Vallis,
*AOFD* 2nd ed., ch. 10 (wave–mean-flow interaction), §10.4 (non-
acceleration); Baldwin et al. (2001), *Rev. Geophys.* **39**, 179–229
(the QBO).
