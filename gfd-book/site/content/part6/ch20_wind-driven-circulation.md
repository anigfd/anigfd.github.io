---
title: "Ch 20. Wind-driven ocean circulation"
weight: 620
part: 6
---

## Overview

The trade winds and westerlies blow steadily over every subtropical ocean
basin, curling clockwise around each gyre. Sverdrup balance — the leading-
order, frictionless response — can only satisfy the no-normal-flow
condition at *one* wall of a closed basin, conventionally the eastern one.
Something has to happen at the western edge to close the circulation. This
chapter builds the two classical answers (Stommel's bottom drag, Munk's
lateral friction) as a genuinely new kind of numerical problem for this
book — a bounded rectangular basin, not a periodic or channel domain — and
shows, from first principles, why western boundary currents like the Gulf
Stream exist and eastern ones don't.

## The model

Notebook: `notebooks/ch20_wind-driven-circulation.py` → exported to
`/nb/ch20_wind-driven-circulation.html`. Governing equations use the symbols in
[notation]({{< relref "notation.md" >}}):

$$\text{Stommel: } \varepsilon\nabla^2\psi+\psi_x=-\sin(\pi y),\qquad
  \text{Munk: } -\delta^3\nabla^4\psi+\psi_x=-\sin(\pi y).$$

`gfdlib.circulation` solves both problems two independent ways: an exact
analytic solution (separating $\psi=f(x)\sin(\pi y)$ reduces each PDE to a
linear ODE for $f(x)$, solved via its characteristic polynomial's roots),
and a general 2D finite-difference direct solve (Kronecker-sum operators,
plain `numpy.linalg.solve`, no SciPy) cross-validated against the analytic
solution with clean 2nd-order convergence. Munk's biharmonic operator is
handled by keeping vorticity as an independent field (ch16/17/19's coupled-
state pattern) closed with Thom's (1933) wall-vorticity formula for no-slip
on the east/west walls. An earlier version enforced no-slip on all four
walls — matching physical intuition but not the separable-solution
assumption — and disagreed with the exact solution by a stubborn ~30% that
did *not* shrink with resolution at any boundary-layer width tested,
immediately distinguishing a genuine discretization error from a
boundary-condition mismatch; the fix (free-slip north/south, the standard
textbook simplification) restored clean convergence.

## What to look for

1. **The asymmetry is the phenomenon.** In either friction model, look at
   the streamfunction: a broad, slow interior drift and one thin, intense
   return current — always on the **west** wall, never the east, no
   matter what you do to the sliders. The interior matches the Sverdrup
   solution $(1-x)\sin(\pi y)$ (overlaid in the mid-basin profile), which
   already "chose" the eastern wall for its boundary condition; friction
   is only allowed to fix what remains, at the west.
2. **Two frictions, two boundary layers.** Toggle Stommel ↔ Munk at
   comparable boundary-layer widths: Stommel's current decays
   monotonically into the interior, while Munk's *overshoots* — a weak
   countercurrent just east of the main jet (the extra structure a
   4th-order operator's oscillatory roots permit). The Gulf Stream's
   observed offshore countercurrent looks distinctly more Munk than
   Stommel.
3. **Squeeze the layer.** Shrink $\varepsilon$ (or $\delta$): the western
   current narrows and its peak velocity grows in inverse proportion,
   because its job is fixed — it must return exactly the transport the
   Sverdrup interior carries equatorward, however thin friction makes it.
   This width–speed trade is why the real Gulf Stream, with its tiny
   effective friction, is so fast and so narrow.
4. **Numerics you can audit.** Flip to the finite-difference solution and
   compare with the analytic curve at increasing resolution: the error
   drops cleanly at 2nd order. The convergence study here is the same
   habit every chapter's test suite applies — and "The model" section's
   story of the all-walls-no-slip bug shows what it catches.

{{< marimo src="/nb/ch20_wind-driven-circulation.html" >}}

## Both fluids

- **Ocean (the primary case):** every subtropical gyre shows exactly this
  asymmetry — broad, slow Sverdrup interior, thin intense western boundary
  current (Gulf Stream, Kuroshio, Agulhas, Brazil Current) — set by the
  same $\beta$-plane balance dialed above, closer to the Munk (lateral-
  friction) limit than the Stommel one for the real ocean.
- **Atmosphere:** the same linear vorticity balance describes the
  depth-averaged, time-mean response to Ekman pumping under the subtropical
  highs, but the atmosphere's turbulent boundary-layer drag and much larger
  deformation radius put it in a different part of parameter space, without
  the ocean's dramatic western intensification — the same equation, a very
  different appearance, depending on which term actually dominates.

## Exercises

1. *(analytic)* Starting from the Sverdrup balance $\psi_x=-\sin(\pi y)$
   with $\psi(1,y)=0$, derive $\psi_{Sv}=(1-x)\sin(\pi y)$ and explain why
   only the eastern boundary condition can be satisfied without friction.
2. *(computational)* Modify `ch20_wind-driven-circulation.py` to plot the
   western boundary layer's peak transport (max of $f(x)$) as a function
   of $\varepsilon$ (or $\delta$), swept over the slider's full range.
   Does it scale the way you'd expect from the boundary-layer width alone?
3. *(exploratory)* Compare the Stommel and Munk profiles at boundary-layer
   widths chosen so their peak transport roughly matches. Which one
   overshoots the Sverdrup interior value more, and can you connect that to
   the order of the underlying ODE (2nd vs. 4th)?

## Further reading

Sverdrup, H. U. (1947), *Proc. Natl. Acad. Sci.* **33**, 318–326 ("Wind-
driven currents in a baroclinic ocean"); Stommel, H. (1948), *Trans. Amer.
Geophys. Union* **29**, 202–206 ("The westward intensification of
wind-driven ocean currents"); Munk, W. H. (1950), *J. Meteorol.* **7**,
79–93 ("On the wind-driven ocean circulation"); Vallis, *AOFD* 2nd ed.,
§19.1-19.4 (wind-driven gyres); Pedlosky, *Geophysical Fluid Dynamics* 2nd
ed., §5.1-5.5.
