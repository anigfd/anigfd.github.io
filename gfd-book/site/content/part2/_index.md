---
title: "Part 2 — Balance & Adjustment"
weight: 200
description: "How a rotating, stratified fluid finds its balance — thermal wind, geostrophic adjustment, and the PV-inversion algorithm that drives nearly everything after it."
---

How a rotating, stratified fluid finds its balance. **Thermal wind** ties
vertical shear directly to a horizontal temperature gradient; **geostrophic
adjustment** shows how an unbalanced disturbance sheds fast gravity waves
and settles into a balanced remnant, predictable in advance from
potential vorticity alone. The two-step **PV-inversion algorithm**
introduced here — invert for velocity, advect PV, invert again — drives
nearly every simulation in the rest of the book.
