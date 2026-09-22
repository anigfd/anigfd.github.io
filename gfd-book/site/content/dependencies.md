---
title: "Chapter dependency map"
weight: 2
---

Chapters are written to be read in order, but nobody reads a textbook in
order. This page says what each chapter actually assumes, so you can jump
to the topic you need and know exactly what to backfill.

## The shape of the book

```
Part I    Foundations           1 → 2 → 3 → 4          (everything starts here)
Part II   Balance & adjustment  4 → 5 → 6 → 7          (the load-bearing part)
Part III  Quasi-geostrophy      7 → 8 → 9, 10
Part IV   Stratification/waves  10 → 11 → 12 → 13
Part V    Instabilities         14; 7 → 15 → 16 → 17
Part VI   Turbulence & circ.    8 → 18 → 19; 7 → 20; 14 → 21
Part VII  Structure (optional)  22, 23, 24 — never a prerequisite for anything
```

Two chapters carry more weight than any others: **Ch. 6** (rotating
shallow water — the source of the deformation radius, PV, and adjustment)
and **Ch. 7** (vorticity and PV inversion — the two-step algorithm behind
every simulation from Ch. 8 onward). If you are lost anywhere in the book,
the odds are good the missing piece is in one of those two.

## What each chapter builds on

"Builds on" means *concepts you need*, not merely code that gets reused —
the shared `gfdlib` machinery travels invisibly. Parenthetical chapters
are helpful but skippable.

| Chapter | Builds on |
|---|---|
| [1. Kinematics & the material derivative]({{< relref "part1/ch01_kinematics.md" >}}) | — (start here) |
| [2. Governing equations & Boussinesq]({{< relref "part1/ch02_governing-equations.md" >}}) | [1]({{< relref "part1/ch01_kinematics.md" >}}) |
| [3. Rotation]({{< relref "part1/ch03_rotation.md" >}}) | [2]({{< relref "part1/ch02_governing-equations.md" >}}) |
| [4. Scaling & the Ro–Bu regime map]({{< relref "part1/ch04_scaling.md" >}}) | [2]({{< relref "part1/ch02_governing-equations.md" >}}), [3]({{< relref "part1/ch03_rotation.md" >}}) |
| [5. Geostrophic & hydrostatic balance]({{< relref "part2/ch05_geostrophic-balance.md" >}}) | [2]({{< relref "part1/ch02_governing-equations.md" >}}), [4]({{< relref "part1/ch04_scaling.md" >}}) |
| [6. Rotating shallow water & adjustment]({{< relref "part2/ch06_geostrophic-adjustment.md" >}}) | [3]({{< relref "part1/ch03_rotation.md" >}}), [5]({{< relref "part2/ch05_geostrophic-balance.md" >}}) |
| [7. Vorticity & potential vorticity]({{< relref "part2/ch07_vorticity-pv.md" >}}) | [1]({{< relref "part1/ch01_kinematics.md" >}}), [6]({{< relref "part2/ch06_geostrophic-adjustment.md" >}}) |
| [8. The QG approximation]({{< relref "part3/ch08_qg-approximation.md" >}}) | [6]({{< relref "part2/ch06_geostrophic-adjustment.md" >}}), [7]({{< relref "part2/ch07_vorticity-pv.md" >}}) |
| [9. Rossby waves & ray tracing]({{< relref "part3/ch09_rossby-waves.md" >}}) | [8]({{< relref "part3/ch08_qg-approximation.md" >}}) |
| [10. Two-layer QG & APE]({{< relref "part3/ch10_two-layer-qg.md" >}}) | [7]({{< relref "part2/ch07_vorticity-pv.md" >}}), [8]({{< relref "part3/ch08_qg-approximation.md" >}}) |
| [11. Stratification & vertical modes]({{< relref "part4/ch11_stratification.md" >}}) | [10]({{< relref "part3/ch10_two-layer-qg.md" >}}) |
| [12. Internal gravity waves]({{< relref "part4/ch12_internal-gravity-waves.md" >}}) | [11]({{< relref "part4/ch11_stratification.md" >}}) |
| [13. Wave–mean interaction]({{< relref "part4/ch13_wave-mean.md" >}}) | [9]({{< relref "part3/ch09_rossby-waves.md" >}}), [12]({{< relref "part4/ch12_internal-gravity-waves.md" >}}) ([1]({{< relref "part1/ch01_kinematics.md" >}})) |
| [14. Convection & Lorenz chaos]({{< relref "part5/ch14_convection-lorenz.md" >}}) | [2]({{< relref "part1/ch02_governing-equations.md" >}}) |
| [15. Barotropic instability]({{< relref "part5/ch15_barotropic-instability.md" >}}) | [7]({{< relref "part2/ch07_vorticity-pv.md" >}}) |
| [16. Baroclinic instability]({{< relref "part5/ch16_baroclinic-instability.md" >}}) | [5]({{< relref "part2/ch05_geostrophic-balance.md" >}}), [10]({{< relref "part3/ch10_two-layer-qg.md" >}}), [15]({{< relref "part5/ch15_barotropic-instability.md" >}}) |
| [17. Symmetric & other instabilities]({{< relref "part5/ch17_other-instabilities.md" >}}) | [5]({{< relref "part2/ch05_geostrophic-balance.md" >}}), [15]({{< relref "part5/ch15_barotropic-instability.md" >}}) |
| [18. Geostrophic turbulence]({{< relref "part6/ch18_geostrophic-turbulence.md" >}}) | [7]({{< relref "part2/ch07_vorticity-pv.md" >}}), [8]({{< relref "part3/ch08_qg-approximation.md" >}}) |
| [19. Eddy transport & mixing]({{< relref "part6/ch19_eddy-transport.md" >}}) | [18]({{< relref "part6/ch18_geostrophic-turbulence.md" >}}) |
| [20. Wind-driven gyres]({{< relref "part6/ch20_wind-driven-circulation.md" >}}) | [7]({{< relref "part2/ch07_vorticity-pv.md" >}}) |
| [21. Overturning circulation]({{< relref "part6/ch21_overturning.md" >}}) | [14]({{< relref "part5/ch14_convection-lorenz.md" >}}) |
| [22. Hamiltonian GFD]({{< relref "part7/ch22_hamiltonian-gfd.md" >}}) | [7]({{< relref "part2/ch07_vorticity-pv.md" >}}) |
| [23. Wave activity & non-acceleration]({{< relref "part7/ch23_wave-activity.md" >}}) | [8]({{< relref "part3/ch08_qg-approximation.md" >}}), [13]({{< relref "part4/ch13_wave-mean.md" >}}) |
| [24. Balanced models & the slow manifold]({{< relref "part7/ch24_balanced-models.md" >}}) | [6]({{< relref "part2/ch06_geostrophic-adjustment.md" >}}), [8]({{< relref "part3/ch08_qg-approximation.md" >}}) |

## Shortest paths to the destinations people actually want

- **Baroclinic instability (Ch. 16, the book's centerpiece):**
  1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 10 → 15 → 16. Ten chapters, and every
  one earns its place — this *is* the book's main sequence.
- **Geostrophic turbulence and jets (Ch. 18):**
  1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 18. You can defer Parts IV and V
  entirely.
- **Chaos and the Lorenz system (Ch. 14):** 1 → 2 → 14. The shortest
  worthwhile detour in the book — no rotation required.
- **Internal waves and mixing (Chs. 12–13):** the Part I–II spine, then
  10 → 11 → 12 → 13 (with Ch. 9 before Ch. 13's ray-tracing half).
- **The ocean circulation chapters (Chs. 20–21):** Ch. 20 needs only the
  spine through Ch. 7 (plus $\beta$ intuition from Ch. 8); Ch. 21 needs
  Ch. 14's convection machinery.

Part VII is the optional capstone: it *revisits* the whole book through
Hamiltonian and wave–mean-flow structure, so it draws on everything — but
nothing draws on it. Skipping it costs no later chapter.
