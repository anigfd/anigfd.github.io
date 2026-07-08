---
title: "Notation — symbols & sign conventions"
---

# NOTATION — symbols & sign conventions

*Linked from every chapter. If the printed equation and the code disagree, the
code is wrong. Fix conventions here first.*

## Coordinates & operators
| Symbol | Meaning |
|---|---|
| $x,y,z$ | eastward, northward, upward |
| $\mathbf{u}=(u,v,w)$ | velocity |
| $\nabla^2$ | horizontal Laplacian (spectral: $-k^2$) |
| $J(A,B)=A_xB_y-A_yB_x$ | horizontal Jacobian |

## Rotation & stratification
| Symbol | Meaning | Convention |
|---|---|---|
| $f=f_0+\beta y$ | Coriolis parameter | $\beta>0$ (increases northward) |
| $\Omega$ | planetary rotation rate | |
| $N^2=-\frac{g}{\rho_0}\frac{\partial\rho}{\partial z}$ | buoyancy freq. squared | $N^2>0$ stable |
| $b=-g\rho'/\rho_0$ | buoyancy | |

## Dynamical fields
| Symbol | Meaning | Sign |
|---|---|---|
| $\psi$ | streamfunction | $\mathbf{u}=(-\psi_y,\psi_x)$ |
| $\zeta=\nabla^2\psi$ | relative vorticity | $\nabla^2\psi=\zeta$; invert with `invert_laplacian(zeta_hat)` |
| $q$ | potential vorticity | Ertel / QG per chapter |

> **Streamfunction sign.** We use $\mathbf u=\hat{\mathbf z}\times\nabla\psi=(-\psi_y,\psi_x)$
> and $\zeta=\nabla^2\psi$. In vorticity-form notebooks solve $\nabla^2\psi=\zeta$
> (so `psi_hat = grid.invert_laplacian(zeta_hat)`). Keep this consistent across all notebooks.

## Shallow water (linearized, beta-plane)
| Symbol | Meaning | Convention |
|---|---|---|
| $H$ | mean layer depth | |
| $\eta$ | free-surface height anomaly | positive = raised surface |
| $c=\sqrt{gH}$ | gravity-wave speed | |
| $L_R=c/f_0$ | Rossby deformation radius | length scale |

Nondimensionalized by $L_R$ (length), $f_0^{-1}$ (time), $H$ (height, so
$\hat\eta=\eta/H$, $\hat c=1$), with $\hat f = 1+\hat\beta(y-y_0)$,
$\hat\beta=\beta L_R/f_0$:

$$\hat\eta_t+\hat u_x+\hat v_y=0,\qquad
\hat u_t-\hat f\hat v=-\hat\eta_x,\qquad
\hat v_t+\hat f\hat u=-\hat\eta_y.$$

**Linear PV anomaly** $q=\zeta-\hat\eta$ (with $\zeta=v_x-u_y$ from the same
$\hat u,\hat v$) obeys $q_t=-\hat\beta\hat v$: exactly conserved pointwise on
the $f$-plane ($\hat\beta=0$). `gfdlib.shallowwater.potential_vorticity`
computes it; `gfdlib.shallowwater.invert_pv` inverts $(\nabla^2-1)\eta_{bal}=q$
for the geostrophically balanced height field carrying that PV.

## Dimensionless numbers
| Number | Definition | Regime |
|---|---|---|
| Rossby $Ro$ | $U/(f_0 L)$ | $\ll1$ rotation-dominated |
| Burger $Bu$ | $(NH/f_0L)^2=(L_R/L)^2$ | balance vs. wave scale |
| Ekman $Ek$ | $\nu/(f_0H^2)$ | friction vs. rotation |
| Froude $Fr$ | $U/(NH)$ | stratification |
| Rayleigh $Ra$ | $g\alpha\Delta T H^3/(\nu\kappa)$ | convective onset |
| Reynolds $Re$ | $UL/\nu$ | inertia vs. viscosity |
| Richardson $Ri$ | $N^2/(\partial_z u)^2$ | shear stability |
