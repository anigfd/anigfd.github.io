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
