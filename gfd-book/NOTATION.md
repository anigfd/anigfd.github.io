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

## Rayleigh-Benard convection (vertical $x$–$z$ channel)
| Symbol | Meaning | Convention |
|---|---|---|
| $\psi(x,z)$ | channel streamfunction | $u=-\psi_z,\ w=\psi_x$ |
| $\zeta=\nabla^2\psi$ | vorticity | $w_x-u_z=\psi_{xx}+\psi_{zz}=\nabla^2\psi$; **same sign as the horizontal case above** |
| $\theta$ | temperature deviation from the linear conduction profile | |
| $Ra_c=27\pi^4/4$ | critical Rayleigh number (free-slip, both boundaries) | onset at $k_c=\pi/\sqrt2$ |

> **A genuine bug to avoid re-introducing.** With $u=-\psi_z,\ w=\psi_x$, vorticity
> is $\zeta=w_x-u_z=\nabla^2\psi$ — the **same** sign as the horizontal
> streamfunction convention above, *not* $\nabla^2\psi=-\zeta$. Getting this
> backward flips the sign of the buoyancy-vorticity feedback and the
> instability never grows at any $Ra$: verified numerically (see ch14 PR) that
> the wrong sign gives pure decay for all $Ra$, while the correct sign
> reproduces subcritical decay ($Ra<Ra_c$) and supercritical growth into
> finite-amplitude convection ($Ra>Ra_c$) as expected.
> `gfdlib.convection.ChannelGrid` (periodic $x$, Dirichlet $z$) and
> `rhs_boussinesq` implement the corrected convention; `lorenz_rhs` is the
> $(X,Y,Z)$ truncation, sharing `gfdlib.timestep.rk4` with every other chapter.

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
