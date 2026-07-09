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

## Geostrophic, hydrostatic & thermal-wind balance
| Symbol | Meaning | Convention |
|---|---|---|
| $\phi$ | kinematic pressure, $p/\rho_0$ | |
| $\mathbf u_g=(u_g,v_g)$ | geostrophic velocity | $u_g=-\phi_y/f,\ v_g=\phi_x/f$ |
| hydrostatic balance | $\phi_z=b$ | |
| thermal wind | $fu_{g,z}=-b_y,\quad fv_{g,z}=b_x$ | vertical shear from the horizontal buoyancy gradient |

> **Geostrophic balance is $\psi=\phi/f$.** $u_g=-\phi_y/f,\ v_g=\phi_x/f$ has
> exactly the streamfunction form above with $\psi=\phi/f$ — geostrophic flow
> always has a streamfunction, and $\phi/f$ is it. `gfdlib.balance` works in a
> single $(y,z)$ cross-section with a closed-form frontal buoyancy field (no
> discretization): `buoyancy_field` and its exact $y$-derivative
> `dbdy_frontal`, and `thermal_wind_u`, the exact vertical integral of the
> thermal-wind relation. The buoyancy front's sign reverses at $z=H_{trop}$
> (an idealized tropopause), so $u_g(y,z)$ has an extremum there at every
> $y$ — the geostrophic jet core sits where the meridional buoyancy gradient
> changes sign.

## PV inversion & staircases (barotropic, $f$- or $\beta$-plane)
| Symbol | Meaning | Convention |
|---|---|---|
| $q=\zeta+\beta y$ | (absolute) potential vorticity | material conservation: $q_t+J(\psi,q)=0$ |

> **The invertibility principle.** Given $q(x,y)$ everywhere, a balance
> condition ($\zeta=\nabla^2\psi$, i.e. nondivergent flow), and boundary
> conditions (here: periodic), $\psi$ — and hence $\mathbf u$ — is uniquely
> determined: `psi_hat = grid.invert_laplacian(zeta_hat)` with
> $\zeta=q-\beta y$. `gfdlib.pv.gaussian_blob` places a periodic-safe PV
> anomaly at a point; `gfdlib.pv.staircase_pv` builds an idealized
> alternating-band PV staircase (a smoothed square wave in $y$). Inverting a
> staircase shows the jet cores ($u$ extrema) sit exactly at the sharp PV
> risers between bands, not within the well-mixed plateaus — verified
> numerically (local $|q_y|$ at each jet core is $\sim$20$\times$ the
> domain-mean $|q_y|$).

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

## Single-layer QG (barotropic PV + deformation-radius stretching)
| Symbol | Meaning | Convention |
|---|---|---|
| $q=\nabla^2\psi-\psi+\beta y$ | quasi-geostrophic PV ($L_R$ nondim to 1) | material conservation: $q_t+J(\psi,q)=0$ |

> **QGPV inversion reuses the shallow-water Helmholtz solve.** With
> $L_R$ nondimensionalized to 1, $(\nabla^2-1)\psi=q-\beta y$ is *the same
> operator* as the shallow-water PV inversion above — `gfdlib.shallowwater.
> invert_pv` computes $\psi$ from the QGPV anomaly directly; no separate QG
> inversion primitive exists or is needed. The linear Rossby-wave part of
> the evolution operator changes accordingly: `gfdlib.qg.dispersion_omega`
> gives $\omega=-\beta k_x/(k_x^2+k_y^2+1)$, bounded by $|\omega|\le\beta/2$
> (attained at $|\mathbf k|=1$) — unlike the unbounded barotropic relation
> in ch18, the deformation radius screens the lowest wavenumbers. An
> isolated vortex (`gfdlib.pv.gaussian_blob`) on this operator sheds a
> trailing Rossby-wave wake and drifts westward — verified numerically.

## Rossby-wave ray tracing (sphere, WKB)
| Symbol | Meaning | Convention |
|---|---|---|
| $\lambda,\phi$ | longitude, latitude | |
| $n,m$ | zonal/meridional wavenumber indices, conjugate to $\lambda,\phi$ | physical wavenumbers $k_x=n/(a\cos\phi),\ k_y=m/a$ |
| $\Omega_s$ | background solid-body rotation rate (added to $\Omega$) | $U(\phi)=\Omega_s a\cos\phi$ |

> **Great circles are a theorem here, not a plotting choice.** For a
> background flow $U=\Omega_s a\cos\phi$, both $U$ and the meridional
> gradient of absolute vorticity scale by $(\Omega+\Omega_s)$ — the same
> factor — so the stationary wavenumber $K_s^2=2(\Omega+\Omega_s)/(\Omega_s
> a^2)$ is **independent of latitude**, and Hamilton's ray equations
> (`gfdlib.rossby.ray_rhs`, via central differences of `dispersion_omega` —
> deliberately not hand-derived symbolically) integrate to exact great
> circles: `great_circle_deviation` checks this directly (fit the plane
> through the ray's first two sampled points on the unit sphere, then check
> every later point lies in it) and finds deviations at the $10^{-9}$
> level, i.e. floating-point roundoff, not approximation error.

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

## Internal gravity waves (vertical $x$–$z$ plane)
| Symbol | Meaning | Convention |
|---|---|---|
| $\psi(x,z)$ | vertical-plane streamfunction | $u=\psi_z,\ w=-\psi_x$ |
| $q=\nabla^2\psi$ | as above, but in the $(x,z)$ plane | $\partial_t^2q=-(N^2(z)\psi_{xx}+f^2\psi_{zz})+S$ |
| $\theta$ | wavevector angle from vertical | $\cos\theta=\sqrt{(\omega^2-f^2)/(N^2-f^2)}$; equals the energy beam's angle from horizontal |

> **A different $\psi$ than the horizontal chapters.** This $(u,w)=(\psi_z,-\psi_x)$
> sign convention is for a *vertical* plane and is unrelated to the horizontal
> $(x,y)$ streamfunction above — the two never appear in the same notebook.
> `gfdlib.internalwaves` reuses `spectral.Grid` with its $y$-axis standing in
> for $z$; `dispersion_omega`/`beam_angle` give $\omega(k,m)$ and $\theta$;
> `step_leapfrog` integrates the $q$ equation and accepts a $z$-dependent
> $N^2$ directly (no assumption of constant stratification). Nondimensionalized
> by a reference buoyancy frequency $N_0$ (time) and the domain width (length):
> $\hat N=N/N_0,\ \hat f=f/N_0,\ \hat\omega=\omega/N_0$.

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
