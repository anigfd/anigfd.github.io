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

## Barotropic instability (Rayleigh-Kuo)
| Symbol | Meaning | Convention |
|---|---|---|
| $U(y)$ | basic-state zonal flow | channel, rigid walls $\phi=0$ |
| $c=c_r+ic_i$ | complex phase speed of a normal mode $\phi(y)e^{ik(x-ct)}$ | unstable iff $c_i>0$ |

> **Rayleigh-Kuo equation:** $(U-c)(\phi''-k^2\phi)+(\beta-U'')\phi=0$ —
> derived independently in this repo (not copied from a text) and checked
> against the classical form at $\beta=0$ (Drazin & Reid). A **necessary**
> condition for instability is that $\beta-U''$ changes sign somewhere in
> the domain (Rayleigh's inflection-point criterion, generalized to a
> $\beta$-plane by Kuo). `gfdlib.instability.growth_rate` solves this as a
> generalized eigenvalue problem $A\phi=cB\phi$, rewritten $B^{-1}A$ (no
> SciPy needed) and diagonalized with plain `numpy.linalg.eigvals`.
> Verified: a profile with no inflection point gives exactly zero growth at
> every $k$ tested; the classical $\tanh(y/\delta)$ shear layer's growth
> rate peaks at $k\delta\approx0.445$, matching the published value
> (Michalke 1964, $k\delta\approx0.4446$) to within grid resolution.
> `gfdlib.instability.double_shear_layer` builds the periodic nonlinear
> initial condition; its evolution reuses `spectral.Grid` and
> `timestep.ifrk4_step` unchanged from ch18.

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

## Baroclinic instability (Eady; 2-layer/Phillips)
| Symbol | Meaning | Convention |
|---|---|---|
| $\mu=kNH/f_0$ | Eady nondimensional horizontal wavenumber | growth rate $\sigma=\mu c_i$ |
| $F=f_0^2/(g'H)$ | layer coupling ("stretching") parameter, 2-layer model | |
| $U_1,U_2$ | upper/lower layer mean flow | $\Delta U=U_1-U_2$ drives the instability |
| $r$ | bottom Ekman drag coefficient | acts on layer 2's relative vorticity only |

> **Both growth-rate problems are solved as small eigenvalue problems
> derived here, not copied from a closed-form formula.** `gfdlib.
> baroclinic.eady_growth_rate` matches the boundary-value problem's
> interior solution to the linearized thermal boundary condition at
> $z=0,H$; verified against the published benchmark
> ($\mu_{max}\approx1.61$, $\sigma_{max}\approx0.31$,
> cutoff $\mu_c\approx2.399$; Vallis/Pedlosky) to within grid resolution.
> `gfdlib.baroclinic.invert_2layer` is the 2-layer QGPV Helmholtz solve
> (closed-form $2\times2$ inverse); `rhs_2layer` is spectral-space, for use
> as the nonlinear RHS in `gfdlib.timestep.ifrk4_step` with hyperviscosity
> as the (scalar, both-layers-identical) linear operator — the mean-flow
> and $\beta$ terms stay in the nonlinear part since they couple the two
> layers ($2\times2$ per wavenumber, not a scalar). Cross-checked against
> an independently-derived linear eigenvalue: a small-amplitude
> single-wavenumber perturbation's measured growth rate converges to the
> predicted value to $<0.1\%$ once the mode purifies (a coincidental choice
> of parameters that placed an integer wavenumber exactly at the
> marginal-stability cutoff was caught and avoided during this check).
>
> **A fixed mean shear with no drag is an unlimited energy source.**
> Nonlinearly, `rhs_2layer` without Ekman drag ($r=0$) was found to grow
> without bound regardless of how much hyperviscosity was added — energy
> conservation was verified exact for the pure nonlinear (Jacobian) terms
> in isolation, so the runaway traced to the mean-flow/$\beta$ terms having
> no saturating mechanism, not a sign error. Adding bottom drag on the
> lower layer's relative vorticity (`r>0`, the standard ingredient in
> Phillips' original 1954 model) produces the expected life cycle —
> exponential growth, a peak, decay, and equilibration — confirmed by a
> multi-hundred-time-unit run at the notebook's actual resolution.

## Symmetric, inertial & Kelvin-Helmholtz instability (survey)
| Symbol | Meaning | Convention |
|---|---|---|
| $\mathrm{Ro}=-\dfrac{\partial U/\partial y}{f}$ | cross-front Rossby number | $f(1+\mathrm{Ro})$ = absolute vertical vorticity |
| $\mathrm{Ri}=N^2/(\partial U/\partial z)^2$ | gradient Richardson number | thermal-wind shear (symmetric); imposed shear (KH) |
| $J$ | Hazel (1972) profile's bulk/minimum Richardson number | $U=\tanh z$, $N^2=J\,\mathrm{sech}^2 z$ |

> **Symmetric & inertial instability: one Ertel-PV sign criterion.** For a
> zonal thermal-wind front, `gfdlib.symmetric` derives (not looks up) the
> Ertel PV $q=(f-\partial U/\partial y)N^2-f(\partial U/\partial z)^2$ and
> shows the front is unstable to slantwise (symmetric) displacements iff
> $qf<0$, i.e. $(1+\mathrm{Ro})<1/\mathrm{Ri}$. Two independently-checkable
> limits confirm the sign convention: at $\mathrm{Ro}=0$ the marginal curve
> reduces to the textbook value $\mathrm{Ri}=1$; for $(1+\mathrm{Ro})<0$
> (absolute vorticity changes sign) the front is unconditionally unstable
> at *every* $\mathrm{Ri}>0$ — pure inertial instability, no stratification
> dependence. `gfdlib.symmetric.classify` labels any $(\mathrm{Ro},
> \mathrm{Ri})$ point stable / symmetric / inertial / gravitational
> ($\mathrm{Ri}<0$) for the Ri-Ro stability map.
>
> **Kelvin-Helmholtz: the Taylor-Goldstein equation, and a real numerical
> pitfall caught along the way.** `gfdlib.instability.
> taylor_goldstein_growth_rate` solves $(U-c)^2(\phi''-k^2\phi)-U''(U-c)
> \phi+N^2\phi=0$ as a natively **linear** $2n\times2n$ generalized
> eigenvalue problem (keeping the buoyancy perturbation as an independent
> unknown rather than eliminating it via division by $(U-c)$, which gives a
> mathematically equivalent but quadratic-in-$c$ pencil). At $N^2=0$ this
> reduces exactly (to $10^{-13}$) to `gfdlib.instability.growth_rate` — a
> strong internal cross-check, not an external one. For the classic Hazel
> (1972) profile, growth was initially found NOT to vanish cleanly above the
> rigorous Miles-Howard threshold $\mathrm{Ri}\geq1/4$ — a hard theorem
> violation is impossible, so this was root-caused, not dismissed: the
> spurious growth scaled with grid spacing $\Delta z$ alone (matched
> $\Delta z$ across different $n$ and domain sizes gave matched spurious
> magnitude) and shrank as $\Delta z\to0$, confirming genuine 2nd-order
> truncation error — worst exactly at the delicate $\mathrm{Ri}=1/4$ point,
> where the equation's critical-layer indicial roots coalesce, a classically
> hard point for finite differences (professional codes use spectral
> methods here). Fixed by tightening the domain/resolution and applying a
> documented noise floor calibrated to this investigation
> (`_TG_NOISE_FLOOR`); the resulting curve decreases smoothly and reaches
> zero by $J\approx0.4$, correctly reproducing the theorem's qualitative
> Ri$\geq1/4$ cutoff (not a knife-edge at exactly $1/4$, an honest
> limitation of a simple solver on this classically delicate problem).
> `gfdlib.instability.hazel_profile` builds $U,N^2$ for this test case.

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

## Buoyancy-driven overturning (abyssal recipes; horizontal convection)
| Symbol | Meaning | Convention |
|---|---|---|
| $w$ | upwelling velocity (Munk's abyssal-recipe balance) | positive = upward |
| $\kappa$ | diapycnal diffusivity | boundary-layer thickness $\sim\kappa/w$ |
| $Q(x,z)$ | differential surface heating | $\propto\cos(2\pi x/L_x)$, warm at $x=0$, cool at $x=L_x/2$ |

> **Munk's (1966) abyssal recipe** $w\,\partial T/\partial z=\kappa\,
> \partial^2T/\partial z^2$ is solved exactly (`gfdlib.overturning.
> abyssal_profile`): a linear 2nd-order ODE with constant coefficients,
> general solution $T=A+Be^{wz/\kappa}$, with $A,B$ fixed by the two BCs
> via a $2\times2$ `numpy.linalg.solve`. Large $w/\kappa$ sweeps the
> bottom value through most of the water column, compressing the
> transition to the surface value into a thin layer near $z=H$ — the real
> ocean's near-uniform abyssal temperature capped by a thin thermocline;
> verified numerically, not just plotted, against both this limit and the
> $w\to0$ pure-diffusion (linear) profile.
>
> **The overturning cell reuses ch14's Boussinesq machinery unchanged.**
> `gfdlib.overturning.rhs_overturning` has the exact same structure as
> `convection.rhs_boussinesq`, but replaces the fixed background-
> stratification advection term ($+\psi_x$, uniform bottom heating) with
> an explicit, horizontally-varying heating field $Q(x,z)$
> (`gfdlib.overturning.surface_heating`) — warm at one end of the
> (periodic) domain, cool at the other, concentrated near the surface.
> This is Rossby's (1965) "horizontal convection" mechanism: differential
> heating along a boundary organizes into a single overturning cell
> (rising at the heated column, sinking at the cooled one) rather than
> the many small convective rolls of uniform bottom heating — verified
> directly, not just asserted: the vertical velocity $w=\psi_x$ comes out
> positive under the heated column and negative under the cooled one in
> every run tested.

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

## Eddy transport & mixing (passive tracer in 2D turbulence)
| Symbol | Meaning | Convention |
|---|---|---|
| $\Gamma=d\bar C/dy$ | imposed mean tracer gradient | $C_{total}=\Gamma y+c'(x,y,t)$, $c'$ periodic |
| $K_{eff}=-\overline{v'c'}/\Gamma$ | effective (eddy) diffusivity | domain average; $K_{eff}>0$ = down-gradient |

> **The mean-gradient trick.** A passive tracer with an unbounded mean
> gradient is handled the same way ch16/ch18 handle a mean shear or
> $\beta y$: split into an imposed linear part and a periodic perturbation
> $c'$, so `gfdlib.mixing.rhs_coupled` can evolve $c'$ on the same doubly-
> periodic `spectral.Grid` as the vorticity that stirs it, with the mean
> gradient appearing as an ordinary source term $-\Gamma v$ (the tracer is
> passive: it's advected by, but does not feed back on, the vorticity).
>
> **$K_{eff}$ is intrinsically $\Gamma$-independent — verified, not
> assumed.** The $c'$ equation is *linear* in $c'$, so $c'\propto\Gamma$ for
> a fixed flow realization and $K_{eff}=-\overline{v'c'}/\Gamma$ must be
> exactly independent of $\Gamma$. Checked numerically: two runs on the
> identical vorticity history with $\Gamma$ differing by a factor of 2 give
> $K_{eff}$ agreeing to machine precision (relative difference $10^{-16}$),
> and a control run with the velocity field identically zero gives
> $K_{eff}=0$ exactly (no possible eddy flux without a flow). For freely-
> decaying 2D turbulence, $K_{eff}$ comes out positive (down-gradient) and
> of the same order as a mixing-length estimate $u_{rms}\times$(domain
> scale) — a real emergent result of the simulation, not tuned to match.

## Wind-driven circulation (Stommel; Munk)
| Symbol | Meaning | Convention |
|---|---|---|
| $\varepsilon$ | Stommel nondimensional bottom-drag parameter | boundary layer width $\sim\varepsilon$ |
| $\delta$ | Munk nondimensional lateral-friction parameter | boundary layer width $\sim\delta$ |
| $f(x)$ | west-east profile, $\psi=f(x)\sin(\pi y)$ | $f(0)=0$ always; $f'(0)=0$ too for Munk (no-slip) |

> **A genuinely new domain shape.** Every earlier chapter's domain was
> doubly-periodic or a channel with two boundaries; this is the first
> *bounded rectangular basin*, four walls. `gfdlib.circulation` solves the
> Stommel and Munk gyre problems two independent ways: an exact analytic
> solution (separating $\psi=f(x)\sin(\pi y)$ reduces the PDE to a linear
> ODE for $f(x)$ with constant coefficients, solved via its characteristic
> polynomial's roots — quadratic formula for Stommel, `numpy.roots` for
> Munk's quartic — and a small `numpy.linalg.solve` for the boundary
> constants), and a general 2D finite-difference direct solve (Kronecker-
> sum operators, `numpy.linalg.solve`, no SciPy) that works for arbitrary
> forcing, cross-validated against the analytic solution with clean 2nd-order
> convergence (`gfdlib` tests).
>
> **Munk's biharmonic operator is handled by keeping $\zeta=\nabla^2\psi$ as
> an independent field** (the same pattern as ch16/17/19's coupled states)
> rather than a hand-derived 4th-order stencil, closed with Thom's (1933)
> wall-vorticity formula $\zeta_{wall}=2\psi_{adjacent}/d^2$ (derived by
> Taylor-expanding $\psi$ from a $\psi=\psi_n=0$ wall) for no-slip on the
> east/west walls, with free-slip ($\zeta=0$, no correction needed) on the
> north/south walls — the standard textbook simplification that keeps the
> problem separable.
>
> **A real bug, caught by a convergence test, not an eyeball check.** An
> earlier version enforced no-slip on *all four* walls (matching physical
> intuition, but not the separable-solution assumption). It disagreed with
> the exact analytic solution by a stubborn $\sim30\%$ — but critically,
> that error did **not** shrink with grid refinement at *any* boundary-layer
> width tested (0.08 to 0.4), immediately distinguishing a genuine
> discretization error (which must vanish as $\Delta x\to0$) from a
> boundary-condition mismatch (which does not). The fix — free-slip on
> north/south, matching the classical Pedlosky/Vallis treatment — restored
> clean 2nd-order convergence.

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
