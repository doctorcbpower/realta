# Wind-capture accretion

Not from Power et al. A physically-motivated wind-capture model for
persistent, wind-fed HMXBs (the Vela X-1-type regime, not requiring
Roche-lobe overflow): the secondary's own radiatively-driven (CAK)
wind, captured by the compact primary via Bondi-Hoyle-Lyttleton
accretion. Based on El Mellah & Casse (2017, MNRAS, arXiv:1609.01532),
with the wind-driving physics informed by Friend & Castor (1982, ApJ
261, 293).

Module split: `stellar/cak_wind.py` computes the donor's wind state
(mass-loss rate, velocity law); `binaries/wind_capture.py` converts
that into a capture fraction, accretion rate, and circularization
radius.

## CAK wind

`Mdot_wind`, terminal velocity, velocity profile `v(r)` — El Mellah &
Casse (2017) eq. (20) (Gayley 1995's Q-parametrization) for the
mass-loss rate, and eq. (7) for the terminal velocity/profile:

```
v_inf = 2.5 * v_esc * alpha/(1-alpha)
v(r)  = v_inf * (1 - R*/r)^0.7
```

Eq. (7) is their own calibrated fit, valid for T > 21,000 K (every
realistic OB-supergiant HMXB donor); NOT their more basic point-source
form (eq. 5, `v_inf = v_esc*sqrt(alpha/(1-alpha))`), which undershoots
observed terminal velocities for typical OB-supergiant parameters.

Calibrated against Friend & Castor's own Vela X-1 donor data (their
Table 1: M=24, R=35 Rsun, L=5e5 Lsun, Γ*=0.50, observed v_inf =
700-1700 km/s, Mdot = 0.6-2e-6 Msun/yr): `eddington_factor` recovers
Γ* = 0.50 to ~8%; `wind_terminal_velocity` (eq. 7) lands inside the
observed range; `wind_mass_loss_rate` (Q ≈ 900, El Mellah & Casse's
OB-supergiant fiducial) is within ~2 orders of magnitude of the
observed rate — both papers themselves flag `Q`/`alpha` as needing
per-star calibration, not a generic value.

Implementation: `stellar/cak_wind.py::wind_mass_loss_rate`,
`wind_terminal_velocity`, `wind_velocity`, `eddington_factor`.
Tests: `tests/test_cak_wind.py`.

## Wind capture

Accretion rate onto the compact object — El Mellah & Casse eq. (18),
their own recommended simplified estimate:

```
beta = 0.77 * (R_acc/a)^2
```

evaluated using the wind velocity AT THE ORBITAL SEPARATION (recovers
their full numerical result to within 6%, per their own text). NOT
their eq. (19), which their own text states underestimates by "at
least a factor of three".

`R_acc = 2*G*M_compact/v_rel^2` (Hoyle & Lyttleton 1939; Bondi & Hoyle
1944); `v_rel` combines the wind and orbital velocity in quadrature.

A plain Bondi-Hoyle-Lyttleton fallback using a caller-supplied wind
speed (e.g. terminal velocity, without the orbital-separation
refinement) is also provided for comparison.

Implementation: `binaries/wind_capture.py::bhl_accretion_fraction`,
`wind_capture_rate`, `bhl_accretion_rate_simple`, `accretion_radius`,
`relative_wind_velocity`.
Tests: `tests/test_wind_capture.py`.

## Circularization radius

Not from either source paper — El Mellah & Casse compute this
numerically only (their Fig. 6), giving no closed-form fit. Uses the
standard wind-accretion literature scaling (Shapiro & Lightman 1976):

```
R_circ/R_acc ≈ (1/4) * (v_orbital/v_rel)^4
```

The `(v_orbital/v_rel)^4` scaling is the standard form in the
wind-accretion literature; the 1/4 prefactor is the commonly-cited
value but not independently verified against a source paper — the
lowest-confidence piece of this module.

Implementation: `binaries/wind_capture.py::circularization_radius_fraction`,
`circularization_radius`.
Tests: `tests/test_wind_capture.py`.

## Wiring into `evolve()`

Opt-in: `config.use_wind_capture`, default `False`. Shares Phase 1.5
with [`post-sn-rlof.md`](post-sn-rlof.md) (`nturn == 1`), applying
while the secondary donor has NOT yet filled its Roche lobe
(`donor_radius < r_l2`) — `use_post_sn_rlof` takes over once it does.
Either flag works independently; together they compose into one
continuous `nturn == 1` pipeline (wind-fed while detached, RLOF once
the donor overflows).

`L_X = eta * Mdot_acc * c^2`, `eta = 0.1`
(`xray/luminosity.py::XRayLuminosity.eta`, the same efficiency
already used for that class's Eddington-luminosity cap), Eddington-
capped, stored `lunit`-normalized like every other `lum_xray` value.

CAK shape parameters: `config.wind_cak_alpha = 0.55` (midpoint of El
Mellah & Casse's cited 0.45-0.65 OB-supergiant range),
`config.wind_cak_q_force = 900.0` (their stated typical value). Both
validated in `__post_init__`.

The wind velocity at the orbital separation
(`cak_wind.wind_velocity(a, donor_radius, v_inf)`) falls as
`donor_radius` grows toward the separation, which raises the capture
fraction (`bhl_accretion_fraction`, `~1/v_rel^2`) continuously right
up to the discrete hand-off to `use_post_sn_rlof` — an emergent
consequence of the existing formulae, not a separate enhancement term.

Activated once and held fixed for the rest of the binary's active
lifetime (matching `f_sur`/`use_post_sn_rlof`'s own convention), not
recomputed every timestep as the donor evolves.

Implementation: `binaries/population.py::evolve`, Phase 1.5.
Tests: `tests/test_wind_capture_wiring.py`.
