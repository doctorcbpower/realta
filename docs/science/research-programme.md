# Realta Science Programme

Realta is intended to become a modular population-synthesis sandbox for exploring how assumptions about stellar populations, binaries, compact objects and their environments propagate into observable and thermodynamic consequences.

The central scientific theme is:

> **How does uncertain stellar and binary physics propagate from individual stars to stellar populations, galaxies and their surrounding environments?**

The initial programme consists of three linked papers plus one open investigation:

1. **X-ray fingerprints of massive-star multiplicity**
2. **Stochasticity and the apparent universality of X-ray/SFR relations**
3. **Massive binaries and X-ray feedback in the early Universe**
4. **Investigation: compact-object accretion and heating in galaxy groups and clusters**

Papers 1-3 should be developed sequentially, with the underlying Realta machinery designed so that the same population, evolution and observable interfaces can support all three.

Item 4 is deliberately kept as an investigation rather than a committed paper. Preliminary (unpublished, student) work suggests the compact-object accretion signal is likely small under standard Bondi-Hoyle assumptions in typical group/cluster gas -- see the "Investigation" section below for the order-of-magnitude argument. It is retained in the programme because the conclusion deserves a more careful second look (in particular, whether special environments such as cool cores or ram-pressure-stripped tails change the picture), not because the result is assumed. Whether it is promoted to a full paper depends on that revisit.

---

# Paper 1 — X-ray fingerprints of massive-star multiplicity

## Working title

**The X-ray fingerprints of massive-star multiplicity**

Alternative:

**From massive binaries to X-ray populations: disentangling the IMF and binary evolution**

## Scientific motivation

Young stellar populations provide an unusually clean laboratory for studying massive-star evolution.

The integrated emission of a young population depends on:

- the IMF;
- metallicity;
- stellar lifetimes;
- binary fraction;
- binary mass ratios;
- orbital-period distribution;
- mass transfer;
- common-envelope evolution;
- mergers;
- compact-object formation;
- supernova kicks.

The optical and UV luminosity of a young population primarily trace the massive-star population itself. X-ray emission, however, is particularly sensitive to the subsequent evolution of binaries and compact objects.

This raises a potentially useful question:

> **Can the time-dependent relationship between stellar/UV and X-ray luminosity distinguish the IMF from massive-star multiplicity and binary evolution?**

The key observable should therefore not simply be `L_X`, but ratios and trajectories such as

\[
R_X(t) = \frac{L_X(t)}{L_{\rm UV}(t)}
\]

and

\[
R_Q(t) = \frac{L_X(t)}{Q_{\rm H}(t)}.
\]

The temporal evolution may contain information that is lost in a single integrated luminosity.

---

## Basic experiment

Construct instantaneous-burst stellar populations with controlled:

- total stellar mass;
- metallicity;
- IMF;
- binary fraction;
- mass-ratio distribution;
- period distribution;
- binary-interaction prescription;
- merger prescription.

Initially compare:

1. single-star populations;
2. non-interacting binaries;
3. standard binary interaction;
4. enhanced interaction;
5. enhanced massive-star mergers.

Keep the IMF fixed while changing binary physics.

Then repeat while varying the IMF.

The aim is to determine whether the effects are observationally degenerate.

---

## Primary questions

### Q1. How rapidly does the X-ray population emerge?

Determine

\[
t_{\rm first\,X}
\]

and the evolution of

\[
L_X(t).
\]

### Q2. Does binary interaction leave a distinctive temporal signature?

Compare

\[
\frac{L_X}{L_{\rm UV}}
\]

between interaction models.

### Q3. Are massive-star mergers especially important?

Test whether mergers produce:

- delayed or enhanced UV emission;
- altered HMXB formation;
- altered compact-object populations;
- distinct X-ray/UV trajectories.

### Q4. Which observables best distinguish IMF and binary physics?

Compare the information contained in:

\[
L_{\rm bol},\quad
L_{\rm UV},\quad
Q_{\rm H},\quad
L_X.
\]

---

## Essential figures

### Figure 1 — Population evolution

For several binary models:

- `L_bol(t)`
- `L_UV(t)`
- `Q_H(t)`
- `L_X(t)`

This establishes the basic evolutionary behaviour.

### Figure 2 — X-ray/UV evolution

Plot

\[
L_X/L_{\rm UV}
\]

against age for different binary prescriptions.

This should be the central figure.

### Figure 3 — Effect of mergers

Compare:

- no mergers;
- standard mergers;
- enhanced mergers.

Show both luminosity evolution and compact-object formation.

### Figure 4 — IMF versus binary degeneracy

Construct a grid in:

\[
(\alpha_{\rm IMF},f_{\rm bin})
\]

and show an observable such as

\[
L_X/L_{\rm UV}
\]

at selected ages.

The goal is to identify degeneracies.

### Figure 5 — Metallicity

Repeat the principal experiment for several metallicities.

Show whether the binary signature survives changes in `Z`.

### Figure 6 — Stochastic realisations

For finite cluster masses, show distributions rather than means.

For example:

\[
P(L_X/L_{\rm UV}|M_{\rm cl},Z,t).
\]

This provides the bridge to Paper 2.

---

## Minimum viable paper

The first version does **not** require a complete MESA integration.

A useful first paper can be based on:

- analytic or tabulated stellar evolution;
- configurable binary prescriptions;
- simplified HMXB prescriptions;
- Monte Carlo population sampling.

MESA/FSPS/POSYDON can subsequently be used as validation or higher-fidelity physics providers.

---

# Paper 2 — Stochasticity and the X-ray/SFR relation

## Working title

**When is X-ray luminosity a reliable tracer of star formation?**

Alternative:

**The stochastic X-ray output of young stellar populations**

## Scientific motivation

Integrated galaxy observables are often interpreted using relations such as

\[
L_X \propto {\rm SFR}.
\]

However, X-ray emission can be dominated by a relatively small number of rare sources.

At low stellar masses, the population may not adequately sample:

- the upper IMF;
- massive binaries;
- compact-object formation;
- HMXB evolutionary phases.

Consequently,

\[
\langle L_X\rangle
\]

may be a poor description of an individual galaxy.

The important question is:

> **How much of the scatter in X-ray/SFR relations arises from IMF sampling, and how much arises from stochastic binary evolution?**

---

## Basic experiment

Generate many Monte Carlo populations spanning:

\[
M_{\rm cl}=10^2-10^6\,M_\odot
\]

and eventually galaxy populations spanning a larger mass range.

For each realisation calculate:

\[
L_X,\quad
L_{\rm UV},\quad
Q_{\rm H}.
\]

Measure:

\[
P(L_X|M,Z,t)
\]

and

\[
P(L_X/{\rm SFR}|M,Z).
\]

---

## Primary questions

### Q1. At what mass does X-ray emission become statistically stable?

Determine the transition between:

- individual-source dominated;
- stochastic population;
- effectively continuous population.

### Q2. How important is binary sampling?

Separate:

- IMF sampling;
- binary sampling;
- evolutionary stochasticity.

### Q3. How does metallicity affect the scatter?

Determine whether low-metallicity populations have both:

- enhanced mean `L_X`;
- enhanced stochasticity.

### Q4. Can UV and X-ray observations distinguish the causes of scatter?

Compare:

\[
L_X/L_{\rm UV}
\]

rather than only `L_X/SFR`.

---

## Essential figures

### Figure 1 — Distribution of `L_X`

Histograms or PDFs of `L_X` for different cluster masses.

### Figure 2 — Scatter versus population mass

Plot

\[
\sigma[\log L_X]
\]

against stellar population mass.

### Figure 3 — IMF sampling versus binary sampling

Run controlled experiments:

- stochastic IMF + fixed binary model;
- fixed IMF + stochastic binaries;
- both stochastic.

This is likely to be one of the most useful figures.

### Figure 4 — `L_X/SFR` distribution

Show the probability distribution rather than a single relation.

### Figure 5 — Metallicity dependence

Mean and scatter of

\[
L_X/{\rm SFR}
\]

as a function of metallicity.

### Figure 6 — Cluster populations forming galaxies

Construct a galaxy from many star clusters with an imposed SFH.

Show convergence from:

\[
{\rm cluster}\rightarrow{\rm galaxy}.
\]

---

## Key conceptual result

The paper should establish that the commonly used population-level X-ray/SFR relation is not necessarily a fundamental relation.

Instead:

\[
L_X =
L_X({\rm IMF},f_{\rm bin},Z,{\rm SFH},M_\star,\mathrm{stochasticity}).
\]

This provides the formal bridge to high-redshift galaxies.

---

# Paper 3 — Massive binaries and X-ray feedback in the early Universe

## Working title

**From massive binaries to the early-Universe X-ray background**

Alternative:

**Binary evolution and the X-ray heating of the early Universe**

## Scientific motivation

High-redshift galaxies contain young, metal-poor stellar populations.

Their X-ray output depends on:

\[
{\rm SFR}(t),
\quad
Z(t),
\quad
{\rm IMF},
\quad
f_{\rm bin},
\quad
\mathrm{binary\ evolution}.
\]

The uncertainty in early-Universe X-ray heating therefore partly originates in unresolved stellar physics.

The key question is:

> **How much uncertainty in early-Universe X-ray feedback is generated by uncertainty in massive-star and binary evolution?**

---

## Core mathematical framework

Construct a population response function:

\[
\Psi_X(\tau,Z,\Theta_{\rm bin})
\]

where `tau` is time since star formation and `Theta_bin` contains the binary parameters.

Then:

\[
L_X(t)
=
\int_0^t
{\rm SFR}(t-\tau)
\Psi_X(\tau,Z(t-\tau),\Theta_{\rm bin})
\,d\tau.
\]

This makes Realta a natural stellar-population kernel generator.

---

## Primary questions

### Q1. How does binary physics alter `L_X/SFR` at high redshift?

### Q2. How much does metallicity amplify the effect?

### Q3. Does uncertainty in binary evolution dominate uncertainty in IMF assumptions?

### Q4. What range of X-ray heating histories follows?

### Q5. Which stellar observables could constrain the relevant uncertainty?

This last question provides a route back to observations.

---

## Essential figures

### Figure 1 — Response functions

\[
\Psi_X(\tau,Z)
\]

for different binary prescriptions.

### Figure 2 — High-redshift `L_X/SFR`

Compare binary models over cosmic time.

### Figure 3 — IMF versus binary uncertainty

Show the relative effect of varying:

- IMF;
- binary fraction;
- merger fraction;
- metallicity.

### Figure 4 — X-ray spectral energy distribution

Where possible, show how the assumed compact-object population changes the spectrum.

### Figure 5 — Heating history

Calculate the resulting X-ray energy deposition/heating history.

### Figure 6 — Parameter uncertainty

Show the envelope of plausible histories arising from stellar/binary uncertainties.

---

# Investigation — Compact-object accretion in galaxy groups and clusters

> **Status: investigation, not a committed paper.** Preliminary (unpublished, student) work suggests this signal is likely negligible under standard assumptions -- consistent with the order-of-magnitude estimate below. Retained here to be revisited more carefully (alternative accretion prescriptions, special environments such as cool cores or stripped tails) before deciding whether it becomes a paper. See docs/science/development-roadmap.md Phase 5 for its place in the development sequence.

## Working title

**Can stellar remnants heat the intragroup and intracluster medium?**

Alternative:

**The contribution of stellar-remnant accretion to intracluster-medium heating**

## Scientific motivation

Galaxies continuously produce:

- white dwarfs;
- neutron stars;
- stellar-mass black holes.

Some compact objects are retained by galaxies, while others can be displaced by:

- supernova kicks;
- dynamical interactions;
- galaxy disruption;
- tidal stripping.

Galaxy groups and clusters also contain an intracluster stellar population.

This creates a population of compact objects embedded in, or moving through, hot diffuse gas.

The question is:

> **Can accretion onto stellar remnants provide a significant source of energy to the intragroup/intracluster medium over Gyr timescales?**

The important comparison is not simply luminosity.

It is:

\[
E_{\rm CO}
=
\int L_{\rm CO}(t)\,dt
\]

against the thermal and radiative energy scales of the surrounding gas.

---

## Model ingredients

### Stellar evolution

Predict:

\[
N_{\rm WD},N_{\rm NS},N_{\rm BH}
\]

as a function of time.

### Dynamical redistribution

Track the fraction of remnants that are:

- retained in galaxies;
- ejected;
- intracluster.

### Gas environment

Specify:

\[
\rho_{\rm gas}(r,t),
\quad
T(r,t),
\quad
v_{\rm rel}(r,t).
\]

### Accretion

Initially use a simple prescription such as Bondi-Hoyle:

\[
\dot M_{\rm B}
=
\frac{4\pi G^2M^2\rho}
{(c_s^2+v_{\rm rel}^2)^{3/2}}.
\]

Then:

\[
L_{\rm acc}=\eta\dot M c^2.
\]

The efficiency should be configurable rather than hard-coded.

---

## Primary questions

### Q1. What is the compact-object population in a group/cluster?

Predict:

\[
N_{\rm BH}(t),N_{\rm NS}(t),N_{\rm WD}(t).
\]

### Q2. What fraction becomes intracluster?

Explore supernova kicks and dynamical ejection.

### Q3. What is the integrated accretion luminosity?

Calculate:

\[
L_{\rm CO}(t).
\]

### Q4. Where does the energy go?

Compare:

\[
L_{\rm CO}
\]

with:

- ICM cooling;
- thermal energy;
- AGN heating;
- stellar feedback.

### Q5. Is the effect more interesting in groups than clusters?

This is particularly important.

The lower temperatures and densities of groups may make the comparison with cooling different from that in rich clusters.

---

## Essential figures

### Figure 1 — Compact-object census

Evolution of:

\[
N_{\rm BH},N_{\rm NS},N_{\rm WD}.
\]

### Figure 2 — Spatial distribution

Radial distribution of:

- galactic remnants;
- ejected remnants;
- intracluster remnants.

### Figure 3 — Accretion luminosity

\[
L_{\rm CO}(t)
\]

over several Gyr.

### Figure 4 — Energy budget

Cumulative:

\[
E_{\rm CO}(t)
\]

compared with:

\[
E_{\rm cool}(t)
\]

and other relevant energy scales.

### Figure 5 — Group versus cluster

Compare several halo masses.

### Figure 6 — Sensitivity analysis

Show dependence on:

- BH mass function;
- kick velocity;
- retention fraction;
- accretion efficiency;
- gas density;
- relative velocity.

### Figure 7 — Spatial heating profile

Calculate something like

\[
\epsilon_{\rm CO}(r)
=
n_{\rm CO}(r)L_{\rm acc}(r).
\]

This determines whether the heating is centrally concentrated or diffuse.

---

## Important caution

The first version should **not** claim to model ICM feedback self-consistently.

Instead:

> Calculate the maximum plausible contribution of compact-object accretion under a controlled set of assumptions.

If the result is negligible, that is still useful.

If it is non-negligible, the investigation can be promoted to a paper, and the next stage can introduce more realistic gas dynamics.

### Why this needs revisiting: an order-of-magnitude check

A quick Bondi-Hoyle estimate for a single 10 M_sun BH in typical group/cluster gas (n ~ 10^-3 cm^-3, kT ~ 1 keV so c_s ~ 500 km/s) gives

\[
\dot M_{\rm B} \sim \frac{4\pi G^2 M^2 \rho}{c_s^3} \sim 10^{-20}\,M_\odot\,{\rm yr}^{-1},
\]

roughly 10^-14 of the Eddington rate for that mass. Even summed over 10^5-10^7 such remnants in a massive cluster, the resulting \(L_{\rm CO}\) is many orders of magnitude below typical ICM cooling luminosities (~10^43-10^44 erg/s). This is consistent with the preliminary student result and is the reason this project is an investigation rather than a committed paper: the diffuse, typical-density case is very likely a null result. The open question is whether the interesting regime is instead confined to special environments -- cool cores, ram-pressure-stripped tails, dense cold clumps -- where densities are orders of magnitude higher and lower relative velocities apply. Figure 6 (sensitivity analysis) is where that question should actually be answered before any claim is made either way.

---

# Relationship between the papers

Papers 1-3 should deliberately build on one another. The compact-object investigation (Item 4) draws on the same machinery but is not part of this sequential chain -- it can be pursued in parallel once Paper 1's compact-object/HMXB formation code exists, and only joins the numbered sequence (as a future paper) if it is promoted.

```text
Paper 1
Massive-star multiplicity
        |
        v
Paper 2
Population stochasticity
        |
        v
Paper 3
Early-Universe X-ray feedback


Investigation (parallel, not yet sequenced)
Compact-object accretion in groups/clusters
```

The same Realta machinery should underpin all three papers, and the investigation reuses the compact-object and accretion abstractions rather than duplicating them.
The scientific progression is:

```text
stellar physics
      |
      v
binary evolution
      |
      v
compact-object population
      |
      v
radiation
      |
      v
galaxy environment
      |
      v
IGM / ICM
```

This is the central narrative of the programme.

