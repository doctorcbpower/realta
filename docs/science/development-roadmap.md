# Realta Development Roadmap

## Purpose

Realta should be developed as a modular population-synthesis sandbox rather than as a monolithic stellar-evolution code.

The design goal is:

> **Make it easy to change a physical assumption, generate a population, evolve it, inspect the resulting events and observables, and reproduce the experiment.**

The development roadmap is therefore driven by the science programme.

---

# 1. Core design principle

A Realta experiment should eventually look approximately like:

```python
population = Population(
    imf=Kroupa(),
    metallicity=0.002,
    binary_model=BinaryModel(...),
    stellar_model=StellarModel(...),
)

history = population.evolve(times)

history.plot("xray")
history.plot("uv")
history.plot("xray_to_uv")

population.compact_objects()
population.events()
```

The user should not need to understand the internal implementation to run an experiment.
2. Development priorities
Priority 1 — Population as the central abstraction
Implement a clean population object containing:

* stars;
* binaries;
* compact objects;
* events;
* model metadata;
* random-number state.

The population should be capable of evolving independently of how individual stellar models are implemented.
Required API

```python
population = Population(...)

population.evolve(times)

population.stars
population.binaries
population.compact_objects
population.events
population.history
```

3. Time-history abstraction
Everything in the science programme depends on time-dependent quantities.
Introduce a standard `PopulationHistory` object.
It should contain:

```text
time
stellar_mass
living_mass
remnant_mass

L_bol
L_UV
L_optical
L_X
Q_H

N_WD
N_NS
N_BH

mergers
SNe
XRBs
ULXs
```

The exact quantities can expand over time.
The important principle is that all models return quantities on a common time grid.
4. Events as first-class objects
Binary evolution naturally produces discrete events.
Introduce an event abstraction:

```python
Event(
    time=...,
    type="merger",
    objects=(...),
    properties={...}
)
```

Possible event types:

```text
birth
mass_transfer
RLOF
common_envelope
merger
SN
BH_formation
NS_formation
WD_formation
XRB_on
XRB_off
ULX_on
ejection
```

This will make later analysis substantially easier.
5. Separate physics from population sampling
Do not mix:

```text
"how a binary evolves"
```

with:

```text
"how binaries are sampled".
```

For example:

```python
binary_population = BinaryPopulation(
    fraction=0.7,
    mass_ratio_distribution=...,
    period_distribution=...,
)
```

should be independent of:

```python
binary_model = BinaryEvolutionModel(...)
```

This separation is essential for Paper 1.
6. Configurable IMFs
Implement a common IMF interface:

```python
imf = Kroupa(...)
imf = Salpeter(...)
imf = Chabrier(...)
```

Required operations:

```python
imf.pdf(mass)
imf.cdf(mass)
imf.sample(n)
imf.integral(...)
```

The IMF should be serialisable into the experiment metadata.
7. Configurable binary distributions
Implement independent distributions for:

```text
binary fraction
primary mass
mass ratio
period
eccentricity
```

Example:

```yaml
binary:
  fraction: 0.7

  mass_ratio:
    model: flat

  period:
    model: log_normal

  eccentricity:
    model: thermal
```

The goal is to make Paper 1 experiments configuration changes rather than code changes.
8. Stellar evolution interface
Create a common interface for stellar evolution providers.
For example:

```python
class StellarModel:
    def evolve(self, mass, metallicity, times):
        ...
```

Potential implementations:

```text
AnalyticStellarModel
InterpolatedStellarModel
MESAStellarModel
FSPSStellarModel
```

Realta should not require MESA or FSPS for basic operation.
They should be optional providers.
9. FSPS integration
FSPS should initially be used to provide population-level spectral information rather than to control the entire Realta evolution.
Potential interface:

```python
spectra = FSPSModel(
    metallicity=...,
    age=...
)

spectra.luminosity("FUV")
spectra.luminosity("NUV")
spectra.luminosity("optical")
```

The interface should allow Realta to ask for:

```text
L_bol
L_UV
ionising photon rate
SED
```

This is primarily required for Papers 1, 2 and 4.
10. Binary evolution interface
Introduce a provider interface:

```python
class BinaryEvolutionModel:
    def evolve(self, binary, times):
        ...
```

Possible providers:

```text
SimpleBinaryModel
MESA/BinaryModel
POSYDONModel
```

The first model can be intentionally simple.
The important requirement is that binary evolution returns standardised outputs:

```text
masses
stellar types
separation
period
mass transfer
merger state
remnant state
events
```

11. Compact-object abstraction
Create a common compact-object class:

```python
CompactObject(
    type="BH",
    mass=10.0,
    natal_kick=...,
    formation_time=...,
    position=...,
    velocity=...
)
```

Types:

```text
WD
NS
BH
```

Properties should include:

```text
formation time
mass
spin (optional)
kick
position
velocity
host galaxy
host cluster
```

This is essential for the compact-object accretion investigation (see docs/science/research-programme.md).
12. Accretion as a modular physics provider
Do not hard-code Bondi accretion.
Create:

```python
class AccretionModel:
    def mdot(self, compact_object, environment):
        ...
```

Initial implementation:

```text
BondiHoyleAccretion
```

Later:

```text
EddingtonLimitedAccretion
RadiativelyInefficientAccretion
DiskAccretion
```

The output should include:

```text
mdot
luminosity
radiative_efficiency
mechanical_power
```

This will make the compact-object accretion investigation much easier to extend.
13. Environment abstraction
Introduce an environment interface:

```python
environment = GasEnvironment(
    density=...,
    temperature=...,
    velocity=...
)
```

It should support:

```python
environment.density(position, time)
environment.temperature(position, time)
environment.velocity(position, time)
```

Initially this can be analytic.
Examples:

```text
UniformICM
BetaModelICM
NFWGasModel
RadialProfileICM
```

This avoids embedding cluster physics inside the compact-object code.
14. Dynamical populations
For the compact-object accretion investigation, Realta eventually needs to distinguish:

```text
objects bound to galaxies
objects ejected from galaxies
intracluster objects
```

Introduce a simple environment hierarchy:

```text
Universe
 └── Cluster
      ├── Galaxy
      │    ├── Stars
      │    └── CompactObjects
      └── IntraclusterPopulation
```

The first implementation can be semi-analytic.
No N-body simulation is required.
15. Population response functions
This is one of the most important developments.
Realta should be able to produce:

```python
response = population.response(
    observable="xray",
    metallicity=0.002
)
```

giving something conceptually like:
[
\Psi_X(\tau,Z,\Theta).
]
This response function becomes the interface between Realta and galaxy-formation models.
It allows:

```text
stellar population synthesis
          |
          v
response kernel
          |
          v
galaxy SFH
          |
          v
galaxy luminosity
```

This is the key development for Paper 3 (early-Universe X-ray feedback). It is also directly useful for Paper 2: convolving a single-realisation response function against a constant SFH gives the quasi-continuous limit that "L_X/SFR" implicitly assumes, so a minimal version of this response-function machinery is worth pulling forward rather than waiting for Phase 5.
16. Convolution with galaxy SFHs
Implement:

```python
galaxy = GalaxyHistory(
    time=...,
    sfr=...,
    metallicity=...
)

luminosity = realta.convolve(
    galaxy,
    response
)
```

This should produce:

```text
L_bol(t)
L_UV(t)
Q_H(t)
L_X(t)
```

This will allow Realta to connect to:

* SHARK;
* semi-analytic galaxy models;
* cosmological simulations;
* observed SFHs.

17. Stochastic experiment framework
Papers 1 and 2 require many realisations.
Create an experiment runner:

```python
experiment = Experiment(
    model=population_model,
    parameters=parameter_grid,
    n_realizations=1000,
    seed=12345,
)

results = experiment.run()
```

Results should retain:

```text
parameter values
random seed
model version
summary observables
```

The random seed must always be recorded.
18. Parameter sweeps
Support:

```yaml
parameters:

  metallicity:
    values: [0.0001, 0.001, 0.002, 0.01]

  binary_fraction:
    values: [0.0, 0.3, 0.7, 1.0]

  merger_efficiency:
    values: [0.0, 0.1, 0.5, 1.0]
```

The output should naturally form a tidy table.
This makes parameter-space plots trivial.
19. Standard analysis products
Realta should provide reusable analysis functions.
For example:

```python
history.xray_to_uv()
history.xray_to_ionising()
history.sfr_indicator()
history.remnant_counts()
history.energy_budget()
```

And:

```python
results.mean(...)
results.median(...)
results.percentile(...)
results.distribution(...)
```

The goal is to avoid writing the same analysis code in every paper notebook.
20. Plotting API
Create standard plotting functions:

```python
history.plot("luminosities")
history.plot("xray_to_uv")
history.plot("compact_objects")
history.plot("energy_budget")
```

But keep plotting separate from the underlying science objects.
The plotting layer should accept standard tables/objects so that publication plots can later be customised.
21. Reproducible experiment files
Every science result should be generated from a configuration file.
For example:

```yaml
experiment:
  name: massive_binary_xray
  seed: 12345

population:
  mass: 1.0e5
  metallicity: 0.002

imf:
  model: kroupa

binary:
  fraction: 0.7
  merger_model: standard

stellar:
  model: analytic

xray:
  model: hmxrb_v1
```

The code should be able to reproduce the run from this file.
22. Experiment provenance
Every output should contain:

```text
Realta version
git commit
configuration
random seed
input data versions
physics model versions
```

This should be automatic.
This is particularly important once collaborators start running their own experiments.
23. Data output
Prefer standard tabular outputs.
For example:

```text
population.parquet
history.parquet
events.parquet
compact_objects.parquet
```

HDF5 can remain available for large datasets.
The important requirement is that outputs are:

* machine-readable;
* portable;
* self-describing;
* versioned.

24. Validation framework
Every physics module should have tests against either:

* analytic expectations;
* published results;
* known limiting cases;
* external codes.

Examples:
IMF
Verify sampled distributions reproduce the analytic IMF.
Stellar evolution
Verify lifetimes and remnant masses against reference values.
Binary evolution
Verify limiting cases:

```text
binary fraction = 0
```

reduces to the single-star population.
Accretion
Verify Bondi scaling:
[
\dot M\propto M^2\rho
(c_s^2+v^2)^{-3/2}.
]
Population synthesis
Verify convergence with increasing population size.
25. Paper-generation notebooks
Each paper should have a small number of notebooks whose only purpose is:

1. load configuration;
2. run or load experiment;
3. generate publication tables;
4. generate publication figures.

For example:

```text
papers/
  paper1_massive_binaries/
      README.md
      configs/
      notebooks/
      figures/
      tables/

  paper2_stochastic_xray/
      ...

  paper3_early_universe_xray/
      ...

investigations/
  compact_object_icm/
      README.md
      configs/
      notebooks/
      figures/
      tables/
```

The `investigations/` directory mirrors the `papers/` structure but signals explicitly that the compact-object accretion work is not yet a committed paper -- promoting it later is just a rename/move, not a restructure.

The notebooks should not contain the underlying physics.
They should call Realta.
26. Recommended development sequence
Phase 1 — Make the existing code a research instrument
Implement:

* clean `Population`;
* `PopulationHistory`;
* standard observables;
* configuration files;
* reproducible random seeds;
* basic plotting;
* experiment runner.

This is the immediate priority.
Phase 2 — Paper 1
Add:

* binary population distributions;
* binary evolution interface;
* merger events;
* compact-object formation;
* HMXB prescription;
* UV/ionising observables.

Goal:
Produce the first `L_X/L_UV` experiment entirely from Realta.
Phase 3 — Paper 2
Add:

* large Monte Carlo experiments;
* population convergence;
* stochastic sampling;
* parameter grids;
* distribution analysis.

Goal:
Make stochastic population synthesis a first-class Realta capability.
Phase 4 — Paper 3
Add:

* response functions;
* SFH convolution;
* metallicity histories;
* galaxy interfaces;
* X-ray SEDs;
* early-Universe heating calculations.

Goal:
Make Realta callable as a stellar/binary feedback module from a galaxy-formation model.
Phase 5 — Investigation: compact-object accretion (exploratory, not yet a committed paper)
Add, if and when revisited:

* spatial compact-object populations;
* natal kicks;
* ejection;
* intracluster populations;
* gas environments;
* accretion models;
* integrated energy budgets.

Goal:
Confirm or refute, with Realta's own machinery, whether compact-object accretion is a non-negligible group/cluster heating channel. Prior (unpublished, student) work suggests it is likely negligible under standard Bondi-Hoyle assumptions in typical gas -- see docs/science/research-programme.md's order-of-magnitude estimate. This phase exists to revisit that conclusion with more care (alternative accretion prescriptions, special environments) rather than to commit to a paper outright; promote it to "Paper 4" only if the revisit finds a genuine effect.
27. What Realta should NOT do
Realta should not become:

* a replacement for MESA;
* a replacement for FSPS;
* a replacement for POSYDON;
* a full N-body code;
* a hydrodynamic simulation;
* a galaxy-formation model.

Instead:

```text
MESA / POSYDON / FSPS
          |
          v
       Realta
          |
          v
 stellar populations
          |
          v
 observables / remnants
          |
          v
 galaxy / group / cluster models
```

Realta's value is the population-level experimental framework.
28. The long-term architecture
The eventual architecture should look roughly like:

```text
                         REALTA
                           |
             ┌─────────────┼─────────────┐
             |             |             |
           IMF          STARS         BINARIES
             |             |             |
             └─────────────┼─────────────┘
                           |
                    COMPACT OBJECTS
                           |
          ┌────────────────┼────────────────┐
          |                |                |
        XRBs             MERGERS        ACCRETION
          |                |                |
          └────────────────┼────────────────┘
                           |
                      OBSERVABLES
                           |
             ┌─────────────┼─────────────┐
             |             |             |
            UV            X-ray       remnants
             |             |             |
             └─────────────┼─────────────┘
                           |
                  ENVIRONMENT MODELS
                           |
          ┌────────────────┼────────────────┐
          |                |                |
        galaxy           group           cluster
          |                |                |
          └────────────────┼────────────────┘
                           |
                    COSMOLOGICAL USE
```

The guiding principle is:
Every new physical model should plug into the same population → history → observable pipeline.
This keeps Realta small enough to understand while allowing it to grow into a genuinely useful research platform.
