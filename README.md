# `realta`

[./realta_banner.jpg|Realta Banner in Celtic Script]

[![CI](https://github.com/doctorcbpower/realta/actions/workflows/ci.yml/badge.svg)](https://github.com/doctorcbpower/realta/actions/workflows/ci.yml)
[![Python Version](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**realta** is a modular stellar and binary population-synthesis framework. It provides a transparent and reproducible way to construct stellar populations from specified formation histories, assign and evolve binary systems under explicit physical prescriptions, and follow the resulting stellar remnants and observables.

The current implementation is based on the coeval stellar-population model developed for high-mass X-ray binary (HMXB) studies in globular clusters, following Power et al. (2009). This provides a well-defined and reproducible baseline from which the code can be developed toward more general stellar and binary populations.

The long-term aim is to make Realta a flexible research tool for exploring the connection between:

```text
stellar population formation
          ↓
initial stellar and binary populations
          ↓
stellar and binary evolution
          ↓
compact-object formation
          ↓
observables and environments
```

This includes the possibility of arbitrary star-formation histories, improved empirical or theoretical models for binary initial conditions, more detailed stellar evolution through interfaces such as MESA, and connections to Galactic dynamics through tools such as `galpy`.

These extensions are intended to remain modular. Realta is **not** intended to become a universal stellar-evolution, gravitational-wave or Galactic-dynamics code. Instead, it aims to provide a controlled population-synthesis layer that can interface with specialist tools where appropriate.

The Power et al. (2009) model remains the reference implementation and should remain reproducible as the code develops.

## Key Features

* **Monte Carlo Population Synthesis:** Constructs and evolves stellar and binary populations using explicitly specified physical prescriptions.

* **Coeval Populations:** The current implementation supports the coeval stellar populations used in the original globular-cluster HMXB studies.

* **Extensible Formation Histories:** The architecture is intended to support more general star-formation histories, including continuous and user-defined formation histories.

* **Binary Population Models:** Samples binary parameters including primary masses, companion masses and orbital properties. The current implementation follows the assumptions of the reference model, with scope for incorporating improved empirical and theoretical distributions.

* **Stellar Evolution:** Uses tabulated stellar lifetimes, remnant masses and ionising photon rates in the current implementation.

* **HMXB and X-ray Modelling:** Implements the HMXB and X-ray prescriptions used in the reference model.

* **Reproducible Simulations:** Uses explicit random-number streams and configurable random seeds.

* **Clean Python API:** Simulations can be run programmatically from Python as well as from the command line.

* **Modular Development:** Physical prescriptions are separated so that alternative models can be introduced without unnecessarily modifying the rest of the population-synthesis machinery.

---

## Scientific Scope

Realta is being developed around a simple principle:

> **Do one scientifically coherent thing well, while allowing different physical models and applications to be connected to the same population framework.**

The core problem is the construction and evolution of stellar and binary populations. Specific applications may include:

* HMXB populations;
* globular-cluster stellar populations;
* compact-object populations;
* alternative binary-evolution prescriptions;
* gravitational-wave progenitor populations;
* Galactic compact-object populations.

Not every application will be appropriate for every physical prescription. In particular, the current HMXB model was developed for massive stellar populations and should not automatically be applied to arbitrary or old stellar populations simply because Realta can generate those populations.

The intention is therefore to keep the **population framework** general while making the **domain of validity of individual physical prescriptions** explicit.

---

## Installation

### Prerequisites

* Python `>= 3.11`
* `pip` package manager

### 1. Basic Installation

Clone the repository and install the package locally:

```bash
git clone https://github.com/doctorcbpower/realta.git
cd realta
pip install .
```

### 2. Developer Installation

To install in editable mode with development dependencies (`pytest`, `ruff`, `mypy`):

```bash
pip install -e ".[dev]"
```

---

## Quickstart & Usage

### 1. Command Line Interface (CLI)

Once installed, a simulation can be run directly using the `realta` executable.
Simulation parameters are set via a YAML config file (see `config.yml` for the
full field list, units, and Power et al. 2009 provenance of each value) --
there are no per-parameter CLI flags:

```bash
# Run with a config file
realta --config config.yml --output ./output

# Or with defaults (see SimulationConfig in src/realta/config.py)
realta --output ./output

# -v/--verbose for debug-level logging
realta --config config.yml --output ./output --verbose
```

### 2. Python API

Realta can also be used as a Python module in analysis scripts:

```python
from realta import ClusterSimulation, SimulationConfig

config = SimulationConfig()

config.ntot = 5000       # Total number of initial binary systems
config.tmax = 12.0       # Maximum simulation time [Myr]
config.iseed = 42        # Random seed

cluster = ClusterSimulation(config)

results = cluster.run(output_dir="data/output")

print("Simulation completed. Output saved to data/output")
```

The current API reflects the original coeval-population model. As the framework develops, formation history will become an explicit part of the population model rather than being implicit in the cluster simulation.

---

## Project Architecture

```text
realta/

├── .github/workflows/    # Continuous integration
├── src/
│   └── realta/
│       ├── binaries/     # Binary population generation and evolution
│       ├── data/         # Tabulated stellar models and remnant data
│       ├── imf/          # Initial Mass Function algorithms
│       ├── io/           # File loading, table parsing, and output
│       ├── simulation/   # Population simulation and orchestration
│       ├── stellar/     # Stellar lifetimes, ionisation, and remnants
│       ├── xray/         # X-ray luminosity and population calculations
│       ├── cli.py        # Command-line interface
│       ├── config.py     # Configuration dataclasses and YAML loading
│       └── random.py     # Reproducible NumPy random streams
├── tests/                # Unit and integration/regression tests
└── pyproject.toml        # Package configuration and dependencies
```

The architecture deliberately separates **population generation**, **stellar and binary evolution**, and **observable models**. This allows the current Power et al. implementation to provide a stable baseline while alternative physical prescriptions can be developed independently.

---

## Development Roadmap

The initial priority is to establish the Power et al. (2009) implementation as a well-tested and reproducible baseline.

Future development may include:

1. General stellar-formation histories beyond coeval populations.
2. Improved empirical and theoretical models for binary initial conditions.
3. More sophisticated stellar and binary evolution.
4. Optional interfaces to detailed stellar-evolution tools such as MESA.
5. Explicit compact-object population outputs.
6. Optional connections to Galactic-dynamics tools such as `galpy`.
7. Alternative physical prescriptions for applications including HMXBs and compact-object populations.

These extensions should be developed incrementally. The baseline model should remain reproducible as the code evolves.

---

## Development & Testing

### 1. Running Tests

Execute the unit and regression test suite with pytest:

```bash
python -m pytest
```

### 2. Code Formatting & Linting

`realta` uses Ruff for linting and formatting:

```bash
# Check code quality
ruff check .

# Check formatting
ruff format --check .
```

---

## Scientific Provenance

The current population and HMXB model is based on:

> Power, C., et al. (2009), MNRAS, 395, 1146.

The long-term development of Realta will maintain explicit provenance for physical prescriptions and distinguish the original reference model from alternative or newly developed prescriptions.

---

## License

Distributed under the MIT License. See `LICENSE` for more information.

