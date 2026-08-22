# `realta`

[![CI](https://github.com/doctorcbpower/realta/actions/workflows/ci.yml/badge.svg)](https://github.com/doctorcbpower/realta/actions/workflows/ci.yml)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`realta` is a Monte Carlo simulation code designed for modeling High-Mass X-ray Binaries (HMXRBs) and stellar populations in globular clusters. In its initial commit, it's a refactoring of the FORTRAN90 code used in [Power et al. (2009, MNRAS, 395, 1146)](https://ui.adsabs.harvard.edu/abs/2009MNRAS.395.1146P/abstract).

---

## Key Features

- **Monte Carlo Population Synthesis:** Simulates binary evolution and cluster dynamics over customized timescales.
- **Custom IMF Models:** Supports standard Initial Mass Functions including Kroupa, Salpeter, and Chabrier models.
- **Physical Data Tables:** Tabulated stellar lifetimes, remnant mass distributions, and ionizing photon rates.
- **Clean Command Line Interface:** Execute cluster runs via simple CLI inputs or run programmatic simulations directly in Python.

---

## Installation

### Prerequisites

- Python `>= 3.8` (Tested on 3.10 through 3.13)
- `pip` package manager

### 1. Basic Installation

Clone the repository and install the package locally:

```bash
git clone [https://github.com/doctorcbpower/realta.git](https://github.com/doctorcbpower/realta.git)
cd realta
pip install .
```

### 2. Developer Installation

To install in editable mode with development dependencies (`pytest`, `ruff`, `mypy`):

```
pip install -e ".[dev]"
```

## Quickstart & Usage

### 1. Command Line Interface (CLI)

Once installed, you can trigger a simulation run directly using the realta executable:

```
# Run a cluster simulation with custom parameters
realta --ntot 1000 --tmax 10.0 --output-dir ./output
```

### 2. Python API 

You can also use realta as a Python module in your analysis scripts:

```
from realta import ClusterSimulation, SimulationConfig

# Initialize configuration
config = SimulationConfig()
config.ntot = 5000         # Total number of initial binary systems
config.tmax = 12.0         # Max simulation time (Gyr)
config.iseed = 42          # Random seed

# Run cluster simulation
cluster = ClusterSimulation(config)
results = cluster.run(output_dir="data/output")

print(f"Simulation completed. Output saved to data/output")
```

## Project Architecture
```
realta/
├── .github/workflows/    # CI/CD pipelines (Ruff & Pytest matrix)
├── src/
│   └── realta/
│       ├── binaries/     # Binary population generation & orbital dynamics
│       ├── data/         # Tabulated stellar models & remnant data
│       ├── imf/          # Initial Mass Function algorithms
│       ├── io/           # File loading, table parsing, and export routines
│       ├── simulation/   # Core cluster simulation loop
│       ├── stellar/      # Stellar lifetimes, ionisation, and remnants
│       ├── xray/         # Luminosity & X-ray population calculations
│       ├── cli.py        # Command Line Interface logic
│       ├── config.py     # Configuration dataclasses and YAML loaders
│       └── random.py     # Reproducible NumPy random stream generator
├── tests/                # Unit and integration regression suites
└── pyproject.toml        # Build configurations & dependencies
```

## Development & Testing

## 1. Running Tests
Execute the unit and regression testing suite with pytest:
```
python3 -m pytest
```
## 2. Code Formatting & Linting
`realta` uses Ruff for fast code linting and formatting:
```
# Check code quality
ruff check .
# Check formatting
ruff format --check .
```

## License
Distributed under the MIT License. See `LICENSE` for more information.
