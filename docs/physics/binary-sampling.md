# Binary sampling distributions

Generalises the fixed formation-time distributions in
`docs/provenance.md` Section 1 into independently configurable ones.
Every default reproduces the ported baseline exactly.

## Binary fraction

`config.binary_fraction` (default `1.0`). At the default, every
`m1 >= mcut` star gets a companion, matching the baseline exactly (no
RNG draw is performed at this setting, so the baseline's random-number
stream is untouched). At `binary_fraction < 1`, a per-star Bernoulli
draw decides whether a companion is assigned; stars without one keep
`m1` populated with `m2 = 0`, `period = 0`, `a = 0` placeholders,
rather than being removed from the arrays (this generalises
`binary_prescription="single"`'s own mechanism, see
[`interaction-prescriptions.md`](interaction-prescriptions.md)).

`has_companion` also gates the pre-SN merger channel's eligibility
(`binaries/population.py::generate_population`), since a placeholder
`period = 0` would otherwise trivially satisfy
`period < p_merge_max_period`.

Implementation: `binaries/population.py::generate_population`.
Tests: `tests/test_binary_sampling_distributions.py`.

## Mass-ratio distribution

`config.mass_ratio_distribution`:
- `"uniform"` (default) — companion mass flat between `mcomp` and
  `m1` (the baseline distribution).
- `"flat_q"` — `m2 = m1 * Uniform(0, 1)`, i.e. uniform in
  `q = m2/m1`.

Implementation: `binaries/population.py::generate_population`.

## Period distribution

`config.period_distribution`:
- `"log_uniform"` (default) — log-flat between `pmin`/`pmax` (the
  baseline distribution).
- `"log_normal"` — truncated log-normal via `scipy.stats.truncnorm`,
  strictly bounded to `[pmin, pmax]`. `mu`/`sigma` are derived
  generically from `pmin`/`pmax` (not literature-sourced) rather than
  fit to any particular population.

Implementation: `binaries/population.py::generate_population`.

## Continuous IMF slope

`config.imf_slope: float | None` overrides `SalpeterIMF`'s default
`alpha = 2.35` (`imf_type = 1` only; silently ignored for
Kroupa/Chabrier, which have no single power-law slope to override).
`None` reproduces the default exactly.

`SalpeterIMF.cdf()`'s denominator (`mmax**beta - mmin**beta`,
`beta = 1 - alpha`) is singular at `alpha = 1.0`; `config.py`'s
`__post_init__` rejects `imf_slope == 1.0` explicitly (a log-form CDF
for that case is not implemented).

Implementation: `imf/factory.py::get_imf`, `slope` parameter.
Tests: `tests/test_imf.py`.
