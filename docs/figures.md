# Figures

Analysis scripts, one config in → one figure out. See each script's
own module docstring for CLI usage.

| Figure | Script | Config | Shows |
|---|---|---|---|
| 1 — Population evolution | `scripts/run_paper1_experiment.py` | `configs/paper1_basic_experiment.yml` | `L_bol`/`L_UV`/`Q_H`/`L_X` vs. time, across every implemented HMXB-activation mode: the five `binary_prescription` values ([`physics/interaction-prescriptions.md`](physics/interaction-prescriptions.md)) plus `post_sn_rlof`/`wind_capture` ([`physics/post-sn-rlof.md`](physics/post-sn-rlof.md), [`physics/wind-capture.md`](physics/wind-capture.md)) — the latter two resolved via `VARIANT_OVERRIDES` (not `binary_prescription` values; layered on the `non_interacting` base). |
| 2 — X-ray/UV evolution | `scripts/run_paper1_experiment.py` | `configs/paper1_basic_experiment.yml` | `L_X/L_UV(t)`, same variant set as Figure 1. `post_sn_rlof` gives the strongest signature of the seven curves; `wind_capture` a clear second; `non_interacting`/`standard_interaction`/`enhanced_interaction` overlap closely (see [`physics/rlof-classifier.md`](physics/rlof-classifier.md)'s donor-selection note). |
| 3 — Effect of mergers | `scripts/figure3_merger_effects.py` | `configs/figure3_merger_effects.yml` | Luminosity evolution and cumulative compact-object/merger counts across `non_interacting`/`standard_interaction`/`enhanced_mergers`. `standard_interaction` and `enhanced_mergers` overlap closely — the same donor-selection property as the Xu et al. cross-check ([`physics/rlof-classifier.md`](physics/rlof-classifier.md)). No compact-object-type (WD/NS/BH) census — `RemnantTable` returns mass only. |
| 4 — IMF vs. binary-fraction degeneracy | `scripts/figure4_imf_binary_grid.py` | `configs/figure4_imf_binary_grid.yml` | `L_X/L_UV` on an `(alpha_IMF, f_bin)` grid at selected ages, using the continuous Salpeter slope and `binary_fraction` ([`physics/binary-sampling.md`](physics/binary-sampling.md)). No interaction prescription varied. Small-`ntot` cells (steep IMF slope) are small-number-statistics noisy. |
| 5 — Metallicity sweep | `scripts/figure5_metallicity_sweep.py` | `configs/figure5_metallicity_sweep.yml` | `L_X/L_UV(t)` for `non_interacting` vs. `standard_interaction`, at `imetal = 1, 2, 3`. At `imetal = 1` (Z=0) the two curves are identical — RLOF classification is skipped entirely at Z=0 ([`physics/rlof-classifier.md`](physics/rlof-classifier.md)). |
| 6 — Stochastic realisations | `scripts/figure6_stochastic_realisations.py` | `configs/figure6_stochastic_realisations.yml` | `P(L_X/L_UV \| M_cl, Z, t)` — distribution, not mean, across repeated realisations (fixed config, varying `iseed`) at several cluster-mass proxies (`ntot`; actual `M_cl` read from `BinaryPopulation.total_mass_msun`). Not Paper 2's parameter-sweep machinery (`Event`/`PopulationHistory` does not exist yet, see `docs/known-gaps.md`) — a standalone repeat-seed runner. A realisation with zero active HMXBs contributes to a reported "quiescent fraction" rather than a finite log-ratio point. |

Every script's tests are smoke tests (config parsing, a small real
run, figure-file existence), not numeric regression pins — these are
comparison/illustration figures built on already-tested underlying
quantities, not new physics paths of their own.
