# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-27

This release corrects errors in 0.1.0, the version registered at DOI
10.17605/OSF.IO/X4T8M. Three rounds of independent AI review found them: one of
0.1.0 and two of drafts of this correction, the first of which itself contained
a wrong recommendation. Each finding was checked against the committed data or new
runs before correction.

### Corrected

- **Main protein claim.** 0.1.0 said `constraints = all-bonds` alone lets villin
  reach dt 6-7 fs and that "the key lever is all-bonds, not HMR". Every
  all-bonds run without HMR at dt 5 fs or more crashed; every passing all-bonds
  run at dt 5-7 fs used HMR.
- **The recommended setting.** dt 6 fs, HMR factor 3, all-bonds, default LINCS
  passes the four checks in 9 of 9 comparisons (three seeds against three
  references), with small systematic shifts above the noise floor: density
  +0.3 %, RDF deviation 0.010-0.014, drift -0.0061. The first draft of 0.2.0
  recommended tightening LINCS (order 8, 2 iterations) because that removes the
  drift bias; with the checks run, tight LINCS passes only 3 of 9. It is not
  recommended.
- **The h-bonds limit.** 0.1.0 said HMR caps villin at dt 3 fs (about 1.5x)
  because of the aspartate CG-OD1 bond. HMR factors 2.5-3.5 pass at dt 4 fs. At
  dt 5 fs, `grompp` names Asp CG-OD1 without HMR, Trp CG-CD1 at factors 2.5 and
  3, Lys CE-NZ at 3.5, and Leu CG-CD1 at 4 and 4.5, which also fails at dt 4 fs
  (`experiments/grompp_limits.py`). The first draft of 0.2.0 named Leu at
  factor 3, which was wrong.
- **Warmup exclusion.** `evaluate()` computed the warmup as
  `WARMUP_STEPS * dt / 1000` with `dt` already in ps, so it skipped 0.004 ps of
  a 2000-step warmup that lasts 4 ps at 2 fs and 14 ps at 7 fs. Temperature,
  density and the RDF included the warmup; for a 20 ps search screen at dt 7 fs
  it was 14 of 34 ps. Fixed
  (`fastmode.warmup_time_ps`, with a test). The committed audit, search and
  selfcheck results carry the error; the re-validation was recomputed with the
  fix.
- **Protein validation.** 0.1.0 reported RMSF correlation 0.90 for the dt 6 fs
  setting, from a run with HMR factor 4.0, not 3, and with no noise floor.
- **Shell speedups.** 0.1.0 reported 7-13x (up to 12.9x) by dividing by a
  46.8 ns/day reference from a different batch. Against the same-batch
  reference they are 4.6x (PME) and 8.8x (reaction field 1.2 nm), with ranges
  from the reference replicas (`results/shell_timing_2026-09.txt`). 0.1.0
  attributed the shell's dt 4 fs ceiling to confined water; the shell runs use
  all-bonds without HMR, an untested alternative cause.
- **MTS.** "Every slow-force interval up to 12 fs passes" holds only for a 20 ps
  small-water probe on drift and temperature; villin PME-3 at dt 4 fs (12 fs)
  crashed.
- **Lyapunov probe.** 0.1.0 gave 1e-70; its formula gives about 1e-696. The
  probe also ran with a stochastic thermostat and a fresh seed per run, so it
  does not measure a Lyapunov time, and path accuracy is not what MD provides.
  The route is open, not closed.
- **Table numbers.** The villin search density column (reported 0.11/0.04/0.09 %,
  measured 0.24/0.34/0.40 %) and the virtual-site dt 4 h-bonds row (reported
  116.4 ns/day, 1.96x, drift 0.0016; measured 114.1, 1.93x, 0.0008).
- **Other statements.** The virtual-site "7 fs wall" was attributed to a bond
  that all-bonds constrains; `all-angles` failures are LINCS non-convergence;
  one water failure was MTS PME-2 with `nstlist 50`, not plain dt 5; the papers
  described an "unresolved" verdict the code never produces and said `audit`
  sweeps constraints and virtual sites; all-bond constraints change the
  configurational distribution, which 0.1.0 called "no model change"; the paper
  listed four water boxes and gave three; the paper PDF printed the author as
  "true".

### Changed

- Default HMR factor 4.0 -> 3.0 (1.008 -> 3.024 u). 0.1.0 labelled 4.0 as
  standard.
- Setting labels always name the HMR factor (`_f3`). In 0.1.0 results, `hmron`
  with no suffix means factor 4.0.
- `gmx` defaults to `$GMX`, then `gmx` on PATH; `GMXLIB` is set only when the
  user sets it. 0.1.0 looked for `build/bin/gmx` above the repository.
- Scripts read inputs from `systems/` instead of sibling directories.
- `hmr-check` writes to the current directory, not the script's directory.
- `valprotein.py` takes `--begin-ps` to skip the warmup.
- `search_dashboard.py` no longer truncates HMR factors to integers.
- `experiments/revalidate.py` no longer prints a speedup: its runs came from
  more than one batch. `--summary-only` rewrites the report from the JSON.

### Added

- `experiments/revalidate.py`: three-seed re-validation of every candidate
  group, with a reference-to-reference noise floor for the four checks and for
  protein observables, production-only drift, and signed density shifts.
- `experiments/timing.py`: interleaved timing rounds, speedup per round.
- `experiments/grompp_limits.py`: which bond `grompp` names for each setting.
- Input systems: the villin structure, topology and build files, and the medium
  and large water boxes, in `systems/`; the virtual-site `eq.gro`.
- Notes on rerunning the `ml/` probes (`ml/README.md`).
- Unit tests for the parsers, the mass repartition transform, the mdp builder,
  the log readers, the four-check validator, the HMR default and the warmup.
- Continuous integration workflow that runs the tests on every push and pull
  request.
- `CONTRIBUTING.md`, this changelog, and issue templates.

### Removed

- Absolute paths, user name and host name of the author's machine from result
  files and GROMACS headers.
- The v0.1.0 DOI from `CITATION.cff`; a DOI for 0.2.0 will follow its
  registration.

## [0.1.0] - 2026-09-12

### Added

- `fastmode audit`: sweep the integration-settings space and validate every
  candidate against a 2 fs reference on four checks.
- `fastmode selfcheck`: measure the verifier's own noise floor with independent
  replicas.
- `fastmode hmr-check`: prove that the hydrogen mass repartition transform
  conserves total mass.
- Hydrogen mass repartition transform in Python.
- Report writer with a ranked result and the number behind every failure.
- Documentation in `README.md` and a methods note in `paper/`.
