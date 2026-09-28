# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-27

This release corrects errors in 0.1.0, the version registered at DOI
10.17605/OSF.IO/X4T8M. An independent review of 0.1.0 found them; each was
checked against the committed data before correction.

### Corrected

- **Main protein claim.** 0.1.0 said `constraints = all-bonds` alone lets villin
  reach dt 6-7 fs and that "the key lever is all-bonds, not HMR". In the search
  data every all-bonds run without HMR at dt 5 fs or more crashed; every
  passing all-bonds run at dt 5-7 fs used HMR. The recommended setting needs
  both.
- **The h-bonds limit.** 0.1.0 said HMR caps villin at dt 3 fs (about 1.5x)
  because of the aspartate CG-OD1 bond. HMR factors 2.5-3.5 pass at dt 4 fs.
  `grompp` names CG-OD1 only without HMR; with HMR it names the leucine CG-CD1
  bond, whose period HMR shortens (22 fs at factor 3, 18 fs at factor 4). The
  dt 4 fs rejection with HMR on came from the default factor of 4.0.
- **Protein validation.** 0.1.0 reported RMSF correlation 0.90 for the
  recommended dt 6 fs setting, but that run used HMR factor 4.0, not factor 3,
  and no noise floor existed. The setting is now re-validated with three seeds
  against three independent references (`results/revalidate_2026-09.*`).
- **Shell speedups.** 0.1.0 reported 7-13x (up to 9.9-12.9x) by dividing by a
  46.8 ns/day reference from a different batch. Against the 64.8 ns/day
  reference timed with the shells, they are 4.6x (PME) to 8.8x (reaction field
  1.2 nm). The 9.9x figure came from a single run never committed.
- **Table numbers.** The villin search density column (reported 0.11/0.04/0.09 %,
  measured 0.24/0.34/0.40 %) and the virtual-site dt 4 h-bonds row (reported
  116.4 ns/day, 1.96x, drift 0.0016; measured 114.1, 1.93x, 0.0008).
- **Lyapunov bound.** 0.1.0 gave 1e-70 for the accuracy a learned propagator
  needs over 200 ps; its own formula gives about 1e-696. The argument also only
  rules out path accuracy, which MD does not have either, so the route is open,
  not closed.
- **Other statements.** The virtual-site "7 fs wall" was attributed to a bond
  that all-bonds constrains; `all-angles` failures are LINCS non-convergence,
  not a step-0 failure; MTS PME-3 is unstable at dt 5 fs, not at every
  timestep; drift includes the warmup; the papers described an "unresolved"
  verdict that the code never produces, and said `audit` sweeps constraints and
  virtual sites, which it does not; the paper listed four water boxes and gave
  three; the paper PDF printed the author as "true".

### Changed

- Default HMR factor 4.0 -> 3.0 (1.008 -> 3.024 u), the common choice. 0.1.0
  labelled 4.0 as standard.
- Setting labels always name the HMR factor (`_f3`). In 0.1.0 results, `hmron`
  with no suffix means factor 4.0.
- `gmx` defaults to `$GMX`, then `gmx` on PATH; `GMXLIB` is set only when the
  user sets it. 0.1.0 looked for `build/bin/gmx` above the repository.
- Scripts read inputs from `systems/` instead of sibling directories.
- `valprotein.py` takes `--begin-ps` to skip the warmup.

### Added

- `experiments/revalidate.py`: three-seed re-validation with a
  reference-to-reference noise floor for the four checks and for protein
  observables, production-only drift, and a LINCS accuracy test.
- `experiments/timing.py`: interleaved timing rounds, speedup per round.
- Input systems: the villin structure, topology and build files, and the medium
  and large water boxes, in `systems/`; the virtual-site `eq.gro`.
- Notes on rerunning the `ml/` probes (`ml/README.md`).
- Unit tests for the parsers, the mass repartition transform, the mdp builder,
  the log readers and the four-check validator.
- Continuous integration workflow that runs the tests on every push and pull
  request.
- `CONTRIBUTING.md`, this changelog, and issue templates.

### Removed

- Absolute paths of the author's machine from result files and headers.

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
