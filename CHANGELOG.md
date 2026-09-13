# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Unit tests for the parsers, the mass repartition transform, the mdp builder,
  the log readers and the four-check validator.
- Continuous integration workflow that runs the tests on every push and pull
  request.
- `CONTRIBUTING.md`, this changelog, and issue templates.

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
