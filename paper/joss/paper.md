---
title: 'fastmode-md: a validator for fast-mode integration settings in GROMACS simulations'
tags:
  - molecular dynamics
  - GROMACS
  - integration settings
  - validation
  - reproducibility
  - Python
authors:
  - name: Daniel Reda
    orcid: 0000-0001-6160-6215
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 12 September 2026
bibliography: paper.bib
---

# Summary

`fastmode-md` is a command-line tool that finds the fastest safe integration
settings for a molecular dynamics (MD) system and proves the result against a
reference run. It is written in Python 3, depends only on NumPy, and drives the
user's existing GROMACS installation [@abraham2015gromacs;
@vanderspoel2005gromacs].

The user gives the tool a coordinate file and a topology file. The tool runs a
2 fs reference, then sweeps the integration settings, and it validates every
candidate against that reference on four checks. It reports the fastest setting
that passes every check, and it reports the measured number behind every
failure. A check that does not run is a failure, never a pass. The tool never
lowers a limit to make a candidate succeed, and it records a `grompp` refusal
as a refusal instead of forcing it with `-maxwarn`.

The tool addresses a failure mode that is common in MD practice. The settings
that give the largest speedup are often the settings that break the physics, and
a broken run still writes a trajectory. `fastmode-md` makes the correctness
comparison explicit and automatic.

# Statement of need

MD users gain speed by raising the timestep, by hydrogen mass repartitioning
(HMR) [@feenstra1999hmr; @hopkins2015hmr], by multiple time stepping (MTS)
[@tuckerman1992mts], and by shortening neighbour-list updates. Each change can
break the physics silently. The failure is hard to see, because the difference
between a thermal fluctuation and an integration error is often smaller than the
run-to-run scatter. GROMACS refuses some settings at `grompp` time, for example
a timestep that exceeds a bond's oscillation limit, but it does not tell the
user whether the settings that did run are correct.

Existing guides and benchmarks describe which settings are fast on a given
machine. They rarely test whether the settings that a user actually selected
reproduce a reference. `fastmode-md` fills that gap. It turns the choice of
integration settings into a measurement on the user's own system, with a
recorded verdict and a recorded number.

# Design and functionality

The tool has three subcommands: `audit`, `selfcheck`, and `hmr-check`.

## audit

`fastmode audit` runs a 2 fs reference and then a sweep of candidates. The sweep
varies the timestep, HMR, MTS, the neighbour-list update interval, the Verlet
buffer tolerance, the constraint set, and the use of virtual sites. Every
candidate is compared with the reference on four checks:

1. conserved energy drift below 0.02 kJ/mol/ps per atom,
2. mean temperature within 1 K of the reference mean,
3. maximum absolute deviation of the O-O radial distribution function below
   0.02, and
4. mean density within 0.5 percent of the reference.

A candidate passes only when all four checks pass, and the report names the
fastest passing candidate.

## selfcheck

A threshold is meaningless without the measurement scatter. `fastmode
selfcheck` runs independent replicas of the reference with different velocity
seeds and compares them. It reports the noise floor of each check. A difference
comparable to the floor is reported as unresolved, not as a pass. The user can
then see whether a check can resolve the effect that it is meant to test.

## hmr-check

HMR is a mass transform on the topology. `fastmode hmr-check` performs the
transform in Python, writes the new topology, re-parses the file that was
written, and asserts that total solute mass is conserved. Water is skipped,
because SETTLE makes water rigid already. The transform debits the added
hydrogen mass from the bonded heavy atom, and it supports several hydrogens on
one heavy atom.

## Timing policy

A run is timed only when the load from other work on the machine is below a
fixed limit. The tool subtracts its own thread count from the one-minute load
average before it applies the limit. When the machine is busy, the run is
recorded as untimed and no speedup is reported, because a shared machine can
change a timing by more than the effect being measured.

# Relation to existing software

`fastmode-md` is not a replacement for GROMACS. It uses GROMACS for every force
evaluation and every integration step, and it depends on GROMACS for the
reference trajectory. The tool adds a verification layer on top of the engine.
It is also not a benchmarking suite: its purpose is a correctness verdict, and
performance is recorded only as supporting evidence for the verdict.

# Limitations

The four checks are necessary but not sufficient. They test energy drift, the
temperature, the solvent structure, and the density. They do not test every
observable, and a user with a different scientific question should add a check
for that question. The default limits and the default production length are set
for the systems in the accompanying study; the `selfcheck` mode exists so that a
user can measure the floors on a different system before trusting a verdict.

# Acknowledgements

The author thanks the GROMACS developers for the simulation engine. The methods
and the validation results are described in a separate methods note, deposited
at <https://osf.io/x4t8m/>.

# References
