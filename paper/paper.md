---
title: "Auditing fast-mode molecular dynamics settings without fooling yourself"
author: "Daniel Reda (ORCID 0000-0001-6160-6215)"
date: "2026-09-27 (version 0.2.0; corrects version 0.1.0 of 2026-09-12)"
---

# Correction notice

Version 0.1.0 of this note (DOI 10.17605/OSF.IO/X4T8M) contains errors. Its
main protein result, that `constraints = all-bonds` alone lets the villin
protein reach a 6-7 fs timestep, is false: without hydrogen mass repartitioning
(HMR), every all-bonds run at 5 fs or more crashed. Its protein-observable
check tested a different HMR factor from the one it recommended, its account of
the h-bonds limit was incomplete, its shell speedups were divided by a
reference from a different batch, and its Lyapunov bound was wrong by more than
600 orders of magnitude. This version corrects each error, adds replicated
re-validation and interleaved timing, and states where the audit method itself
fell short. The corrections came from an independent review of the released
version.

# Abstract

Integration settings can make a molecular dynamics simulation several times
faster, but the fastest setting that runs is not always correct. We present an
open-source tool that searches the integration-settings space and validates
every candidate against a 2 fs reference of the same system with four checks:
conserved energy drift, temperature, the oxygen-oxygen radial distribution
function, and density. We use it on pure water and on the villin HP-35 domain on
a commodity CPU. Water reaches dt 5 fs (about 2.4x). For villin, h-bond
constraints stop at dt 4 fs, because `grompp` rejects heavy-atom bonds that are
too fast for a larger step; hydrogen mass repartitioning cannot lift this limit,
and at a mass factor of 4 it lowers it. All-bond constraints remove the limit
but crash without HMR; with HMR factor 3 they reach dt 6 fs, which passes all
four checks and a protein-observable comparison within run-to-run noise in
three seeds, at 2.56x (range 2.09-2.60x) with tightened LINCS settings. Default LINCS
settings add a systematic energy drift at that timestep. A spherical solvent
shell reaches 4.6-8.8x but changes the model. We also report machine-learning
routes we tested, most of which failed on measured grounds.

# Statement of need

Molecular dynamics users raise the timestep, repartition hydrogen masses, switch
on multiple time stepping, and shorten neighbour-list updates to gain speed.
Each change can break the physics without an error message, and the difference
between a thermal fluctuation and an integration error is often smaller than the
run-to-run scatter. GROMACS refuses some settings at `grompp` time, but it does
not tell the user whether the settings that did run are correct.

Systematic tests of physical validity exist: the `physical_validation` package
(Merz and Shirts, PLOS ONE 2018) tests integrator convergence and ensemble
validity. The tool described here is narrower. It searches over integration
settings and compares every candidate with a same-system reference on four
fixed checks, recording `grompp` refusals as refusals rather than forcing them
with `-maxwarn`.

# Method

## Systems

We use three TIP3P water boxes of 2652, 8916 and 21087 atoms and the villin
HP-35 domain in TIP3P, 14528 atoms (582 protein atoms, 4648 waters, 2 chloride
ions, `amber99sb`). We also build a virtual-site variant of the protein with
`pdb2gmx -vsite h -heavyh`. Runs use a development build of GROMACS
(2027.0-dev) on an ARM CPU with four OpenMP threads and no GPU. No dispersion
correction is applied, which lowers the absolute TIP3P density; reference and
candidates share the choice. The thermostat is v-rescale with a strong
coupling (`tau-t = 0.1 ps`).

## The four checks

A candidate must satisfy all of the following against the 2 fs reference of the
same system:

1. Conserved energy drift below 0.02 kJ/mol/ps per atom, from the GROMACS log
   line `Conserved energy drift` (whole run, including 2000 warmup steps; the
   production-only drift agrees to within 0.0005 in our re-validation).
2. Mean temperature within 1 K of the reference mean.
3. Maximum absolute deviation of the O-O radial distribution function below 0.02.
4. Mean density within 0.5 percent of the reference.

A missing check is a failure. The RDF and density describe the water; for the
protein we add a separate comparison of the radius of gyration, backbone RMSD
and per-residue RMSF.

## Noise floor

The tool's `selfcheck` mode runs three independent 2 fs replicas with different
velocities and reports the largest pairwise difference for each check. At 200 ps:

| system | temperature (K) | density (percent) | O-O RDF |
|---|---|---|---|
| water, 2652 atoms | 0.931 | 0.172 | 0.0059 |
| villin, 14528 atoms | 0.722 | 0.103 | 0.0064 |

The tool reports these floors but does not use them in the verdict; a verdict is
a fixed-threshold comparison. On small water the temperature floor nearly
equals the 1 K limit, so a correct setting can fail by chance. A maximum over
three samples is a rough estimate of the floor.

## Timing

A run is timed only when the external machine load is below 2. This does not
make timings comparable between batches: the same 2 fs villin reference ran at
46.8, 64.8 and 69.6 ns/day in three batches. We therefore report each speedup
against a reference timed in the same batch, and for the recommended setting we
time all settings in five interleaved rounds and report the median ratio.

# Results

## Water

All three water boxes reach dt 5 fs with no HMR, no MTS, `nstlist 10` and a
buffer tolerance of 0.005, for 2.34x, 2.40x and 2.36x. MTS with PME every third
step at dt 5 fs (a 15 fs slow-force interval) is unstable on every box (drift
about 0.5 against the 0.02 limit); every slow-force interval of 12 fs or less
passes. MTS was never faster than plain dt 5 fs.

## Villin: where the limits are

**h-bond constraints stop at dt 4 fs.** `grompp` warns when a bond's
oscillation period is shorter than five timesteps. Without HMR it names the
aspartate CG-OD1 bond (22 fs) at dt 5 fs. With HMR it names the leucine CG-CD1
bond instead: HMR moves mass from each carbon to its hydrogens, which lightens
the methyl carbons and shortens that bond's period to 22 fs at factor 3,
18 fs at factor 4 and 14 fs at factor 4.5. So HMR cannot raise the h-bonds
limit, and at factor 4 or more `grompp` rejects even dt 4 fs. These are
rule-of-thumb rejections; the runs were not attempted. Plain dt 4 fs passes all
four checks at 1.94x.

**All-bond constraints need HMR.** Constraining every bond removes the
heavy-atom bond limit. Without HMR, every all-bonds run at dt 5 fs or more
crashed. With HMR factors 2.5 to 4.5, all-bonds runs certified at dt 5, 6 and
7 fs, except factor 3 at dt 7 fs, which failed the RDF check (0.0206); every
setting at dt 8 fs crashed. We have not identified what sets the dt 8 fs limit.

**Default LINCS settings bias the drift.** At dt 6 fs with all-bonds, the
default LINCS settings (order 4, one iteration) gave a drift of -0.0061 in all
three seeds; `lincs-order = 8` with `lincs-iter = 2` gave +0.0010, the same as
plain dt 4 fs. At dt 7 fs the drift moved from -0.0068 to +0.0030.

**All-angles constraints fail.** Every all-angles run stopped with too many
LINCS warnings (1458 at dt 4 fs); coupled angle constraints do not converge in
LINCS.

**Virtual sites do not stack with all-bonds.** A virtual-site rebuild reaches
2.74x at dt 6 fs, about the same as HMR plus all-bonds on the plain system. Its
dt 7 fs candidate fails the RDF check (0.0299).

## Villin: re-validation of the recommended setting

We ran three seeds each of the 2 fs reference, plain dt 4 fs, and dt 6 fs with
HMR factor 3 and all-bonds, 200 ps each, and compared every candidate seed with
every reference seed. The reference-to-reference differences give a noise floor
for each quantity, including the protein observables.

| quantity (min / mean / max) | 2 fs vs 2 fs | dt 4 fs, h-bonds | dt 6 fs, HMR 3, all-bonds |
|---|---|---|---|
| RDF max deviation | 0.002 / 0.004 / 0.006 | 0.007 / 0.011 / 0.016 | 0.010 / 0.012 / 0.013 |
| Rg difference (nm) | 0.004 / 0.005 / 0.008 | 0.000 / 0.003 / 0.008 | 0.001 / 0.008 / 0.020 |
| RMSF correlation | 0.84 / 0.87 / 0.91 | 0.68 / 0.83 / 0.95 | 0.82 / 0.89 / 0.98 |

Both candidates pass the four checks in 9 of 9 comparisons, and their protein
observables fall inside the reference-to-reference spread, except one dt 6 fs
seed with an Rg about 2 percent high. Their RDF deviations exceed the noise floor
while staying under the limit: the larger timesteps shift water structure by a
small, measurable amount.

In five interleaved timing rounds, plain dt 4 fs runs 1.93x (range 1.54-1.94x) and dt 6 fs
runs 2.67x (range 2.13-2.70x) with default LINCS and 2.56x (range 2.09-2.60x) with the tightened settings
(median, with range, relative to the same round's 2 fs reference).

Version 0.1.0 reported one protein comparison for this setting (RMSF
correlation 0.90). That run used HMR factor 4, not factor 3, and there was no
noise floor to interpret it; the reference-to-reference correlation alone
ranges from 0.84 to 0.91.

# The solvent shell: a model change

Removing bulk water and keeping only a solvent shell gives the largest numbers.
A droplet of mobile water inside a box, with an outer water layer frozen or
restrained, runs at dt 4 fs. Against a periodic 2 fs reference timed in the same
batch (64.8 ns/day, 5 replicas of 0.3 ns):

| configuration | n | Rg (nm) | RMSF corr. | O-O order q | speedup |
|---|---|---|---|---|---|
| periodic, dt 2 fs | 5 | 0.9614 | 1.000 | 0.557 | 1.0x |
| shell, PME, dt 4 fs | 3 | 0.9627 | 0.969 | 0.403 | 4.6x |
| shell, reaction field 1.2 nm, dt 4 fs | 5 | 0.9579 | 0.988 | 0.406 | 8.8x |
| shell, reaction field 1.5 nm, dt 4 fs | 5 | 0.9583 | 0.968 | 0.405 | 6.9x |

The shell reproduces the protein observables, but the mobile water is
interfacial (tetrahedral order about 0.40 against a bulk 0.557), and in the
frozen variant the protein approaches the boundary. A position-restrained
boundary moves the water toward bulk order (q 0.469 in a single run) but does
not raise the dt 4 fs ceiling. The shell is valid for protein observables, not
for bulk-solvent thermodynamics. Version 0.1.0 reported 7-13x for these
configurations by dividing by a slower reference from a different batch.

# Machine-learning probes

We tested GPU-free machine-learning routes to a general speedup:

| route | outcome |
|---|---|
| learned next-position propagator | open: the Lyapunov time is 0.12 ps, so no propagator stays on one exact path beyond about 1 ps; MD itself does not either, and statistical accuracy was not tested |
| position-only force model | closed: R^2^ 0.81; the residual is PME plus the constraint force |
| learned MTS splitting | closed: the limit is a resonance, not a smooth bias |
| local water-patch proposal | closed: the Hastings penalty exceeds the gain |
| Boltzmann generator, Cartesian | closed: unphysical bond geometry; proposal energy about 1e18 |
| Boltzmann generator, rigid body | closed: never learns excluded volume; acceptance 0/2000 |
| one-site coarse-grained water | closed: iterative Boltzmann inversion diverges; the first RDF peak is 25 percent too low |

Version 0.1.0 called the propagator route closed, citing a required accuracy of
1e-70 for a 200 ps path. The correct bound from the same formula is about
1e-696, and the argument only rules out path accuracy, which the MD integrator
does not provide either.

On this hardware a neural network is slower than optimized C per force
evaluation, so a learned method can only win by reducing work or by generating
samples, and the sample generators we tried failed on the sharp structure of a
dense liquid.

# Discussion

Every speed ingredient here is standard: HMR, all-bond constraints, virtual
sites, reaction field and spherical boundaries. The contribution is an audit
procedure and a record of what it found, including its own failures.

Those failures are instructive. Version 0.1.0 passed its own checks and still
drew a wrong conclusion, because the write-up credited all-bond constraints
without checking the all-bonds runs without HMR, which had crashed; the
protein check ran on a
different setting from the one recommended, and single runs were compared
without a noise floor. Replicated runs, a reference-to-reference floor for every
reported quantity, and same-batch timing corrected the result. We recommend
them for any claim of this kind.

# Conclusion

For villin on a CPU, dt 6 fs with HMR factor 3, `constraints = all-bonds`,
`lincs-order = 8` and `lincs-iter = 2` passes the four checks and a
protein-observable comparison within noise in three seeds, at 2.56x (range 2.09-2.60x).
Plain dt 4 fs with h-bonds is the conservative choice at 1.93x (range 1.54-1.94x). Both are
results for one small protein on one machine and need validation on other
systems. A spherical shell reaches 4.6-8.8x but must be labelled a model change.

# Data and code availability

The tool, the search harness, the validation scripts, the input systems and the
reduced results are at `https://github.com/dragon-str/fastmode-md` under the MIT
license. Trajectories are not stored. The core tool, the re-validation and the
timing run from the repository; the machine-learning probes expect the original
working-tree layout (see `ml/README.md`). Version 0.1.0 is registered on OSF at
DOI `10.17605/OSF.IO/X4T8M`.

# Use of AI

AI coding agents (Claude models from Anthropic, and DeepSeek V4.1 Flash served by Fireworks AI) wrote most of the code, ran the simulations and drafted the text, working from the author's written specification. An AI reviewer found the errors that version 0.2.0 corrects. The author directed the work and is responsible for its content.

# Acknowledgements

The author thanks the GROMACS developers for the simulation engine.
