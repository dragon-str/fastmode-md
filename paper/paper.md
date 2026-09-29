---
title: "Auditing fast-mode GROMACS settings against a same-system reference"
author: "Daniel Reda (ORCID 0000-0001-6160-6215)"
date: "2026-09-27 (version 0.2.0; corrects version 0.1.0 of 2026-09-12)"
---

# Correction notice

Version 0.1.0 of this note (DOI 10.17605/OSF.IO/X4T8M) contains errors. Its
main protein result, that `constraints = all-bonds` alone lets the villin
protein reach a 6-7 fs timestep, is false: without hydrogen mass repartitioning
(HMR), every all-bonds run at 5 fs or more crashed. The setting it recommended
(6 fs, HMR factor 3, all-bonds) does pass the checks, but its protein
validation tested a different HMR factor, its account of the h-bonds limit was
wrong, its shell speedups were divided by a reference from a different batch, a
units error meant the warmup was never excluded from the checks, and its
Lyapunov figure was both miscalculated and confounded by a stochastic
thermostat. This version corrects each error and re-validates the villin
settings with three seeds each. Three rounds of independent AI review, of the
released version and of two drafts of this correction, found these errors; the
first draft itself contained a wrong recommendation, which the second round
caught.

# Abstract

Integration settings can make a molecular dynamics simulation several times
faster, but the fastest setting that runs is not always correct. We present an
open-source tool that searches the integration-settings space and validates
every candidate against a 2 fs reference of the same system with four checks:
conserved energy drift, temperature, the oxygen-oxygen radial distribution
function, and density. We use it on pure water and on the villin HP-35 domain on
a commodity CPU. Water reaches dt 5 fs (about 2.4x). For villin, h-bond
constraints stop at dt 4 fs (about 1.9x), because `grompp` rejects heavy-atom
bonds that are too fast for a larger step; hydrogen mass repartitioning cannot
lift this limit, and at a mass factor of 4 it lowers it. All-bond constraints
remove the limit but crash without HMR. With HMR factor 3 they reach dt 6 fs
(about 2.7x), which passes all four checks in three seeds but shifts density,
water structure and energy drift by small amounts that exceed the run-to-run
noise. Tightening the LINCS constraint solver removes the drift bias but moves
the water structure further, and that setting fails the checks in 6 of 9
comparisons. A spherical solvent shell
reaches 4.6-8.8x but changes the model.

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
   line `Conserved energy drift` (whole run; in the re-validation the
   production-only drift agrees to within 0.0005).
2. Mean temperature within 1 K of the reference mean.
3. Maximum absolute deviation of the O-O radial distribution function below 0.02.
4. Mean density within 0.5 percent of the reference.

The first 2000 steps (4 ps at 2 fs, up to 14 ps at 7 fs) are meant to be
excluded from checks 2-4; each run is that warmup plus the production time. Up
to version 0.2.0 a units error excluded only the first frame. The audit, search
and selfcheck results reported below carry that error, and the re-validation
does not. A missing check
is a failure. The RDF and density describe the water; for the protein we add a
separate comparison of the radius of gyration, backbone RMSD and per-residue
RMSF.

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
three samples is a rough estimate of the floor: the re-validation below, with
the warmup correctly excluded, measures a villin density floor of up to 0.26 %,
about 2.5 times the value in this table.

## Timing

A run is timed only when the external machine load is below 2. This does not
make timings comparable between batches: the same 2 fs villin reference has run
at between 45.1 and 69.6 ns/day across batches. We report each speedup against a
reference timed in the same batch. For the villin settings we timed every
setting in five interleaved rounds and report the median ratio; those rounds
ran under an external load of up to about 3.4, so the speedups are approximate.

# Results

## Water

All three water boxes reach dt 5 fs with no HMR, no MTS, `nstlist 10` and a
buffer tolerance of 0.005, for 2.34x, 2.40x and 2.36x. MTS with PME every third
step at dt 5 fs (a 15 fs slow-force interval) is unstable on every box (drift
about 0.5 against the 0.02 limit). A 20 ps probe on small water, checking drift
and temperature only, passed every slow-force interval up to 12 fs and failed
every interval of 15 fs or more. MTS was never faster than plain dt 5 fs.

## Villin: where the limits are

**h-bond constraints stop at dt 4 fs.** `grompp` warns when a bond's
oscillation period is shorter than five timesteps. Without HMR it names the
aspartate CG-OD1 bond (22 fs) at dt 5 fs. HMR moves mass from each
hydrogen-bearing heavy atom to its hydrogens, which lightens that atom and
speeds up its bonds to other heavy atoms. With HMR, `grompp` names tryptophan
CG-CD1 (22 fs) at factors 2.5 and 3, lysine CE-NZ (21 fs) at 3.5, and leucine
CG-CD1 at 4 (18 fs) and 4.5 (14 fs), where it rejects even dt 4 fs. The
aspartate bond carries no hydrogens, so HMR cannot raise the limit above its
22 fs period. These are rule-of-thumb rejections; the runs were not attempted.

**All-bond constraints need HMR.** Constraining every bond removes the
heavy-atom bond limit. Without HMR, every all-bonds run at dt 5 fs or more
crashed. With HMR factors 2.5 to 4.5, all-bonds runs certified at dt 5, 6 and
7 fs, except factor 3 at dt 7 fs, which failed the RDF check (0.0206); every
setting at dt 8 fs crashed. We have not identified what sets the dt 8 fs limit.

**All-angles constraints fail.** Every all-angles run stopped with too many
LINCS warnings (1458 at dt 4 fs); coupled angle constraints do not converge in
LINCS.

**Virtual sites do not stack with all-bonds.** A virtual-site rebuild reaches
2.74x at dt 6 fs, about the same as HMR plus all-bonds on the plain system. Its
dt 7 fs candidate fails the RDF check (0.0299).

## Villin: re-validation

We ran three seeds each of the 2 fs reference, plain dt 4 fs (h-bonds), and
dt 6 fs with HMR factor 3 and all-bonds with default LINCS (order 4, one
iteration) and with tight LINCS (order 8, two iterations), 200 ps each with the
warmup excluded. We compared every candidate seed with every reference seed; the
reference-to-reference differences give the noise floor (min / mean / max):

| quantity | 2 fs vs 2 fs | dt 4 fs, h-bonds | dt 6 fs, default LINCS | dt 6 fs, tight LINCS |
|---|---|---|---|---|
| density (%) | 0.03 / 0.17 / 0.26 | 0.04 / 0.16 / 0.31 | 0.11 / 0.30 / 0.41 | 0.17 / 0.38 / 0.55 |
| RDF max deviation | 0.002 / 0.004 / 0.006 | 0.008 / 0.011 / 0.016 | 0.010 / 0.012 / 0.014 | 0.017 / 0.021 / 0.025 |
| Rg difference (nm) | 0.004 / 0.006 / 0.008 | 0.000 / 0.004 / 0.008 | 0.001 / 0.008 / 0.021 | 0.000 / 0.004 / 0.009 |
| RMSF correlation | 0.83 / 0.86 / 0.90 | 0.67 / 0.83 / 0.95 | 0.81 / 0.89 / 0.98 | 0.72 / 0.86 / 0.97 |
| energy drift | +0.0008 | +0.0011 | -0.0061 | +0.0010 |
| four checks passed | | 9 of 9 | 9 of 9 | 3 of 9 |

Plain dt 4 fs and dt 6 fs with default LINCS pass every comparison. Both shift
the water structure above the noise floor and below the limit. dt 6 fs also
raises the density by about 0.3 % in all three seeds and carries a systematic
drift of -0.0061 from LINCS error, and one of its seeds has an Rg 1.7 % above
the reference mean. Tight LINCS removes the drift bias but moves the water
structure further (RDF deviation 0.017-0.025 against 0.010-0.014), and the
setting fails the checks in 6 of 9 comparisons, with its mean RDF deviation
(0.021) at the limit; its larger density shift is within the noise. We have not
established why. Single runs at dt 7 fs pass 2 of 3 comparisons with default
LINCS and none with tight LINCS. Every candidate group has RMSF correlations
below the lowest reference pair (0.83): 0.67 for plain dt 4 fs, 0.81 for dt 6 fs
with default LINCS and 0.72 with tight LINCS. Three runs per setting cannot
resolve RMSF differences of this size, and the nine comparisons per setting
share runs, so they are not independent.

In five interleaved timing rounds, plain dt 4 fs runs 1.93x (range 1.54-1.94x)
and dt 6 fs with default LINCS 2.67x (range 2.13-2.70x), each relative to the
same round's 2 fs reference.

Version 0.1.0 reported one protein comparison for dt 6 fs all-bonds (RMSF
correlation 0.90). That run used HMR factor 4, not factor 3, and there was no
noise floor to interpret it.

# The solvent shell: a model change

Removing bulk water and keeping only a solvent shell gives the largest numbers.
A droplet of mobile water inside a box, with an outer water layer frozen or
restrained, runs at dt 4 fs with all-bond constraints and no HMR. Against a
periodic 2 fs h-bonds reference timed in the same batch (5 replicas of 0.3 ns,
55.8-69.5 ns/day, mean 64.8):

| configuration | n | Rg (nm) | RMSF corr. | O-O order q | speedup |
|---|---|---|---|---|---|
| periodic, dt 2 fs | 5 | 0.9614 | 1.000 | 0.557 | 1.0x |
| shell, PME, dt 4 fs | 3 | 0.9627 | 0.969 | 0.403 | 4.6x (4.3-5.3x) |
| shell, reaction field 1.2 nm, dt 4 fs | 5 | 0.9579 | 0.988 | 0.406 | 8.8x (8.2-10.3x) |
| shell, reaction field 1.5 nm, dt 4 fs | 5 | 0.9583 | 0.968 | 0.405 | 6.9x (6.4-8.0x) |

The ranges come from the spread of the reference replicas. The shell
reproduces the protein observables, but the mobile water is interfacial
(tetrahedral order about 0.40 against a bulk 0.557), and in the frozen variant
the protein approaches the boundary. A position-restrained boundary moves the
water toward bulk order (q 0.469 in a single run) but does not raise the dt 4 fs
ceiling. Because the shell runs use all-bonds without HMR, and the periodic
system in that state also crashes at dt 5 fs, the missing HMR is an untested
candidate cause of that ceiling. The shell is valid for protein observables, not
for bulk-solvent thermodynamics. Version 0.1.0 reported 7-13x for these
configurations by dividing by a slower reference from a different batch.

# Machine-learning probes

We tested GPU-free machine-learning routes to a general speedup:

| route | outcome |
|---|---|
| learned next-position propagator | open: the divergence probe ran with a stochastic thermostat and does not measure a Lyapunov time; statistical accuracy was not tested |
| position-only force model | closed: R^2^ 0.81; the residual is PME plus the constraint force |
| learned MTS splitting | closed: the limit is a resonance, not a smooth bias |
| local water-patch proposal | closed: the Hastings penalty exceeds the gain |
| Boltzmann generator, Cartesian | closed: unphysical bond geometry; proposal energy about 1e18 |
| Boltzmann generator, rigid body | closed: never learns excluded volume; acceptance 0/2000 |
| one-site coarse-grained water | closed: iterative Boltzmann inversion diverges; the first RDF peak is 25 percent too low |

Version 0.1.0 called the propagator route closed, citing a required accuracy of
1e-70 for a 200 ps path; the same formula gives about 1e-696. Both numbers rest
on a divergence rate measured with a stochastic thermostat and a fresh random
seed per run, which mixes thermostat noise into the growth, so the probe does
not measure a Lyapunov time. And path accuracy is not what MD provides, so even
a clean measurement would not close the route.

On this hardware a neural network is slower than optimized C per force
evaluation, so a learned method can only win by reducing work or by generating
samples, and the sample generators we tried failed on the sharp structure of a
dense liquid.

# Discussion

Every speed ingredient here is standard: HMR, all-bond constraints, virtual
sites, reaction field and spherical boundaries. The contribution is an audit
procedure and a record of what it found, including its own failures.

Those failures are the useful part. Version 0.1.0 passed its own checks and
still drew a wrong conclusion: the write-up credited all-bond constraints
without checking the all-bonds runs without HMR, which had crashed; the protein
check ran on a different setting from the one recommended; single runs were
compared without a noise floor; and a units error kept the warmup in three of
the four checks. The first draft of this correction repeated the pattern: it recommended
tighter LINCS settings on the strength of the drift alone, without running the
four checks on them, and those settings fail. Replicated runs, a noise floor for
every reported quantity, running every check on every recommended setting, and
same-batch timing are what caught these errors. We recommend all four for any
claim of this kind.

Constraining every bond is not a neutral choice: holonomic constraints change
the protein's configurational distribution, not only its dynamics. That does not
explain the shifts we measured, which are in water density and structure; water
is rigid in every run, and the search shows no extra density shift from all-bond
constraints at dt 4 fs. The cause of the dt 6 fs shifts is not established.

# Conclusion

For villin on a CPU, plain dt 4 fs with h-bonds (about 1.9x) passes the four
checks with the smallest density and drift shifts. dt 6 fs with HMR factor 3,
`constraints = all-bonds` and default LINCS settings (about 2.7x) also passes,
with small measurable shifts in density (+0.3 %) and drift (-0.006); both
settings shift the water structure by similar amounts. It suits questions those
shifts do not affect. With tightened LINCS settings it fails the checks in 6 of
9 comparisons. Both
are results for one small protein on one machine and need validation on other
systems. A spherical shell reaches 4.6-8.8x but must be labelled a model change.

# Data and code availability

The tool, the search harness, the validation and re-validation scripts, the
input systems and the reduced results are at
`https://github.com/dragon-str/fastmode-md` under the MIT license. Trajectories
are not stored. The core tool, the re-validation, the timing and the `grompp`
diagnosis run from the repository; the shell validation needs a shell system
that is not included, and the machine-learning probes expect the original
working-tree layout (see `README.md` and `ml/README.md`). Version 0.1.0 is
registered on OSF at DOI `10.17605/OSF.IO/X4T8M`.

# Use of AI

AI coding agents (Claude models from Anthropic, and DeepSeek V4.1 Flash served
by Fireworks AI) wrote most of the code, ran the simulations and drafted the
text, working from a specification that an AI assistant drafted and the author
approved. AI reviewers found the
errors that version 0.2.0 corrects. The author directed the work and is
responsible for its content.

# Acknowledgements

The author thanks the GROMACS developers for the simulation engine.
