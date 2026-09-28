# Fast-mode MD audit: findings and recommendations

> **Corrected in v0.2.0 (2026-09-27).** Version 0.1.0, the version registered at
> DOI 10.17605/OSF.IO/X4T8M, stated that `constraints = all-bonds` alone lets the
> villin protein reach dt 6-7 fs. That is false: without hydrogen mass
> repartitioning (HMR), every all-bonds run at dt 5 fs or more crashes. The
> working setting needs both. Several other numbers and explanations were also
> wrong. `CHANGELOG.md` lists every correction.

## Goal

Audit GROMACS integration settings and validate every candidate against a 2 fs
reference of the same system. The settings are timestep, HMR, multiple time
stepping (MTS), neighbour-list frequency, Verlet buffer tolerance, constraints,
and virtual sites. Four checks must pass: conserved energy drift, temperature,
O-O RDF, and density. Every failure is reported with its number, and no
threshold is lowered.

## What the tool provides

- `fastmode.py` - the `audit` tool, plus `hmr-check` and `selfcheck` (noise
  floor).
- `stage2.py` - virtual-site rebuild and sweep.
- `search.py` + `search_dashboard.py` - a two-stage search over the
  integration-settings space, with a live dashboard.
- `spherical.py`, `softshell.py` - solvent-boundary variants.
- `validate10x.py`, `valprotein.py`, `experiments/revalidate.py`,
  `dashboard.py`, `run_timed.sh`.
- `ml/` - a documented set of machine-learning probes, most of them negative.

## Speed results (ARM Mac, 4 OpenMP threads)

Each speedup is measured against a 2 fs reference timed in the same batch. The
same 2 fs villin reference ran at 46.8 ns/day in one batch and 69.6 in another,
so a speedup is only meaningful against its own batch's reference.

| system | setting | speedup | what changes |
|---|---|---|---|
| water small/medium/large | dt 5 fs, no HMR | 2.34x / 2.40x / 2.36x | timestep only |
| villin HP-35 | dt 4 fs, h-bonds, no HMR | 1.94x | timestep only |
| villin HP-35 | dt 6 fs, HMR factor 3, all-bonds | 2.56x with tight LINCS; 2.67x with default LINCS | HMR and all-bond constraints |
| villin HP-35 | dt 7 fs, HMR factor 4.5, all-bonds | 3.01x | marginal; not recommended |
| villin, spherical shell, reaction field 1.2 nm | dt 4 fs, all-bonds | 8.8x | model change |
| villin, spherical shell, PME | dt 4 fs, all-bonds | 4.6x | model change |

`constraints = all-bonds` freezes every bond vibration, including heavy-atom
bonds. That is a change to the model's dynamics, even though it is a standard
one. HMR changes the atomic masses. Neither changes the equilibrium ensemble in
principle, but both are setting changes that need validation, which is what the
four checks and the protein-observable comparison provide.

The dt 6 fs setting was re-validated in September 2026 with three seeds against
three independent 2 fs references (`experiments/revalidate.py`,
`results/revalidate_2026-09.txt`). It passes the four checks in 9 of 9
comparisons. Its backbone RMSF correlates with the references at 0.82 to 0.98,
while two 2 fs references correlate with each other at only 0.84 to 0.91, so the
difference is within run-to-run noise. Its radius of gyration differs from the
references by 0.001 to 0.020 nm, against 0.004 to 0.008 nm between references;
one seed sits about 2 percent high. Its O-O RDF deviation (0.010 to 0.013) is
above the reference noise floor (0.002 to 0.006) and below the 0.02 limit, so
the setting shifts water structure by a small, measurable amount. Plain dt 4 fs
behaves the same way (RDF deviation 0.007 to 0.016).

Version 0.1.0 reported an RMSF correlation of 0.90 for this setting. That run
used HMR factor 4.0, not the factor 3.0 setting it recommended, and no noise
floor existed to judge it against.

The shell results are a model change: a droplet in a periodic box with a frozen
(or restrained) outer water shell. Against a periodic 2 fs reference timed in the
same batch (64.8 ns/day, 5 replicas), the shell runs 4.6x with PME (3 replicas)
and 8.8x with a 1.2 nm reaction field (5 replicas). They reproduce the protein
ensemble (Rg within 0.4 percent, RMSF correlation 0.97-0.99, 0.3 ns replicas),
but the mobile water is interfacial, not bulk (tetrahedral order 0.40-0.47
against a bulk 0.557), and in the frozen variant the protein approaches the
0 K wall. They are validated for protein observables, not for solvent
thermodynamics. Version 0.1.0 reported 9.9-12.9x for the shell, by dividing
by a slower reference from a different batch.

## Where the limits are

1. **With h-bond constraints, the protein stops at dt 4 fs.** `grompp` rejects a
   bond whose oscillation period is shorter than 5 timesteps. Without HMR the
   fastest heavy-atom bond is the aspartate CG-OD1 (period 22 fs), which fails
   at dt 5. HMR cannot help: it moves mass from carbon to hydrogen, which makes
   the leucine methyl carbons lighter, so their CG-CD1 bond speeds up. At HMR
   factor 3 that bond has the same 22 fs period; at factor 4 it drops to 18 fs
   and `grompp` rejects even dt 4. These are `grompp`'s rule-of-thumb
   rejections, not observed instabilities.
2. **All-bonds constraints need HMR to go past dt 4.** Constraining every bond
   removes the heavy-atom bond limit. Without HMR, every all-bonds run at dt 5 fs
   or more crashes. With HMR (factors 2.5 to 4.5), all-bonds runs certify at
   dt 5, 6 and 7 fs, and every setting at dt 8 fs crashes. We have not
   identified which motion sets the dt 8 limit.
3. **Default LINCS settings bias the energy drift at large timesteps.**
   With the default LINCS settings (order 4, one
   iteration), the dt 6 fs all-bonds setting drifts at -0.0061 kJ/mol/ps per atom
   in all three seeds, so the drift is systematic, not noise. `lincs-order = 8`
   with `lincs-iter = 2` brings it to +0.0010, the same as plain dt 4 fs. At
   dt 7 fs the same change moves the drift from -0.0068 to +0.0030. Every value
   is inside the 0.02 limit, but the tighter settings remove the bias.
4. **Virtual sites do not stack with all-bonds.** A virtual-site rebuild with
   all-bonds reaches 2.74x at dt 6 fs, about the same as the HMR plus all-bonds
   route, for more setup work.
5. **MTS with PME every third step is unstable at dt 5 fs** (a 15 fs slow-force
   interval): drift about 0.5 against the 0.02 limit on every water box. Every
   MTS setting with a slow-force interval of 12 fs or less passes, and every
   setting at 15 fs or more fails.
6. **The shell timestep ceiling is 4 fs.** dt 5 fails in every boundary and
   electrostatics combination tested, which points at the confined surface
   water rather than the boundary.
7. **`all-angles` constraints fail on every setting.** Each run stops with too
   many LINCS warnings (1458 at dt 4 fs). Angle constraints couple many
   constraints together, and LINCS does not converge on them.

## Machine-learning probes

All probes were GPU-free and cheap.

| route | result |
|---|---|
| learned next-position propagator | open: the Lyapunov time (0.12 ps) rules out path accuracy, but MD itself is not path-accurate; statistical accuracy was not tested |
| position-only force model | closed: R^2 0.81; the residual is PME plus the constraint force |
| learned MTS correction | closed: the instability is a resonance, not a bias |
| local water-patch proposal | closed: the Hastings penalty exceeds the gain |
| Boltzmann generator (Cartesian) | closed: unphysical bond geometry, U ~ 1e18 |
| Boltzmann generator (rigid-body) | closed: never learns excluded volume, acceptance 0 |
| one-site CG water (pair + 3-body IBI) | closed: first peak 25 percent too low |
| one-site CG water (pair only) | closed: IBI diverges |

A neural network is slower than optimized C per force evaluation on this
hardware, so a learned method can only win by reducing work or by generating
samples. The sample-generating routes we tried failed on the sharp structure of
a dense liquid: stiff bonds, excluded volume, and the first RDF peak.

## Limits of the checks

- The four checks test energy drift, temperature, water structure (O-O RDF) and
  density. The RDF and density are water properties; the protein is checked
  only by the separate protein-observable comparison.
- The drift check reads GROMACS's whole-run drift, which includes the 2000
  warmup steps. Temperature, density and the RDF exclude them.
- On small water, the temperature noise floor (0.93 K) nearly equals the 1 K
  limit, so a correct setting can fail by chance.
- The systems use no long-range dispersion correction (`DispCorr`), which is why
  the TIP3P density sits near 970-987 kg/m^3. Reference and candidate share
  the choice, so the comparison is fair, but the absolute density is low.
- The thermostat coupling (`tau-t = 0.1 ps`) is strong, which damps temperature
  errors and weakens the temperature check.
- Runs used a development build of GROMACS (2027.0-dev).

## Scope and novelty

Every speed ingredient here is standard: HMR, virtual sites, all-bonds
constraints, reaction field, and spherical boundaries. Systematic tests of
physical validity in MD already exist (Merz and Shirts, PLOS ONE 2018,
`physical_validation`). The contribution is narrower: a search over
integration settings in which every candidate is checked against a same-system
reference, and a documented set of failures, including this version's own
corrections.

## Recommendations

For a small solvated protein on a CPU: dt 6 fs with HMR factor 3 and
`constraints = all-bonds`, with `lincs-order = 8` and `lincs-iter = 2`, which runs 2.56x faster than dt 2 fs (median of five interleaved rounds) and passes the four checks and the protein comparison within noise in three seeds. Keep dt 4 fs with
h-bonds as the conservative choice. Use the shell only as an explicitly labelled
model change for protein-only questions.
