# Fast-mode MD audit: findings and recommendations

> **Corrected in v0.2.0 (2026-09-27), registered at DOI 10.17605/OSF.IO/AQ9ZV.**
> Version 0.1.0, the version registered at
> DOI 10.17605/OSF.IO/X4T8M, stated that `constraints = all-bonds` alone lets the
> villin protein reach dt 6-7 fs. That is false: without hydrogen mass
> repartitioning (HMR), every all-bonds run at dt 5 fs or more crashes. The
> setting it recommended (dt 6 fs, HMR factor 3, all-bonds) does pass, but for a
> different reason than it gave, and its validation was not sound. Several other
> numbers and explanations were also wrong. `CHANGELOG.md` lists every
> correction.

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
- `experiments/revalidate.py`, `experiments/timing.py`,
  `experiments/grompp_limits.py` - the September 2026 re-validation.
- `validate10x.py`, `valprotein.py`, `dashboard.py`, `run_timed.sh`.
- `ml/` - a documented set of machine-learning probes, most of them negative.

## Speed results (ARM Mac, 4 OpenMP threads)

Each speedup is against a 2 fs reference timed in the same batch. The same 2 fs
villin reference has measured between 45.1 and 69.6 ns/day across batches, so a
speedup is only meaningful against its own batch's reference.

| system | setting | speedup | what changes |
|---|---|---|---|
| water small/medium/large | dt 5 fs, no HMR | 2.34x / 2.40x / 2.36x | timestep |
| villin HP-35 | dt 4 fs, h-bonds, no HMR | 1.9x | timestep |
| villin HP-35 | dt 6 fs, HMR factor 3, all-bonds, default LINCS | 2.7x | timestep, masses, constraints |
| villin, spherical shell, reaction field 1.2 nm | dt 4 fs, all-bonds | 8.8x (8.2-10.3x) | model change |
| villin, spherical shell, PME | dt 4 fs, all-bonds | 4.6x (4.3-5.3x) | model change |

The villin speedups come from interleaved timing on a machine carrying other
work (external load up to about 3.4, above the tool's 2.0 limit), so they are
approximate. The shell ranges come from the spread of the five reference
replicas.

## The villin settings, re-validated

Three seeds each, against three independent 2 fs references, 200 ps with the
warmup excluded (`results/revalidate_2026-09.txt`):

- **dt 4 fs, h-bonds:** passes the four checks in 9 of 9 comparisons. Water
  structure shifts slightly (RDF deviation 0.008-0.016, above the 0.002-0.006
  reference-to-reference floor, below the 0.02 limit).
- **dt 6 fs, HMR 3, all-bonds, default LINCS:** passes 9 of 9. It shows three
  small systematic shifts, each above the floor and inside the limit: density
  +0.3 % in all three seeds, RDF deviation 0.010-0.014, and a drift of -0.0061
  from LINCS error. One seed's radius of gyration is 1.7 % above the reference
  mean.
- **dt 6 fs, HMR 3, all-bonds, LINCS order 8 and 2 iterations:** removes the
  drift bias (+0.0010) but passes only 3 of 9, failing on RDF (up to 0.025) or
  density (up to 0.55 %). Why tighter constraints move the structure further is
  not established.
- **dt 7 fs, HMR 3, all-bonds (one seed each):** 2 of 3 with default LINCS, 0 of
  3 with tight LINCS.

Every candidate group has RMSF correlations below the lowest reference pair
(0.83): plain dt 4 fs reaches 0.67, dt 6 fs with default LINCS 0.81, and with
tight LINCS 0.72. Three runs per setting cannot resolve RMSF differences of this
size, and the nine comparisons per setting are not independent.

`constraints = all-bonds` freezes every bond, which changes the protein's
configurational distribution, not only its dynamics; HMR changes only masses.
That does not explain the measured shifts: they are in water density and
structure, water is rigid in every run, and the search shows no extra density
shift from all-bonds at dt 4 fs. The cause of the dt 6 fs shifts is not
established.

Version 0.1.0 reported an RMSF correlation of 0.90 for the dt 6 fs setting. That
run used HMR factor 4.0, not factor 3, and there was no noise floor to judge it.

## The solvent shell

The shell is a model change: a droplet with a frozen (or restrained) outer water
layer in a larger box. Against the periodic 2 fs reference timed in the same
batch, it runs 4.6x with PME and 8.8x with a 1.2 nm reaction field. It
reproduces the protein observables (Rg within 0.4 percent, RMSF correlation
0.97-0.99, 0.3 ns replicas), but the mobile water is interfacial (tetrahedral
order 0.40-0.47 against a bulk 0.557). It is valid for protein observables, not
for solvent thermodynamics. The shell runs use all-bonds without HMR against an
h-bonds reference, so their speedup combines three changes. Version 0.1.0
reported up to 12.9x by dividing by a slower reference from a different batch.

## Where the limits are

1. **With h-bond constraints, villin stops at dt 4 fs.** `grompp` rejects a bond
   whose period is under 5 timesteps. Without HMR it names the aspartate CG-OD1
   bond (22 fs) at dt 5. HMR lightens every hydrogen-bearing heavy atom, so their
   bonds to other heavy atoms speed up: `grompp` names tryptophan CG-CD1 (22 fs)
   at factors 2.5 and 3, lysine CE-NZ (21 fs) at 3.5, and leucine CG-CD1 at 4
   (18 fs) and 4.5 (14 fs), where it rejects even dt 4. The aspartate bond has
   no hydrogens, so HMR cannot raise the limit above it
   (`results/grompp_limits_2026-09.txt`).
2. **All-bond constraints need HMR to go past dt 4.** Without HMR every
   all-bonds run at dt 5 fs or more crashed. With HMR factors 2.5 to 4.5 they
   certified at dt 5, 6 and 7 fs, except factor 3 at dt 7 fs, which failed the
   RDF check (0.0206); every setting at dt 8 fs crashed. What sets the dt 8
   limit is not identified.
3. **LINCS accuracy trades drift against structure** at dt 6 fs, as above.
4. **Virtual sites do not stack with all-bonds.** A virtual-site rebuild with
   all-bonds reaches 2.74x at dt 6 fs, about the same as HMR plus all-bonds, for
   more setup work.
5. **MTS.** PME every third step at dt 5 fs (15 fs) is unstable on every water
   box (drift about 0.5 against 0.02). A 20 ps probe on small water, checking
   drift and temperature only, passed every slow-force interval up to 12 fs and
   failed every interval of 15 fs or more. On villin, PME every third step at
   dt 4 fs (12 fs) crashed. MTS was never faster than plain dt 5 fs on water.
6. **The shell runs stop at dt 4 fs.** They use all-bonds without HMR, and the
   periodic system in that state also crashes at dt 5 fs, so the missing HMR is
   an untested candidate cause.
7. **`all-angles` constraints fail on every setting.** Each run stops with too
   many LINCS warnings (1458 at dt 4 fs); coupled angle constraints do not
   converge in LINCS.

## Machine-learning probes

All probes were GPU-free and cheap.

| route | result |
|---|---|
| learned next-position propagator | open: the divergence probe ran with a stochastic thermostat, so it does not measure a Lyapunov time; statistical accuracy was not tested |
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

- The RDF and density are water properties; the protein is checked only by the
  separate protein-observable comparison.
- The drift check reads GROMACS's whole-run drift. With the warmup excluded, the
  re-validation production drift agrees with it to within 0.0005.
- Up to v0.2.0 a units error skipped only the first frame of the warmup. The
  warmup is 2000 steps, so its length grows with the timestep (4 ps at 2 fs, 8
  ps at 4 fs, 14 ps at 7 fs), and each run is that warmup plus the production
  time. The committed audit, search and selfcheck results therefore average over
  the warmup as well; for a 20 ps search screen at dt 7 fs it was 14 of 34 ps.
  `fastmode.warmup_time_ps` holds the fix, and the re-validation uses it.
- On small water the temperature noise floor (0.93 K) nearly equals the 1 K
  limit, so a correct setting can fail by chance.
- No long-range dispersion correction is used, so TIP3P densities sit at
  970-994 kg/m^3 across the systems. Reference and candidate share the choice.
- The thermostat coupling (`tau-t = 0.1 ps`) is strong, which damps temperature
  errors and weakens the temperature check.
- Runs used a development build of GROMACS (2027.0-dev).

## Scope and novelty

Every speed ingredient here is standard: HMR, virtual sites, all-bonds
constraints, reaction field, and spherical boundaries. Systematic tests of
physical validity in MD already exist (Merz and Shirts, PLOS ONE 2018,
`physical_validation`). The contribution is narrower: a search over integration
settings in which every candidate is checked against a same-system reference,
and a record of what that found, including this version's own corrections.

## Recommendations

For a small solvated protein on a CPU, dt 4 fs with h-bonds (about 1.9x) has the
smallest density and drift shifts. dt 6 fs with HMR factor 3 and
`constraints = all-bonds` at default LINCS settings (about 2.7x) passes every
check but shifts density (+0.3 %) and drift (-0.006) by small, measurable
amounts; both settings shift the water structure by similar amounts. Use dt 6 fs
where those shifts do not matter. Do not tighten LINCS on that
setting without re-validating. Use the shell only as an explicitly labelled
model change for protein-only questions. Validate any of these on your own
system.
