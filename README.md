# fastmode

[![DOI](https://img.shields.io/badge/DOI-10.17605%2FOSF.IO%2FX4T8M-blue)](https://doi.org/10.17605/OSF.IO/X4T8M)

Find fast GROMACS integration settings for a system, and check each one against
a 2 fs reference of that same system.

> **Correction notice (v0.2.0, 2026-09-27).** The DOI above points to v0.1.0, which
> contains errors. Its main protein claim, that `constraints = all-bonds` alone
> lets villin reach dt 6-7 fs, is false: without hydrogen mass repartitioning
> (HMR), every all-bonds run at dt 5 fs or more crashes. The setting it
> recommended (dt 6 fs, HMR factor 3, all-bonds) does pass, but the protein check
> it reported tested a different HMR factor, several table numbers were wrong,
> the shell speedups were divided by a reference from a different batch, and a
> units error meant the warmup was never excluded. This version corrects them,
> re-validates the villin settings with three seeds each, and changes the default
> HMR factor from 4.0 to 3.0. `CHANGELOG.md` lists every correction.

A production run usually uses `dt = 2 fs`. Larger timesteps can give a 2 to 3x
speedup, but a bad choice can fail without an error message. This tool runs a
2 fs reference, sweeps the settings, validates every candidate against the
reference on four checks, and reports the fastest setting that passes all four,
with the number behind every failure.

For a protein, two settings work together. With the usual `constraints =
h-bonds`, `grompp` stops villin at dt 4 fs, because heavy-atom bonds become
too fast for a larger step. `constraints = all-bonds` removes that limit, but
all-bonds alone crashes at dt 5 fs and above. Adding HMR (factor 3) carries it
to dt 6 fs, which passes every check in three seeds with small, measurable
shifts in density and water structure. See "Where the limits are" and
"Re-validation" below.

## Run it

```sh
# from the repository root; GMX defaults to `gmx` on PATH
python3 fastmode.py audit \
    --gro systems/npt_3.0.gro \
    --top systems/topol_3.0.top \
    --out out/water_small \
    --ref-ps 200
```

Outputs, in `--out`:

- `report.txt`   the human-readable report
- `results.json` every run and every metric, machine-readable
- `runs/`        one directory per candidate, with the mdp, tpr, log, edr, xtc

Options:

- `--gmx`          path to the `gmx` binary (default: `$GMX`, else `gmx` on PATH)
- `--ntomp`        OpenMP threads (default 4)
- `--ref-ps`       production length in ps for each run (default 200)
- `--strategy`     `staged` (default) or `grid`

Set `GMXLIB` only if your force fields live outside the GROMACS installation.

`hmr-check` checks the mass transform in isolation:

```sh
python3 fastmode.py hmr-check --top systems/villin/topol.top
```

`selfcheck` measures the verifier itself. It runs three or more independent 2 fs
samples of the same system, each with a different velocity seed, and reports the
largest difference between any two. The command reports the floor; it never
moves a threshold, and the audit does not use it automatically. Compare a
candidate's difference with the floor yourself before you trust a close verdict.

```sh
python3 fastmode.py selfcheck \
    --gro systems/npt_3.0.gro \
    --top systems/topol_3.0.top \
    --out out/selfcheck \
    --ref-ps 200
```

## Systems included

- `systems/npt_3.0.gro`, `npt_4.5.gro`, `npt_6.0.gro` and their topologies:
  TIP3P water boxes of 2652, 8916 and 21087 atoms (OPLS-AA water model files).
- `systems/villin/`: villin HP-35 in TIP3P, 14528 atoms (582 protein atoms,
  4648 waters, 2 chloride ions), `amber99sb` force field, with
  the equilibrated structure `eq.gro`, the cleaned `prot.pdb` and the mdp files
  used to build it.
- `vsites/build/`: the virtual-site rebuild of villin used by `stage2.py`.

The force fields are the ones distributed with GROMACS.

## What is swept

For each of these, the tool runs a full MD run and validates it.

- timestep `dt` in {2, 3, 4, 5} fs
- HMR on or off (solute only; see below)
- multiple time stepping (MTS): off, PME every 2nd step, PME every 3rd step
- neighbour-list interval `nstlist` in {10, 25, 50, 100}
- Verlet buffer tolerance in {0.005, 0.05}

The `audit` sweep keeps `constraints = h-bonds`. All-bond constraints, larger
timesteps and virtual sites are covered by `search.py` and `stage2.py`.

`--strategy staged` sweeps the timestep ladder for each topology (HMR off, and
HMR on for a solute), keeps the faster topology, then runs a coordinate descent
over MTS, `nstlist` and the tolerance. `--strategy grid` runs the whole cross
product, which is up to several hundred multi-minute runs.

## The four checks

Every candidate is compared with a 2 fs reference of the same system, run with
the same mdp template. A candidate fails if any check fails. The thresholds are
fixed.

1. Conserved energy drift below 0.02 kJ/mol/ps/atom, read from the log line
   `Conserved energy drift`. This covers the whole run, including warmup; in the
   re-validation runs the production-only drift agreed with it to within 0.0005.
2. Mean temperature within 1 K of the reference mean.
3. O-O radial distribution function: maximum absolute deviation from the
   reference below 0.02, from `gmx rdf -ref "name OW" -sel "name OW"`.
4. Mean density within 0.5 % of the reference mean.

Temperature, density and the RDF discard the first 2000 steps, and `mdrun` gets
`-resetstep 2000` so the timing excludes them too. Up to v0.2.0 a units error
(`fastmode.warmup_time_ps`) skipped only the first frame, so the committed audit,
search and selfcheck results include the first 4 ps of each 200 ps run. The
re-validation results use the corrected code.

What the checks do not cover: the RDF and density are water properties, so on a
protein system they test the solvent, not the protein. Use `valprotein.py` or
`experiments/revalidate.py` for protein observables. The thermostat coupling
(`tau-t = 0.1 ps`) is strong, which damps temperature errors. The systems use no
dispersion correction, which lowers the absolute TIP3P density; reference and
candidate share the choice.

## Hydrogen mass repartition

HMR is a mass transform on the topology, done here in Python. For each solute
hydrogen, the mass is multiplied by the HMR factor (default 3.0, so
1.008 -> 3.024 u; Hopkins et al. 2015) and the added mass is subtracted from its
bonded heavy atom, read from `[ bonds ]`. The transform writes the file, then
re-parses it and asserts that total mass is conserved to 1e-6 relative. Pure
water is skipped, since SETTLE already makes it rigid.

Version 0.1.0 used a default factor of 4.0 and labelled it standard. On villin
it leaves the leucine methyl carbons light enough that `grompp` rejects even
dt 4 fs. In v0.1.0 result files, a label with `hmron` and no `_f` suffix means
factor 4.0; from v0.2.0 every label names its factor.

## Virtual sites (stage 2)

`stage2.py` builds a virtual-site system and sweeps it.

```sh
python3 stage2.py build                 # pdb2gmx, solvate, ionise, em, 300 ps eq
python3 stage2.py sweep --ref-ps 200    # dt {4, 5, 6, 7} fs vs its own 2 fs ref
```

The rebuild uses `systems/villin/prot.pdb`, the cleaned HP-35 domain, with
`pdb2gmx -vsite h -heavyh -ff amber99sb -water tip3p -ignh`: 293 hydrogens
become virtual sites, 313 sites in total, mass 4083.784, charge +2. Then
`editconf -d 1.2 -bt cubic`, `solvate`, `genion -neutral`, minimisation, and a
300 ps equilibration.

The rebuilt system is a different system, so its numbers compare only with its
own 2 fs reference. Both that reference and every candidate use
`constraints = all-bonds`, so this sweep measures the timestep effect at fixed
constraints and does not measure the effect of the constraints themselves.

## Timing policy

A run is timed only when the load from other work is at or below 2.0, after
subtracting the sweep's own OpenMP threads from the 1-minute load average. If
any contributing run was untimed, the report withholds the speedup.

This rule does not make timings comparable across batches. The same 2 fs villin
reference measured 69.6 ns/day in the audit batch, 46.8 in the search batch and
64.8 in the shell batch. Compare a speedup only with the reference timed in the
same batch. `experiments/timing.py` times settings in interleaved rounds for
that reason; its own runs exceeded the load limit, as the re-validation section
records.

## Verification noise floor

`selfcheck` measures the smallest difference the checks can resolve. Three
independent 2 fs samples, 200 ps each:

| system      | atoms | temperature | density | RDF max dev |
|-------------|-------|-------------|---------|-------------|
| water small | 2652  | 0.931 K     | 0.172 % | 0.0059      |
| villin      | 14528 | 0.722 K     | 0.103 % | 0.0064      |
| check limit |       | 1.0 K       | 0.5 %   | 0.02        |

Each floor is the largest difference among three samples, which is a rough
estimate. The temperature floor on small water (0.93 K) nearly equals the 1 K
limit, so a correct setting can fail the temperature check by chance, and a
close temperature failure is not proof of a bad setting. These floors were
measured with the warmup-exclusion error described above. The re-validation
below adds a second, independent villin floor measured with the corrected code.

## Results: the audit

Idle 10-core ARM machine, 4 OpenMP threads, 200 ps per run, speedups against the
2 fs reference of the same batch. Full tables are in `results/out_*_report.txt`.

| system      | atoms | fastest passing   | speedup | drift  | RDF dev | density |
|-------------|-------|-------------------|---------|--------|---------|---------|
| water small | 2652  | dt 5 fs, no HMR   | 2.34x   | 0.0031 | 0.0124  | 0.04 %  |
| water med.  | 8916  | dt 5 fs, no HMR   | 2.40x   | 0.0011 | 0.0122  | 0.05 %  |
| water large | 21087 | dt 5 fs, no HMR   | 2.36x   | 0.0007 | 0.0156  | 0.20 %  |
| villin      | 14528 | dt 4 fs, no HMR   | 1.94x   | 0.0011 | 0.0117  | 0.00 %  |

The water boxes land on dt 5 fs with no HMR, no MTS, `nstlist 10`, tolerance
0.005. Villin stops at dt 4 fs under h-bond constraints (see "Where the limits
are").

Failures worth recording, each with its number:

- MTS with PME every 3rd step at dt 5 fs (a 15 fs slow-force interval) is
  unstable on every water box: drift 0.51 to 0.53 kJ/mol/ps/atom, temperature
  +4.4 to +5.3 K, RDF 0.028 to 0.030. A separate 20 ps probe on small water,
  checking drift and temperature only (`ml/mtsres.txt`), passed every
  slow-force interval up to 12 fs and failed every interval of 15 fs or more.
  On villin, PME every 3rd step at dt 4 fs (12 fs) crashed, and on large water
  PME every 2nd step at dt 5 fs with `nstlist 25` failed the RDF check (0.0201).
- MTS is never faster here: the best `mts 2` water run is 661 ns/day against 721
  for plain dt 5.
- `nstlist 100` on the small box is refused by `grompp`: the pair-list cut-off
  (1.516 nm) exceeds half the shortest box vector (1.4991 nm).
- Noise-sized failures appear at the thresholds: water small dt 3 fails density
  at 0.55 % (floor 0.172 %), and water small dt 5 with MTS PME every 2nd step,
  `nstlist 50` and tolerance 0.05 fails temperature at 298.90 K (floor 0.931 K).
  These are single samples.

`run_timed.sh` reproduces the audit.

## Where the limits are

**With h-bond constraints, villin stops at dt 4 fs.** `grompp` warns when a
bond's oscillation period is shorter than 5 timesteps, and the audit treats the
warning as a rejection. `experiments/grompp_limits.py` reproduces the bond it
names from the shipped villin files (`results/grompp_limits_2026-09.txt`):

| dt   | HMR factor        | `grompp` | bond named        | period |
|------|-------------------|----------|-------------------|--------|
| 4 fs | off, 2.5, 3, 3.5  | accepts  |                   |        |
| 4 fs | 4.0               | rejects  | Leu42 CG-CD1      | 18 fs  |
| 4 fs | 4.5               | rejects  | Leu42 CG-CD1      | 14 fs  |
| 5 fs | off               | rejects  | Asp44 CG-OD1      | 22 fs  |
| 5 fs | 2.5, 3            | rejects  | Trp64 CG-CD1      | 22 fs  |
| 5 fs | 3.5               | rejects  | Lys48 CE-NZ       | 21 fs  |
| 5 fs | 4.0, 4.5          | rejects  | Leu42 CG-CD1      | 18, 14 fs |

HMR moves mass from each hydrogen-bearing heavy atom to its hydrogens, which
lightens that atom and speeds up its bonds to other heavy atoms. The aspartate
CG-OD1 bond carries no hydrogens, so HMR leaves its 22 fs period unchanged, and
the fastest bond can never get slower than that. HMR therefore cannot raise the
h-bonds limit above dt 4 fs, and at factor 4 or more it lowers it. These are
`grompp`'s rule-of-thumb rejections; the runs were not attempted.

**All-bond constraints remove the bond limit, and need HMR.** From the villin
search (`results/search_villin_search.json`):

| constraints | HMR               | dt 4   | dt 5   | dt 6   | dt 7   | dt 8+  |
|-------------|-------------------|--------|--------|--------|--------|--------|
| h-bonds     | off, or 2.5-3.5   | pass   | grompp | grompp | grompp | grompp |
| all-bonds   | off               | screen only | crash  | crash  | crash  | crash  |
| all-bonds   | 2.5 to 4.5        | pass   | pass   | pass   | pass*  | crash  |

"screen only": the run passed the 20 ps stability screen but was not selected
for certification. \* dt 7 fs with factor 3 fails the RDF check on
certification (0.0206); the other factors pass. We have not identified which
motion sets the dt 8 fs limit.

**LINCS accuracy trades drift against structure.** At dt 6 fs with all-bonds,
the default LINCS settings (order 4, one iteration) give a drift of -0.0061 in
all three re-validation seeds, a systematic bias that plain dt 4 fs does not
show (+0.0011). `lincs-order = 8` with `lincs-iter = 2` removes it (+0.0010),
but moves the water structure and density further from the reference, far
enough to fail the RDF or density check in 6 of 9 comparisons (see below). We
have not established why. Default LINCS is the setting that passes.

**`all-angles` constraints fail on every setting.** Each run stops with too many
LINCS warnings (1458 at dt 4 fs). Angle constraints couple many constraints
together, and LINCS does not converge on them.

## Re-validation of the villin settings (`experiments/revalidate.py`)

Three seeds each of a 2 fs reference, plain dt 4 fs (h-bonds), and dt 6 fs with
HMR factor 3 and all-bonds, with default and with tight LINCS; 200 ps each, the
first 2000 steps excluded. Each candidate seed is compared with each reference
seed; the first column compares the references with each other, which is the
noise floor. Values are min / mean / max over pairs.
`results/revalidate_2026-09.txt` has the full output.

| quantity         | 2 fs vs 2 fs          | dt 4 fs, h-bonds      | dt 6 fs, default LINCS | dt 6 fs, tight LINCS  |
|------------------|-----------------------|-----------------------|------------------------|-----------------------|
| temperature (K)  | 0.18 / 0.44 / 0.67    | 0.03 / 0.28 / 0.63    | 0.03 / 0.33 / 0.63     | 0.03 / 0.33 / 0.75    |
| density (%)      | 0.03 / 0.17 / 0.26    | 0.04 / 0.16 / 0.31    | 0.11 / 0.30 / 0.41     | 0.17 / 0.38 / 0.55    |
| O-O RDF max dev  | 0.002 / 0.004 / 0.006 | 0.008 / 0.011 / 0.016 | 0.010 / 0.012 / 0.014  | 0.017 / 0.021 / 0.025 |
| Rg (nm)          | 0.004 / 0.006 / 0.008 | 0.000 / 0.004 / 0.008 | 0.001 / 0.008 / 0.021  | 0.000 / 0.004 / 0.009 |
| RMSF correlation | 0.83 / 0.86 / 0.90    | 0.67 / 0.83 / 0.95    | 0.81 / 0.89 / 0.98     | 0.72 / 0.86 / 0.97    |

| setting                          | four-check verdicts | drift   | mean density shift |
|----------------------------------|---------------------|---------|--------------------|
| dt 4 fs, h-bonds                 | 9 of 9 pass         | +0.0011 | +0.11 %            |
| dt 6 fs, HMR 3, default LINCS    | 9 of 9 pass         | -0.0061 | +0.30 %            |
| dt 6 fs, HMR 3, tight LINCS      | 3 of 9 pass         | +0.0010 | +0.38 %            |
| dt 7 fs, HMR 3, default LINCS (1 seed) | 2 of 3 pass   | -0.0068 |                    |
| dt 7 fs, HMR 3, tight LINCS (1 seed)   | 0 of 3 pass   | +0.0030 |                    |

What this shows, and what it does not:

- dt 4 fs (h-bonds) and dt 6 fs (HMR 3, all-bonds, default LINCS) pass every
  check against every reference.
- Both shift the water structure measurably: their RDF deviations sit above the
  reference-to-reference floor and below the 0.02 limit. dt 6 fs also raises
  the density by about 0.3 % in all three seeds, above the floor, and one of its
  seeds has an Rg 1.7 % above the reference mean. Its drift carries the LINCS
  bias described above. These are small, systematic effects inside the limits.
- The protein RMSF correlations of all candidates overlap the
  reference-to-reference range, but plain dt 4 fs reaches as low as 0.67, below
  the lowest reference pair (0.83). Three runs per setting cannot resolve RMSF
  differences of this size.
- The nine comparisons per setting share three candidate runs and three
  references, so they are not nine independent tests.

Version 0.1.0 reported one protein comparison for dt 6 fs all-bonds (RMSF
correlation 0.90). That run used HMR factor 4.0, not the factor 3 it
recommended, and had no noise floor to judge it against.

**Timing** (`experiments/timing.py`, `results/timing_2026-09.txt`): five
interleaved rounds, speedup = ns/day over the same round's 2 fs reference,
median and range.

| setting                                        | speedup                  |
|------------------------------------------------|--------------------------|
| dt 4 fs, h-bonds                               | 1.93x (range 1.54-1.94x) |
| dt 6 fs, HMR 3, all-bonds, default LINCS       | 2.67x (range 2.13-2.70x) |
| dt 6 fs, HMR 3, all-bonds, tight LINCS         | 2.56x (range 2.09-2.60x) |
| dt 7 fs, HMR 3, all-bonds, tight LINCS         | 2.97x (range 2.37-3.01x) |

These runs broke the tool's own timing rule: the machine carried other work,
with an external load of up to about 3.4 against the 2.0 limit. Interleaving
exposes every setting to the same background, so the ratios are usable, but
the absolute ns/day are low and the ranges are wide. Treat the speedups as
approximate.

**Recommendation.** For villin on a CPU, dt 4 fs with h-bonds (about 1.9x) is the
setting with the smallest systematic shifts. dt 6 fs with HMR factor 3 and
`constraints = all-bonds` at default LINCS settings (about 2.7x) passes every
check but shifts density, water structure and drift by small, measurable
amounts; use it if those shifts do not matter for your observables. Validate
either on your own system.

## Brute-force search (`search.py`)

A two-stage funnel searches the integration-setting space: screen many settings
at 20 ps for stability (drift, temperature, density), then certify the top ones
at 200 ps against all four checks. `search_dashboard.py` shows it live on port
8780. State is in `<out>/search.json` and resumes after a stop.

Villin search, 126 screen settings, 24 certified, 23 pass. Speedups are against
this batch's 2 fs reference (46.8 ns/day); densities are relative to that
reference.

| setting | ns/day | speedup | drift | RDF dev | density |
|---|---|---|---|---|---|
| dt 7 fs, HMR 4.5, all-bonds | 141.1 | 3.01x | -0.0147 | 0.0173 | 0.24 % |
| dt 6 fs, HMR 3, all-bonds   | 127.1 | 2.71x | -0.0061 | 0.0100 | 0.34 % |
| dt 7 fs, HMR 4.0, all-bonds | 137.8 | 2.94x | -0.0106 | 0.0152 | 0.40 % |

The dt 7 fs settings are close to the drift and RDF limits and are not
recommended.

## Virtual-site search (`search_vsites/`)

The same funnel on the virtual-site rebuild (`vsites/build/eq.gro`), to test
whether virtual sites and all-bonds stack. Reference 59.3 ns/day (200 ps
certify).

| setting | ns/day | speedup | drift | RDF dev | pass |
|---|---|---|---|---|---|
| dt 7, all-bonds | 190.1 | 3.21x | 0.0001 | 0.0299 | FAIL (RDF) |
| dt 6, all-bonds | 162.1 | 2.74x | -0.0011 | 0.0167 | PASS |
| dt 5, all-bonds | 138.9 | 2.34x | -0.0005 | 0.0109 | PASS |
| dt 4, all-bonds | 113.7 | 1.92x | 0.0001 | 0.0069 | PASS |
| dt 4, h-bonds | 114.1 | 1.93x | 0.0008 | 0.0094 | PASS |

They do not stack: the best passing setting is dt 6 fs all-bonds at 2.74x, about
the same as HMR plus all-bonds on the plain system. Above dt 4 the h-bonds
topology is refused by `grompp`, and dt 8 and up crash. The density of the
virtual-site runs rises steadily with the timestep (0.12 %, 0.21 %, 0.35 % at
dt 4, 5, 6 fs), a systematic trend below the 0.5 % limit.

## Spherical solvent shell and reaction field (`spherical.py`)

This is a model change, not an integration setting. It keeps the protein plus a
mobile water shell, freezes the outer water, and puts the droplet in a larger
box, which removes most of the bulk water. `spherical.py build` makes the
system; `spherical.py run` runs it. `--coulomb reaction-field` replaces PME.

Replicated validation (`validate10x.py`, 0.3 ns replicas). The periodic 2 fs
reference ran in the same batch: 5 replicas, 55.8 to 69.5 ns/day, mean 64.8
(`results/shell_timing_2026-09.txt`). The range in the speedup column comes from
that reference spread; the shell runs themselves varied by under 2 percent.

| configuration | n | Rg (nm) | RMSF corr | O-O order q | speedup |
|---|---|---|---|---|---|
| periodic, dt 2 fs, h-bonds | 5 | 0.9614 | 1.000 | 0.557 | 1.0x |
| shell, PME, dt 4 fs | 3 | 0.9627 | 0.969 | 0.403 | 4.6x (4.3-5.3x) |
| shell, reaction field 1.2 nm, dt 4 fs | 5 | 0.9579 | 0.988 | 0.406 | 8.8x (8.2-10.3x) |
| shell, reaction field 1.5 nm, dt 4 fs | 5 | 0.9583 | 0.968 | 0.405 | 6.9x (6.4-8.0x) |

The shell runs use `constraints = all-bonds` without HMR; the reference uses
h-bonds. So the speedup combines the model change, the larger timestep and the
constraint change.

The shell reproduces the protein observables, but the mobile water is
interfacial (tetrahedral order about 0.40 against a bulk 0.557), and in the
frozen variant the protein approaches the 0 K wall. It is valid for protein
observables, not for solvent thermodynamics.

The shell runs stop at dt 4 fs; dt 5 and 6 crash with a SETTLE error. Because
they use all-bonds without HMR, and the periodic system with all-bonds and no
HMR also crashes at dt 5 fs, the missing HMR is a candidate cause that was not
tested. Version 0.1.0 attributed the limit to the confined surface water.

`validate10x.py` reads the shell system from `spherical_v3/`, which is not
included: 3506 atoms (the protein, 974 waters, 2 chloride ions) in a 7 nm box.
Version 0.1.0 records it as `spherical.py build` with a 0.9 nm mobile shell and
a 0.1 nm frozen layer.

Version 0.1.0 reported single-run shell speedups of up to 12.9x, dividing by the
46.8 ns/day reference from a different batch. Against the same-batch reference
they are about 1.4 times lower.

## Soft boundary (`softshell.py`)

The frozen outer water is a rigid 0 K wall. `softshell.py build` replaces it with
a restrained water moleculetype (`BWA`, harmonic position restraints, in a local
`boundary.itp`), so the boundary water is thermostatted and can fluctuate: 974
mobile, 627 restrained, 5387 atoms, box 8 nm. Run with
`softshell.py run --coulomb reaction-field`.

A single 200 ps run at dt 4 fs gave 325.8 ns/day, RMSF correlation 0.922, and
mobile-water tetrahedral order 0.469 (std 0.261), closer to bulk than the frozen
shell (0.403, std 0.314). It does not raise the timestep ceiling: dt 5 still
fails with a SETTLE error at restraint stiffness k 1000 and k 200, with both
reaction field and PME. It was not timed against a same-batch reference.

## Machine-learning probes (`ml/`)

A set of GPU-free probes for a learned speedup. Most closed with a measured
failure; `ml/README.md` has the details. The Lyapunov probe does not measure
what it set out to: it ran with a stochastic thermostat and a fresh random seed
per run, so its 0.12 ps figure mixes thermostat noise with chaotic divergence.
Whether a learned propagator can be statistically accurate was not tested.

## Environment

The results were produced with a development build of GROMACS (2027.0-dev,
ARM NEON SIMD, no GPU) on an Apple M4, Python 3 with NumPy. Absolute ns/day will
differ on other machines and GROMACS versions.
