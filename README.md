# fastmode

[![DOI](https://img.shields.io/badge/DOI-10.17605%2FOSF.IO%2FX4T8M-blue)](https://doi.org/10.17605/OSF.IO/X4T8M)

Find the fastest safe GROMACS setting for a system, and prove it on that system.

A production run usually uses `dt = 2 fs` and no virtual sites. A correct setup
reaches 4 to 6 fs, which is about a 2 to 3x speedup. The settings are fiddly, and
a bad choice fails silently, so this tool measures the claim instead of asserting
it. It runs a 2 fs reference, sweeps the settings, validates every candidate
against the reference on four checks, and reports the fastest setting that
passes every check. It reports every failure with its number.

For a protein, the lever is `constraints = all-bonds`, not hydrogen virtual
sites. On this build, `grompp` refuses `dt >= 4 fs` with plain hydrogen mass
repartition (HMR) alone, because the carboxylate bond CG-OD1 has an oscillational
period of 2.2e-02 ps and no hydrogen transform can lengthen it. That caps HMR at
3 fs, or about 1.5x. Adding `constraints = all-bonds` constrains the remaining
heavy-atom bonds and lets the plain protein reach `dt` 6 to 7 fs with no warning.
Hydrogen virtual sites are an alternative that removes every hydrogen from the
dynamics, but a rebuilt virtual-site system and a plain all-bonds system reach the
same speed, so virtual sites are not required. Stage 2 builds the virtual-site
system and sweeps it; see below.

## Run it

```sh
# from the repository root
export GMXLIB=/path/to/gromacs/share/top

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

- `--gmx`          path to the `gmx` binary (default `build/bin/gmx`)
- `--ntomp`        OpenMP threads (default 4; the fastmode track uses 4)
- `--ref-ps`       production length in ps for each run (default 200)
- `--strategy`     `staged` (default) or `grid`

`hmr-check` proves the mass transform in isolation:

```sh
python3 fastmode.py hmr-check --top ../phase2b/topol.top
```

`selfcheck` measures the verifier itself. It runs three or more independent 2 fs
samples of the same system, each with a different velocity seed, and reports the
largest difference between any two. That difference is the smallest change the
four checks can resolve. The command reports the floor; it never moves a
threshold. Run it before you trust any RDF verdict.

```sh
python3 fastmode.py selfcheck \
    --gro systems/npt_3.0.gro \
    --top systems/topol_3.0.top \
    --out out/selfcheck \
    --ref-ps 200
```

## What is swept

For each of these, the tool runs a full MD run and validates it.

- timestep `dt` in {2, 3, 4, 5} fs
- hydrogen mass repartition (HMR) on or off (solute only; see below)
- multiple time stepping (MTS): off, PME every 2nd step, PME every 3rd step
- neighbour-list interval `nstlist` in {10, 25, 50, 100}
- Verlet buffer tolerance in {0.005, 0.05}

`--strategy staged` runs the search in two parts. First it sweeps the timestep
ladder for each topology (HMR off, and HMR on for a solute). Hydrogen
repartition does not speed a run up at a fixed timestep; it enables a larger
one, so the two must be compared over the whole ladder. It keeps the faster
topology. Then it is a coordinate descent over MTS, `nstlist` and the
tolerance, swept from that base so the gains compose. `--strategy grid` runs
the whole cross product. The grid is up to several hundred multi-minute runs,
so the staged search is the default.

## The four checks

Every candidate is compared with a 2 fs reference of the same system, run with
the same mdp template. A candidate fails if any check fails. The thresholds are
fixed and are never lowered.

1. Conserved energy drift below 0.02 kJ/mol/ps/atom, read from the log line
   `Conserved energy drift`.
2. Mean temperature within 1 K of the reference mean.
3. O-O radial distribution function: maximum absolute deviation from the
   reference below 0.02, from `gmx rdf -ref "name OW" -sel "name OW"`.
4. Mean density within 0.5 % of the reference mean.

The first 2000 steps are discarded for every measurement, and `mdrun` gets
`-resetstep 2000`.

## Hydrogen mass repartition

HMR is a mass transform on the topology, done here in Python. For each solute
hydrogen, the mass is multiplied by 4.0 and the added mass is subtracted from
its bonded heavy atom, read from `[ bonds ]`. The transform asserts that total
mass is conserved to 1e-6 relative. Pure water is skipped: SETTLE already makes
it rigid, so repartition buys nothing. For a pure-water system the tool records
HMR as not applicable.

An HMR copy of the topology is written into the run directory, and any local
`#include` files beside the original topology are copied with it. The transform
writes the file first, then re-parses the file it wrote and checks conservation
on that artifact, not on the arithmetic that produced it.

## Virtual sites (stage 2)

`stage2.py` builds a virtual-site system and sweeps it. It is the path to more
than 3 fs on a protein.

```sh
python3 stage2.py build                 # pdb2gmx, solvate, ionise, em, 300 ps eq
python3 stage2.py sweep --ref-ps 200    # dt {4, 5, 6, 7} fs vs its own 2 fs ref
```

The rebuild uses `prot.pdb`, the cleaned HP-35 domain. The related `1yrf.pdb`
carries an extra N-terminal segment (residues 35 to 41), which makes `pdb2gmx`
split the chain and stop. Build it with
`pdb2gmx -vsite h -heavyh -ff amber99sb -water tip3p -ignh`: 293 hydrogens
become virtual sites, 313 sites in total, mass 4083.784, charge +2. Then
`editconf -d 1.2 -bt cubic`, `solvate`, `genion -neutral`, minimisation, and a
300 ps equilibration.

The rebuilt system is a different system. Its numbers compare only with its own
2 fs reference, never with the `phase2b` numbers. Both the reference and every
candidate use `constraints = all-bonds`; the reference uses the same constraints
so the comparison is fair.

The sweep covers `dt` in {4, 5, 6, 7} fs. At 200 ps, `dt` 4, 5 and 6 pass and
`dt` 7 fails on density (0.59 %) and RDF (0.0295). `dt` 6 is the fastest pass at
2.80x but is close to both limits; `dt` 5 is the safer choice. See Results.

## Timing policy

The machine is shared. A run is timed only when the work of OTHER users and jobs
is quiet. The rule subtracts the sweep's own OpenMP threads from the 1-minute
load average first, then requires the remainder to be at or below 2.0. Without
that subtraction the sweep's own four threads would push the load above 2.0 and
mark every run after the first as untimed, which is not the intent of the rule.

If any contributing run was untimed, the report withholds the speedup and prints
every performance number as untimed. The winner is still selected by measured
ns/day, but a loaded number is never presented as a result.

## Assumptions

- Every run covers the same physical time as the reference, so a larger
  timestep does not simply shorten the run.
- The thermostat is `v-rescale` with `tau-t = 0.1`, matching
  `bench/base.mdp`. At `tau-t = 0.5`, `grompp` treats the
  `verlet-buffer-tolerance = 0.05` warning as fatal, which would make that
  sweep point untestable.
- The pressure coupling is `C-rescale` at 1 bar, isotropic.
- The reference itself carries finite-run statistical error. The window is
  sized so that error sits below the thresholds; the measured floors are in
  the calibration table below.

## Verification noise floor

`selfcheck` measures the smallest difference the four checks can resolve. Three
independent 2 fs samples, 200 ps each, one velocity seed per sample:

| system      | atoms | temperature | density | RDF max dev |
|-------------|-------|-------------|---------|-------------|
| water small | 2652  | 0.931 K     | 0.172 % | 0.0059      |
| villin      | 14528 | 0.722 K     | 0.103 % | 0.0064      |
| check limit |       | 1.0 K       | 0.5 %   | 0.02        |

At 200 ps every floor sits below its limit, so a correct candidate passes and a
wrong one fails. The temperature floor is the tightest at 0.93 K against a 1 K
limit. A single 200 ps sample can therefore fail the temperature check by
chance, so a marginal temperature failure is not proof of a bad setting. Read
`selfcheck.txt` for the full table. The floors fall as the run gets longer.

Before this calibration, 50 ps runs gave an RDF floor near 0.018 against the
0.02 limit, which is why short runs are not trusted for an RDF verdict.

## Results

Timed on an idle 10-core machine, 4 OpenMP threads, 200 ps per run, every run
timed (other load 0.00). Full tables are in each `out/<system>/report.txt`.

| system      | atoms | fastest passing                        | speedup | drift | RDF dev | density |
|-------------|-------|----------------------------------------|---------|-------|---------|---------|
| water small | 2652  | dt 5 fs, no HMR                        | 2.34x   | 0.0031 | 0.0124 | 0.04 %  |
| water med.  | 8916  | dt 5 fs, no HMR                        | 2.40x   | 0.0011 | 0.0122 | 0.05 %  |
| water large | 21087 | dt 5 fs, no HMR                        | 2.36x   | 0.0007 | 0.0156 | 0.20 %  |
| villin      | 14528 | dt 4 fs, no HMR                        | 1.94x   | 0.0011 | 0.0117 | 0.00 %  |
| villin vsite| 14571 | dt 6 fs, all-bonds                     | 2.80x   | -0.0011| 0.0189 | 0.40 %  |

The three plain-water systems all land on `dt 5 fs`, no HMR, no MTS,
`nstlist 10`, tolerance 0.005, for about 2.3 to 2.4x. The villin protein stops
at `dt 4 fs`, 1.94x, because `grompp` refuses `dt 5 fs` on the CG-OD1 carboxylate
bond and HMR cannot move it. The virtual-site rebuild clears that block and
reaches `dt 6 fs`, 2.80x.

This audit sweeps `constraints` only as `h-bonds`. The brute-force search below
adds `constraints = all-bonds` and reaches `dt 6 fs` at 2.71x with no virtual
sites. The two virtual-site figures differ because the stage-2 sweep and the
later search used different references (2.80x and 2.74x respectively).

The virtual-site result is the headline for proteins, but it is marginal: the
`dt 6` RDF deviation is 0.0189 against the 0.02 limit, and its density is 0.40 %
against 0.5 %. `dt 7` fails both (RDF 0.0295, density 0.59 %). Given the villin
RDF floor of 0.0064, the `dt 6` deviation is real and close to the limit. Treat
`dt 5 fs` on the virtual-site system as the safer choice; it measures 2.37x at
RDF 0.0164 and density 0.13 %.

Failures worth recording, each with its number:

- HMR on a protein does not speed a run at a fixed `dt`, and `grompp` rejects
  `dt 4` and `dt 5` with HMR on (villin, "Too many warnings"). The transform is
  correct but cannot exceed `dt 3`.
- MTS with PME every 3rd step is unstable: drift about 0.51 to 0.53
  kJ/mol/ps/atom, temperature +4.5 to +5.3 K, RDF 0.028 to 0.030, every water
  system. Rejected.
- MTS is never faster here: the best `mts 2` water run is 661 ns/day against 721
  for plain `dt 5`.
- `nstlist 100` on the small box is refused by `grompp`: the pair-list cut-off
  (1.516 nm) exceeds half the shortest box vector (1.4991 nm).
- Noise-sized failures appear at the thresholds: water small `dt 3` fails density
  at 0.55 % (floor 0.172 %), and water small `dt 5` with tolerance 0.05 fails
  temperature at 298.90 K (floor 0.931 K). These are single samples, not proof
  of a bad setting.

The timed sweep over water small, medium and large, villin, and the
virtual-site system is complete. `run_timed.sh` reproduces it.

## Brute-force search (`search.py`)

A two-stage funnel searches the integration-setting space: screen many settings
at 20 ps for stability (drift, temperature, density), then certify the top ones
at 200 ps against all four checks. `search_dashboard.py` shows it live on port
8780. State is in `<out>/search.json` and resumes after a stop.

Villin search result (126 screen settings, 24 certified, 23 pass). The reference
is the 2 fs run.

| setting | ns/day | speedup | drift | RDF dev | density |
|---|---|---|---|---|---|
| dt 7 fs, HMR factor 4.5, all-bonds | 141.1 | 3.01x | -0.0147 | 0.0173 | 0.11 % |
| dt 6 fs, HMR on, all-bonds | 127.1 | 2.71x | -0.0061 | 0.0100 | 0.04 % |
| dt 7 fs, HMR on, all-bonds | 137.8 | 2.94x | -0.0106 | 0.0152 | 0.09 % |

The key lever is `constraints = all-bonds`, not HMR: it removes the stiff heavy-atom
bond that stops the plain protein at 4 fs, so the plain topology reaches 6 to 7 fs.
The `dt 7` peak is marginal (drift and RDF both near their limits). `all-angles`
fails at step 0 on every setting (LINCS constraint error; protein rings make the
constraint graph unsolvable), and 8 to 10 fs all fail. So 7 fs is the ceiling.

Protein-observable validation (`valprotein.py`, single 200 ps runs):

| run | Rg nm | backbone RMSD nm | RMSF corr | RMSF max dev |
|---|---|---|---|---|
| 2 fs reference | 0.9535 | 0.1290 | 1.000 | 0.0000 |
| dt 6, all-bonds | 0.9649 | 0.0803 | 0.903 | 0.0386 |
| dt 7, HMR 4.5, all-bonds | 0.9596 | 0.0822 | 0.801 | 0.0676 |
| dt 4, HMR 3, h-bonds | 0.9648 | 0.0863 | 0.900 | 0.0374 |

`dt 6` all-bonds reproduces the reference RMSF (r = 0.903). `dt 7` is less faithful
(r = 0.801, max deviation 0.068 nm). The ~1 % Rg offset is shared by the h-bonds
`dt 4` candidate, so it is run-to-run variation, not a constraint effect. Single
200 ps runs lack the power to separate constraint effects from thermal noise; the
safe pick is `dt 6 fs` all-bonds at 2.71x, not the marginal `dt 7` at 3.01x.

## Virtual-site search (`search_vsites/`)

The same funnel on the virtual-site rebuild (`vsites/build/eq.gro`), to test
whether virtual sites and all-bonds stack. Reference 59.3 ns/day (200 ps
certify).

| setting | ns/day | speedup | drift | RDF dev | pass |
|---|---|---|---|---|---|
| dt 7, all-bonds | 190.1 | 3.21x | 0.0001 | 0.0299 | FAIL (RDF) |
| dt 6, all-bonds | 162.1 | 2.74x | -0.0011 | 0.0167 | PASS |
| dt 5, all-bonds | 138.9 | 2.34x | | 0.0109 | PASS |
| dt 4, all-bonds | 113.7 | 1.92x | | 0.0069 | PASS |
| dt 4, h-bonds | 116.4 | 1.96x | 0.0016 | | PASS |

They do **not** stack. The best passing setting is `dt 6` all-bonds at 2.74x,
which is the same as the plain all-bonds route (2.71x). Virtual sites remove the
hydrogen degrees of freedom; they do not remove the heavy-atom bond that sets the
7 fs wall. Above dt 4 the h-bonds topology is still refused by grompp, and dt 8
and up crash. The `dt 7` setting gains speed but its RDF fails (0.0299).

The search space of this engine is now bounded: the integration settings
(dt, HMR, constraints, mass scale) top out near 3x, and the safe setting is
about 2.6 to 2.7x. More search on these axes will not reach 10x.

## Spherical solvent shell and reaction field (`spherical.py`)

This is a model change, not an integration setting. It keeps the protein plus a
mobile water shell, freezes the outer water, and puts the droplet in a larger
box. It removes the bulk water, which is most of the atoms. `spherical.py build`
makes the system; `spherical.py run` runs it. `--coulomb reaction-field` replaces
PME. Compare every variant against the periodic 2 fs reference (46.8 ns/day) with
`valprotein.py`.

| variant | atoms | dt | ns/day | vs periodic | RMSF corr | speedup realistic? |
|---|---|---|---|---|---|---|
| periodic all-bonds | 14528 | 6 fs | 127.1 | 2.71x | 0.903 | yes, no model change |
| shell 0.8 freeze 0.25, PME | 3011 | 2 fs | 205.9 | 4.4x | 0.650 | no, ensemble changed |
| shell 0.8 freeze 0.25, PME | 3011 | 4 fs | 366.0 | 7.8x | 0.810 vs its dt2 | no, ensemble changed |
| shell 0.9 freeze 0.1, PME | 3506 | 2 fs | 159.6 | 3.4x | 0.914 | marginal |
| shell 0.9 freeze 0.1, PME | 3506 | 4 fs | 266.3 | 5.7x | 0.960 vs its dt2 | marginal |
| shell 0.9 freeze 0.1, RF | 3506 | 2 fs | 431.4 | 9.2x | 0.793 | no, ensemble changed |
| shell 0.9 freeze 0.1, RF | 3506 | 4 fs | 605.8 | 12.9x | 0.918 vs its dt2 | no, ensemble changed |

Key points:

- The first shell must stay mobile. Freezing at 0.55 nm (v1) drops the RMSF
  correlation to 0.650; freezing only beyond 0.8 nm (v2) recovers it to 0.914.
- The shell ceiling is dt 4 fs. dt 5 and 6 crash with a SETTLE error.
- Reaction field is a further model change. It is fast (it removes the expensive
  PME grid of the large box) but its RMSF correlation against the periodic
  reference drops to 0.793.
- The honest picks: 2.71x with no model change, or 5.7x with a spherical shell at
  a small measured ensemble cost. The 9 to 13x figures include a reaction-field
  model change and need their own validation before use.

Single 200 ps runs have limited statistical power. The correlations are
indicative, not certified.

## Soft boundary (`softshell.py`)

The frozen outer water is a rigid 0 K wall. `softshell.py build` replaces it with
a restrained water moleculetype (`BWA`, harmonic position restraints, in a local
`boundary.itp`), so the boundary water is thermostatted and can fluctuate. Same
size and geometry as the v3 shell: 974 mobile, 627 restrained, 5387 atoms, box
8 nm. Run with `softshell.py run --coulomb reaction-field`.

| variant | dt | ns/day | vs periodic | RMSF corr | mobile q | q std |
|---|---|---|---|---|---|---|
| frozen shell v3, RF | 4 fs | 605.8 | 12.9x | 0.918 | 0.403 | 0.314 |
| soft shell, RF | 4 fs | 325.8 | 7.0x | 0.922 | 0.469 | 0.261 |
| periodic bulk reference | 2 fs | 46.8 | 1.0x | 1.000 | 0.557 | 0.202 |

Key points:

- The soft boundary fixes the unphysical 0 K wall. The mobile water moves toward
  bulk order (tetrahedral q 0.403 -> 0.469, toward the bulk 0.557) with a narrower
  distribution (std 0.314 -> 0.261). This is the validation improvement.
- The soft boundary does NOT raise the timestep ceiling. dt 5 still fails with a
  SETTLE error, at restraint stiffness k 1000 and k 200, with both reaction field
  and PME. The limiter is the confined mobile water at the protein surface, not
  the wall. dt 4 remains the shell ceiling.
- It is slower (7.0x vs 12.9x) because it keeps 974 mobile waters instead of 809
  and uses rcoulomb 1.2 nm.
