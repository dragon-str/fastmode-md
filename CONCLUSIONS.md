# Fast-mode MD audit: findings and recommendations

## Goal

Audit GROMACS integration settings and validate every candidate against a 2 fs
reference of the same system. The settings are timestep, hydrogen mass
repartitioning (HMR), multiple time stepping (MTS), neighbour-list frequency,
Verlet buffer tolerance, constraints, and virtual sites. The physics must be
checked, not assumed. Four checks must pass against a measured noise floor:
conserved energy drift, temperature, O-O RDF, and density. Every failure is
reported with its number, and no threshold is lowered.

## What the tool provides

- `fastmode.py` - the `audit` tool, plus `hmr-check` and `selfcheck` (noise
  floor).
- `stage2.py` - virtual-site rebuild and sweep.
- `search.py` + `search_dashboard.py` - a two-stage search over the
  integration-settings space, with a live dashboard.
- `spherical.py`, `softshell.py` - solvent-boundary variants.
- `validate10x.py`, `valprotein.py`, `dashboard.py`, `run_timed.sh`.
- `ml/` - a documented set of machine-learning probes, most of them negative.
- READMEs that record every number and every failure.

## Speed results (this CPU, 4 OpenMP threads, ARM Mac)

| system | setting | speedup | model change |
|---|---|---|---|
| water small/medium/large | dt 5 fs, no HMR | 2.34x / 2.40x / 2.36x | no |
| villin HP-35 | dt 6 fs, `constraints=all-bonds` | 2.71x | no |
| villin HP-35 | dt 7 fs, HMR f4.5, all-bonds | 3.01x | no, marginal |
| villin, spherical shell + reaction field | dt 4 fs, all-bonds | 9.9-12.9x | yes |
| villin, soft boundary shell | dt 4 fs, all-bonds | 7.0x | yes |

The drop-in setting (2.71x) is the safest useful result. It needs no model
change, and it passed the four checks plus a protein-observable comparison
(backbone RMSF correlation 0.90). The 3.01x dt 7 setting is marginal on drift and
RDF and is not recommended.

The 9.9-12.9x shell results are real, but they are a model change: a droplet in a
periodic box with a frozen (or restrained) outer water shell and a reaction
field. They reproduce the protein ensemble well (Rg within 0.4 percent, RMSF
correlation 0.97-0.99 at dt 4), but the mobile water is interfacial, not bulk
(tetrahedral order 0.40-0.47 against a bulk 0.557), and in the frozen variant the
protein approaches the 0 K wall. They are validated for protein observables, not
for solvent thermodynamics.

## Where the limits are

1. The protein timestep limit is the carboxylate CG-OD1 bond (oscillation period
   0.022 ps), not a hydrogen bond. HMR cannot change a C-O bond, so HMR alone
   caps the speedup near 1.5x.
2. `constraints = all-bonds` is the real lever. It removes the grompp warning and
   lets the plain protein reach dt 6-7 fs, which gives the 2.71x.
3. Virtual sites do not stack with all-bonds. A virtual-site rebuild reaches the
   same 2.7x as the plain all-bonds route, for much more work.
4. MTS with PME every third step is unstable on water (drift 0.5, not 0.02). The
   stability boundary is a resonance at a 12 fs slow-force interval: every
   setting at or below 12 fs passes, and every setting at 15 fs or more fails.
5. The shell timestep ceiling is 4 fs. It is set by the confined surface water,
   not by the boundary stiffness or the electrostatics (dt 5 fails in every
   combination tested).
6. `all-angles` constraints are impossible here. Protein rings make the LINCS
   constraint graph unsolvable, so the run fails at step 0.

## Negative machine-learning results

All probes were GPU-free and cheap.

| route | result |
|---|---|
| learned next-position propagator | closed: Lyapunov time 0.12 ps; a 200 ps path needs 1e-70 accuracy |
| position-only force model | closed: R^2 0.81; the residual is PME plus the constraint force |
| learned MTS correction | closed: the instability is a resonance, not a bias |
| local water-patch proposal | closed: the Hastings penalty exceeds the gain |
| Boltzmann generator (Cartesian) | closed: unphysical bond geometry, U ~ 1e18 |
| Boltzmann generator (rigid-body) | closed: never learns excluded volume, acceptance 0 |
| one-site CG water (pair + 3-body IBI) | closed: first peak 25 percent too low |
| one-site CG water (pair only) | closed: IBI diverges |

An ML surrogate did not help on this hardware, because a neural network is slower
than optimized C per force evaluation. A learned method can only win by reducing
work or by generating samples. Every sample-generating route failed on the same
obstacle: the sharp, high-dimensional structure of a dense liquid (stiff bonds,
excluded volume, and the first RDF peak). This is a clean, reusable negative
result.

## Scope and novelty

Every speed ingredient here is standard: HMR, virtual sites, all-bonds
constraints, reaction field, and spherical boundaries. The 10x is therefore a
stack of known approximations, not a new method. The contribution is the
methodology: a systematic integration-space search with four observable checks
against a measured noise floor, the quantitative result that the carboxylate bond
rather than hydrogen sets the protein timestep, the negative result that virtual
sites do not stack with all-bonds, and a documented set of closed ML routes with
numbers.

## Recommendations

Ship the drop-in 2.71x setting. Keep the shell plus reaction field as an
explicitly labelled model change for protein-only questions. Stop the GPU-free ML
search; the remaining honest path to a general 10x is a GPU or a research-grade
coarse-grained model with a proper force model.
