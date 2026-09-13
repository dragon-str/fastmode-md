# Final conclusions - fast-mode MD audit

## What was asked

The handoff asked for a tool that searches GROMACS integration settings (timestep,
hydrogen mass repartitioning, multiple time stepping, neighbour-list frequency,
Verlet buffer tolerance, virtual sites) and validates every candidate against a
2 fs reference of the same system. The physics must be checked, not assumed:
conserved energy drift, temperature, O-O RDF, and density, each against a
measured noise floor, with every failure reported and no threshold lowered.

## What was delivered

- `fastmode.py` - the `audit` tool plus `hmr-check` and `selfcheck` (noise floor).
- `stage2.py` - virtual-site rebuild and sweep.
- `search.py` + `search_dashboard.py` - a two-stage brute-force funnel over the
  integration-settings space, with a live dashboard.
- `spherical.py`, `softshell.py` - solvent-boundary variants.
- `validate10x.py`, `valprotein.py`, `dashboard.py`, `run_timed.sh`.
- `ml/` - a complete, documented set of ML probes, most of them negative.
- READMEs that record every number and every failure.

## Speed results (this CPU, 4 OpenMP threads, ARM Mac)

| system | setting | speedup | model change |
|---|---|---|---|
| water small/medium/large | dt 5 fs, no HMR | 2.34x / 2.40x / 2.36x | no |
| villin HP-35 | dt 6 fs, `constraints=all-bonds` | 2.71x | no |
| villin HP-35 | dt 7 fs, HMR f4.5, all-bonds | 3.01x | no, marginal |
| villin, spherical shell + reaction field | dt 4 fs, all-bonds | 9.9-12.9x | yes |
| villin, soft boundary shell | dt 4 fs, all-bonds | 7.0x | yes |

The drop-in setting (2.71x) is the safest useful result. It needs no model change
and it passed the four checks plus a protein-observable comparison (backbone RMSF
correlation 0.90). The 3.01x dt 7 setting is marginal on drift and RDF and is not
recommended.

The 9.9-12.9x shell results are real but they are a model change: a droplet in a
periodic box with a frozen (or restrained) outer water shell and a reaction field.
They reproduce the protein ensemble well (Rg within 0.4%, RMSF correlation
0.97-0.99 at dt 4), but the mobile water is interfacial, not bulk (tetrahedral
order 0.40-0.47 against a bulk 0.557), and in the frozen variant the protein
approaches the 0 K wall. Validated for protein observables, not for solvent
thermodynamics.

## Physics findings (the useful part)

1. The protein timestep limit is the carboxylate CG-OD1 bond (oscillation period
   0.022 ps), not a hydrogen bond. HMR cannot change a C-O bond, so HMR alone
   caps the speedup near 1.5x.
2. `constraints = all-bonds` is the real lever: it removes the warning and lets
   the plain protein reach dt 6-7 fs. This is what gives the 2.71x.
3. Virtual sites do NOT stack with all-bonds. A virtual-site rebuild reaches the
   same 2.7x as the plain all-bonds route, for much more work.
4. MTS with PME every 3rd step is unstable on water (drift 0.5, not 0.02). The
   stability boundary is a resonance at a 12 fs slow-force interval: everything
   at or below 12 fs passes, everything at 15 fs or more fails.
5. The shell timestep ceiling is 4 fs, set by the confined surface water, not by
   the boundary stiffness or the electrostatics (dt 5 fails at every combination
   tested).
6. `all-angles` constraints are impossible here (protein rings make the LINCS
   graph unsolvable).

## ML exploration (all GPU-free, all negative, all cheap)

| route | result |
|---|---|
| learned next-position propagator | closed: Lyapunov time 0.12 ps; a 200 ps path needs 1e-70 accuracy |
| position-only force model | closed: R^2 0.81; residual is PME plus the constraint force |
| learned MTS correction | closed: the instability is a resonance, not a bias |
| local water-patch proposal | closed: the Hastings penalty exceeds the gain |
| Boltzmann generator (Cartesian) | closed: unphysical bond geometry, U ~ 1e18 |
| Boltzmann generator (rigid-body) | closed: never learns excluded volume, acceptance 0 |
| one-site CG water (pair + 3-body IBI) | closed: first peak 25% too low |
| one-site CG water (pair only) | closed: IBI diverges |

The unifying lesson: an ML surrogate does not help on this hardware because a
neural network is slower than optimized C per force evaluation. ML can only win
by reducing work or by generating samples, and every sample-generating route we
tried fails on the same obstacle - the sharp, high-dimensional structure of a
dense liquid (stiff bonds, excluded volume, the first RDF peak). This is a
clean, reusable negative result.

## Did we achieve anything meaningful?

Yes, as engineering, not as new physics.

- A working, validated audit tool with a measured noise floor. It rejects unsafe
  settings and never reports a speedup it cannot back with physics.
- An honest 2.71x drop-in speedup for a small protein, and a validated
  (model-change) path to about 10x for protein observables.
- A precise map of where the limits are: the CG-OD1 bond, the 12 fs MTS
  resonance, the 4 fs shell ceiling.
- A set of ML dead ends ruled out with measurements rather than opinion.

We did not find a new algorithm or a new force field, and we did not reach a
general 10x without a model change.

## Is it publishable?

Not as a high-impact result. Every speed ingredient is standard: HMR, virtual
sites, all-bonds constraints, reaction field, spherical boundary. The 10x is a
stack of known approximations, not a discovery.

There is a real but modest paper here, best framed as a methods and
reproducibility note: "how to choose fast-mode MD settings without fooling
yourself." Its contributions would be (1) the systematic integration-space search
with four observable checks against a measured noise floor; (2) the quantitative
result that the carboxylate bond, not hydrogen, sets the protein timestep;
(3) the negative result that virtual sites do not stack with all-bonds;
(4) a documented set of closed ML routes with numbers. Suitable for a
software/tools venue or arXiv, not a physics journal.

## Recommendation

Ship the drop-in 2.71x setting. Keep the shell + reaction field as an explicitly
labelled model change for protein-only questions. Stop the GPU-free ML search;
the remaining honest path to a general 10x is a GPU or a research-grade coarse
grained model with a proper force model.
