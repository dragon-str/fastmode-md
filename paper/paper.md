---
title: "Auditing fast-mode molecular dynamics settings without fooling yourself"
author:
  - name: Daniel Reda
    orcid: 0000-0001-6160-6215
date: 2026-09-12
---

# Abstract

Integration settings can make a molecular dynamics simulation several times
faster, but the fastest setting that runs is often the wrong one. We present an
open-source audit tool that searches the integration-settings space and validates
every candidate against a 2 fs reference of the same system with four physics
checks: conserved energy drift, temperature, the oxygen-oxygen radial
distribution function, and density. The tool measures its own noise floor, so a
check is only reported as passed when the difference is larger than the
measurement error, and it reports every failure with its number. We use it to map
the speed limit on a commodity CPU for pure water and for the villin HP-35
domain. The safe, drop-in setting is a 2.7x speedup with no model change. A
spherical solvent shell with a reaction field reaches 7-13x but is a model change
that reproduces protein observables while leaving the mobile water interfacial.
We identify the physics that sets each limit: the protein timestep is bounded by
a carboxylate carbon-oxygen bond, not by hydrogen, so hydrogen mass repartitioning
caps near 1.5x; multiple time stepping has a stability resonance at a 12 fs
slow-force interval; and the shell ceiling is set by confined surface water. We
also report a set of machine-learned speedup routes that we tested and closed,
each with its measured failure mode. The contribution is methodological: a way to
choose fast-mode settings that does not mistake a broken simulation for a fast
one.

# Statement of need

Molecular dynamics users routinely raise the timestep, repartition hydrogen
masses, switch on multiple time stepping, and shorten neighbour-list updates to
gain speed. Each of these changes can silently break the physics. The failure is
hard to see because a broken run still produces a trajectory, and because the
difference between a thermal fluctuation and an integration error is often
smaller than the run-to-run scatter. Software packages warn about some changes
(for example, GROMACS refuses a timestep that exceeds a bond's oscillation
limit), but they do not tell the user whether the settings that did run are
correct.

The audit tool addresses that gap. It makes the comparison explicit. It runs a
reference at 2 fs, then runs candidates that vary the timestep, hydrogen mass
repartitioning (HMR), multiple time stepping (MTS), the neighbour-list update
frequency, the Verlet buffer tolerance, the constraint set, and virtual sites. It
applies four checks against the reference. It measures the noise floor of each
check with independent replicas before it trusts a verdict. It never lowers a
threshold, and it records grompp refusals as such rather than forcing them with
`-maxwarn`. The result is a ranked, evidence-backed list of settings instead of a
single unverified number.

# Method

## Systems

We use four water boxes and one protein. The water boxes contain 2652, 8916 and
21087 atoms (TIP3P at three densities) and the villin HP-35 domain
contains 14528 atoms (582 protein atoms, 4648 waters, 2 chloride ions). We also
build a virtual-site variant of the protein with `pdb2gmx -vsite h -heavyh`.
Runs use GROMACS 2027.0-dev on an ARM CPU with four OpenMP threads and no GPU.

## The four checks

A candidate must satisfy all of the following against the 2 fs reference of the
same system, evaluated over the same production window:

1. Conserved energy drift below 0.02 kJ/mol/ps per atom, read from the GROMACS
   log line `Conserved energy drift`.
2. Mean temperature within 1 K of the reference mean.
3. Maximum absolute deviation of the O-O radial distribution function below 0.02.
4. Mean density within 0.5 percent of the reference.

A missing check is a failure, never a pass. Every failure is reported with its
measured number.

## Noise floor

A threshold is meaningless without the measurement scatter. The tool's
`selfcheck` mode runs three independent replicas of the reference with different
velocities and compares the three results. At the production length of 200 ps the
floors are:

| system | temperature (K) | density (percent) | O-O RDF |
|---|---|---|---|
| water, 2652 atoms | 0.931 | 0.172 | 0.0059 |
| villin, 14528 atoms | 0.722 | 0.103 | 0.0064 |

The RDF floor is the binding constraint. At 50 ps it rises to about 0.018, close
to the 0.02 limit, which is why the default production window is 200 ps. The
checks are applied against these measured floors, and a difference comparable to
the floor is reported as unresolved rather than as a pass or a failure.

## Timing policy

A run is timed only when the external machine load is below 2. External load is
the one-minute load average minus the tool's own thread count. Otherwise the run
is recorded as untimed and no speedup is reported, because a shared machine can
change a timing by more than the effect being measured.

# Results

## Drop-in speedups

The fastest settings that pass all four checks with no model change are:

| system | setting | speedup |
|---|---|---|
| water, 2652 atoms | dt 5 fs, MTS off, nstlist 10 | 2.34x |
| water, 8916 atoms | dt 5 fs, MTS off, nstlist 10 | 2.40x |
| water, 21087 atoms | dt 5 fs, MTS off, nstlist 10 | 2.36x |
| villin HP-35 | dt 6 fs, `constraints = all-bonds` | 2.71x |
| villin HP-35, virtual sites | dt 6 fs, `constraints = all-bonds` | 2.74x |

The villin drop-in winner also passed a protein-observable check: its backbone
root-mean-square fluctuation profile correlates with the reference at 0.90, and
its radius of gyration agrees to within about 1 percent.

## Where the limits are

**Hydrogen is not the protein limit.** HMR alone cannot exceed dt 3 fs on the
protein: grompp flags the carboxylate CG-OD1 bond with an oscillation period of
0.022 ps. HMR does not change a carbon-oxygen reduced mass, so it caps the
speedup near 1.5x. The real lever is `constraints = all-bonds`, which removes the
warning and lets the plain protein reach dt 6-7 fs.

**Virtual sites do not stack with all-bonds.** A rebuilt virtual-site system
reaches the same speed as the plain all-bonds route (2.74x against 2.71x) for
much more setup work. The virtual-site dt 7 candidate fails the RDF check
(deviation 0.0299). The reason is the same carbon-oxygen bond: virtual sites
remove hydrogen degrees of freedom, not the heavy-atom bond that sets the wall.

**All-angles constraints are unusable.** Every all-angles candidate failed at
step 0 with a LINCS constraint error, independent of timestep and of HMR. The
protein rings make the constraint graph unsolvable. This is structural, not a
timestep problem.

**Multiple time stepping has a resonance.** Evaluating the long-range force
every second step at a 4-5 fs base timestep can pass, but at a slow-force
interval of 15 fs or more it fails or crashes; every interval at or below 12 fs
passes. MTS never beat the plain dt 5 setting in our sweep, and the PME-every-3
variant was unstable on every water box (drift about 0.5 against a 0.02 limit,
temperature 4.5-5.3 K high).

# The solvent shell: a model change, not a speedup

Removing bulk water and keeping only a thin solvent shell gives the largest
number. A spherical droplet of 974 mobile waters inside a box, with an outer
water layer restrained or frozen and a reaction field for the electrostatics,
runs at dt 4 fs and 7-13x the periodic reference depending on the boundary
variant. We validated it with up to five replicas of 0.3 ns against a periodic 2 fs
reference:

| configuration | Rg (nm) | RMSF corr. | O-O order q | drift |
|---|---|---|---|---|
| periodic, dt 2 fs | 0.9614 | 1.000 | 0.557 | 0.0008 |
| shell, PME, dt 4 fs (n=3) | 0.9627 | 0.969 | 0.403 | -0.0126 |
| shell, reaction field 1.2 nm, dt 4 fs | 0.9579 | 0.988 | 0.406 | -0.0119 |
| shell, reaction field 1.5 nm, dt 4 fs | 0.9583 | 0.968 | 0.405 | -0.0120 |

A soft (position-restrained rather than frozen) boundary was validated with a
single 200 ps run: Rg 0.9644, RMSF correlation 0.922, mobile-water order 0.469.

The shell reproduces the protein ensemble well: the radius of gyration agrees to
within 0.4 percent and the RMSF correlation is 0.97-0.99. But it is a different
model. The mobile water is interfacial, not bulk (tetrahedral order 0.40-0.47
against a bulk 0.557), and in the frozen variant the protein approaches the
boundary (closest approach 0.28-0.30 nm, and up to 8.5 percent of frames within
0.35 nm for the largest cutoff). Replacing the frozen wall with a
position-restrained boundary moves the water toward bulk order (q 0.469) and
removes the unphysical zero-temperature wall, but it does not raise the 4 fs
ceiling and it is slower (7.0x).

The shell result is therefore valid for protein observables and is not valid for
bulk-solvent thermodynamics. The headline number must be reported with that
qualifier.

# Machine-learned speedups: a documented negative result

We tested every GPU-free machine-learning route we could identify that would give
a general speedup, and we closed each one with a measurement:

| route | outcome |
|---|---|
| learned next-position propagator | closed: Lyapunov time 0.12 ps; a 200 ps trajectory needs 1e-70 accuracy |
| position-only force model | closed: R^2 0.81; the residual is PME plus the constraint force |
| learned MTS splitting | closed: the limit is a 12 fs resonance, not a smooth bias |
| local water-patch proposal | closed: the Hastings penalty exceeds the gain |
| Boltzmann generator, Cartesian | closed: unphysical bond geometry; proposal energy about 1e18 |
| Boltzmann generator, rigid body | closed: the flow never learns excluded volume; acceptance 0/2000 |
| one-site coarse-grained water | closed: iterative Boltzmann inversion diverges; the first RDF peak is 25 percent too low |

The unifying reason is that an ML surrogate does not help on this hardware: a
neural network is slower than optimized C per force evaluation. ML can only win
by reducing work or by generating samples. Every sample generator we tried failed
on the same obstacle, the sharp high-dimensional structure of a dense liquid:
stiff bonds, excluded volume, and the first RDF peak. This is a useful and
reusable negative result, and it sharpens the remaining options to a GPU or a
research-grade coarse-grained force model.

# Discussion

The value of this work is methodological, and it is worth stating the limits
plainly. Every speed ingredient we use is standard: hydrogen mass
repartitioning, virtual sites, all-bonds constraints, reaction field, and
spherical solvent boundaries. The 2.7x drop-in setting is an assembly of known
options, and the 7-13x shell result is a stack of approximations rather than a new
method. We did not discover a new algorithm, and we did not reach a general 10x
without changing the model.

What is new is the audit discipline. Four checks, a measured noise floor, and a
timing policy turn a search over integration settings into an evidence-backed
recommendation. The same discipline produced two results that we would not have
trusted from a single run: that the protein timestep is set by a carboxylate
bond rather than by hydrogen, and that virtual sites do not stack with all-bonds.
It also produced a clean set of closed machine-learning routes with numbers
rather than opinions.

# Conclusion

For a small protein on a commodity CPU, the recommended setting is dt 6 fs with
`constraints = all-bonds`, which gives about 2.7x with no model change and passes
the four checks plus a protein-observable comparison. A spherical shell with a
reaction field can reach 7-13x, but it must be labelled a model change and
validated for the observables of interest. Hydrogen mass repartitioning alone is
not worth much, multiple time stepping is not worth the risk at these timesteps,
and a GPU-free learned surrogate is not a path to speed here.

# Data and code availability

The tool, the search harness, the validation scripts, and the reduced results
are at `https://github.com/dragon-str/fastmode-md` under the MIT license. The
trajectories are not stored because every reported number is reduced to a report
file in `results/` and every run is reproducible from the tool. The archived
version is registered on OSF at `https://osf.io/x4t8m/` (DOI
`10.17605/OSF.IO/X4T8M`).

# Acknowledgements

The author thanks the GROMACS developers for the simulation engine.
