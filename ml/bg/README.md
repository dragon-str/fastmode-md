# Boltzmann generator prototype (`bg/`)

Goal: generate decorrelated equilibrium samples of water without stepping MD,
using a normalizing flow as an independence proposal in Metropolis-Hastings,
scored by the exact GROMACS potential. Equilibrium properties only.

## Setup

- `build/` - 128 flexible TIP3P waters, box 1.565 nm, 384 atoms. `em.mdp`,
  `prod.mdp` (dt 0.5 fs, `-DFLEXIBLE`, no constraints), 500 ps reference
  (`prod.xtc`, 1001 frames).
- `small/` - 31 flexible TIP3P waters, box 0.98 nm, 93 atoms, 4 ns reference
  (4001 frames). MD reference potential: mean -1193.4, std 39.8 kJ/mol.
- `bg.py` - representation (min-image displacement from a fixed anchor), RealNVP
  (affine coupling, tanh-capped scale), training by maximum likelihood, and a
  batched Metropolis independence sampler using GROMACS `-rerun` energies.
  CLI: `bg.py --build <dir> --frames <traj.gro> --work <dir> [--epochs] [--clip]
  [--load]`.

## Result (negative, with diagnosis)

Both systems failed. After 400 epochs the flow fits the training set
(`nll` near 0) but its samples give potential energies of order 1e10-1e20 and
acceptance is exactly 0. A round-trip test (`forward` then `forward_inv`) has
error 3e-6, so the inverse is correct; this is not an implementation bug.

The cause is representation, not capacity. The flow samples Cartesian
coordinates, so hydrogens land at unphysical distances from their oxygens. The
flexible O-H bond has k = 502416 kJ/mol/nm^2, so a 0.05 nm error costs ~600
kJ/mol and a 0.5 nm error costs ~6e4 kJ/mol; the flow's tails reach several nm,
so the proposal energy explodes. A generic Cartesian flow cannot represent the
stiff bond geometry that MD enforces with constraints.

## Rigid-body representation (`bgrigid.py`)

Each water is an oxygen position (3) plus an orientation rotation vector (3);
hydrogens are rebuilt from the fixed TIP3P geometry, so bond geometry is exact
by construction (verified: O-H 0.095720, H-O-H 104.520, O-O relative error
6e-17). The orientation Jacobian `J(omega) = 2(1-cos|omega|)/|omega|^2` enters
the Metropolis ratio. Reference: `rigid/`, 31 rigid TIP3P waters, box 0.98 nm,
4 ns at dt 2 fs (`rigid.xtc`, 976 training frames after a 100 ps warm-up).

Still fails. After 600 epochs the flow fits loosely (nll -80) but acceptance is
0/2000 and proposal energies are ~1e21 kJ/mol. The diagnostic:

| quantity | MD frames | flow samples |
|---|---|---|
| min O-O distance (nm) | 0.257 mean, 0.242 min | 0.039 mean, 0.004 min |
| mean log q | 88.2 | 9.4 |

The flow does not learn the hard-core O-O repulsion: it proposes oxygens at
0.04 nm where the liquid keeps 0.26 nm, and it assigns *higher* density to those
overlaps (9.4) than to real liquid frames (88.2). This is the real wall.

Two reasons, both fundamental for this setting. First, an independence sampler's
acceptance falls exponentially with dimension; the proposal q must be extremely
close to the target p, and a RealNVP trained on ~1000 frames in 183 dimensions
is nowhere close. Second, the liquid's excluded volume is a sharp constraint
that a generic density model does not learn from so little data.

## Verdict

A GPU-free flow-based Boltzmann generator does not work for a 31-water liquid.
The rescues do not help: a system small enough to sample (a few waters) gives no
useful speedup, and conditioning the flow on a local environment is exactly the
water-patch route already closed (Experiments 10-12). The remaining credible
GPU-free ML routes are a learned 3-body coarse-grain potential, learned
collective variables for enhanced sampling, and large-timestep reweighting.

