# ML fast-mode probe (`phase0/ml/`)

Question: can a model change give ~10x over the classical all-atom run, on this
CPU, without a GPU? This directory holds the small experiments that answer it.

Status: experiments 1 to 3 are done. The verdict is below. No result here changes
any threshold in `phase0/fastmode.py`.

## Experiment 1 - step cost against net cost (`ceiling.py`)

Measured on the water small box (2652 atoms, 4 OpenMP threads):

| item | value |
|---|---|
| GROMACS step at dt 5 fs | 0.68 ms/step (631 ns/day, 300-step mdrun) |
| NumPy MLP forward, hidden 16 | 0.13 ms/call |
| NumPy MLP forward, hidden 64 | 0.35 ms/call |
| NumPy MLP forward, hidden 128 | 0.68 ms/call |
| NumPy MLP forward, hidden 256 | 1.40 ms/call |

The compute is not the blocker. A hidden-128 network costs about one MD step, so
a learned propagator that covers 5 to 20 steps could be faster in principle. Two
caveats: this timing excludes the neighbour list and the descriptor, which a real
model must also build, and it says nothing about accuracy.

## Experiment 2 - step cost against particle count (`ceiling.py`)

Fit of the measured step cost across the four systems: `t_step ~ N^0.873`.

| coarse-grain ratio | predicted ms/step reduction |
|---|---|
| 2:1 | 1.8x |
| 4:1 | 3.3x |
| 8:1 | 6.0x |
| 16:1 | 11.2x |

A 10x particle-reduction gain needs about a 16:1 model. That is very coarse.

## Experiment 3 - force fit (`forcefit.py`)

A rotationally equivariant pairwise model, `F_i = sum_j c(r_ij) rhat_ij`, with a
radial basis expansion of `c`. Trained on 35 frames (1 ps apart) of the water
box, tested on 5 held-out frames. Forces come from `gmx dump` on a 1 ps rerun.

| cutoff (nm) | basis | features | force RMS error | R squared |
|---|---|---|---|---|
| 0.5 | 16 | 48 | 250.9 | 0.815 |
| 0.8 | 24 | 72 | 256.1 | 0.807 |
| 1.0 | 32 | 96 | 255.1 | 0.808 |
| 1.2 | 40 | 120 | 254.7 | 0.809 |

Force magnitude is 582.6 kJ/mol/nm. The fit plateaus at R squared 0.81 for every
cutoff and basis. The missing 19 percent is therefore not radial resolution. It
is the many-body part: three-body angle terms and the SETTLE constraint force,
which no pairwise radial form can represent.

## Experiment 4 - angular descriptor (`threebody.py`)

Added a local anisotropy scalar to the coefficient, so the model becomes
`F_i = sum_j c(r_ij, q_ij) u_ij` with `q_ij = sum_k cos(theta_jik) g(r_ik) / n_i`.
This is a cheap three-body term, and it stays equivariant.

| model | features | force RMS error | R squared |
|---|---|---|---|
| pairwise | 48 | 252.6 | 0.812 |
| pairwise + angular | 384 | 251.5 | 0.814 |

Eight times the features move R squared by 0.002. The residual is therefore not a
lack of angular resolution either.

Decomposition by atom type (pairwise model):

| atoms | force RMS error | magnitude | R squared |
|---|---|---|---|
| all | 252.6 | 582.5 | 0.812 |
| oxygen only | 308.5 | 795.2 | 0.849 |
| hydrogen only | 160.4 | 439.2 | 0.867 |

Even the oxygen forces stop at 0.85. The unrepresented part is not only the
hydrogen constraint force. It is also the long-range PME force, which is not a
function of the local geometry, and the velocity-dependent constraint Lagrange
multiplier. Neither is a function of the instantaneous local positions, so no
local position-only model can reach them.

## Experiment 5 - where can a coarse map be selective (`cgprobe.py`, `cgvillin.py`)

Question: if a model replaces some waters and keeps others atomistic, where is the
boundary? Answer from the two reference trajectories.

Tetrahedral order `q` in bulk water has a correlation length of about one
molecule. The spatial correlation of the `q` deviation is positive only in the
first shell (0.137 at r 0.23 nm, 0.055 at 0.33 nm) and is zero by 0.4 nm. Pure
water has no coherent regions to partition. Order statistics: `q` mean 0.571,
std 0.192; coordination mean 5.20, std 1.17; 12 percent of molecules above 0.8.

Villin perturbs water structure only inside a thin shell. `q` against the
distance to the nearest protein atom:

| distance (nm) | q | std | samples |
|---|---|---|---|
| 0.0-0.2 | 0.371 | 0.202 | 808 |
| 0.2-0.4 | 0.514 | 0.203 | |
| 0.4-0.6 | 0.574 | | |
| 0.6 and out | 0.572 | 0.191 | |

Beyond 0.6 nm the water is statistically bulk. The water fraction beyond a
distance from the protein: 0.942 at 0.4 nm, 0.887 at 0.6, 0.733 at 1.0, 0.455 at
1.5, 0.185 at 2.0. Mean nearest-protein distance is 1.413 nm.

Conclusion: a distance-adaptive map has a thin atomistic shell and a uniform
bulk. Spatial selectivity gains little. The useful levers are site elimination
and a larger timestep, not spatial partitioning.

## Experiment 6 - one-site water by iterative Boltzmann inversion (`cg_ibi.py`)

Standalone NumPy probe, no GROMACS change. One site per water, effective pair
potential fitted by `u <- u + kT ln(g / g_target)`, with the target from the
reference O-O RDF. Langevin velocity-Verlet, 512 particles, dt 0.005 ps, 25
iterations, about 8.5 minutes.

The fit diverges monotonically. The maximum `|g - g_target|` grows from 0.448 at
iteration 0 to 1.863 at iteration 24, so the initial guess is the best result. The
failure is localized: the fit places the first peak at r 0.33 nm with g 2.59,
where the target peak is at 0.28 nm with g 2.64, and it overshoots the second
shell. The rest of the curve matches within about 0.05.

Conclusion: a one-site isotropic pair potential cannot place water's sharp
tetrahedral first peak and the second shell at the same time. Real coarse-grain
water needs orientational or multi-site terms, and that erodes the particle-count
saving.

## Experiment 7 - Lyapunov horizon (`lyapunov.py`)

Two trajectories start from the same coordinates and velocities, with one water
molecule rigidly translated by 0.001 nm. Both run 5 ps at dt 2 fs. The RMS
per-atom separation grows from 3.4e-5 nm to 0.1 nm in 0.57 ps (285 steps at
dt 2 fs). The fitted growth rate is 8.1 /ps, a Lyapunov time of 0.12 ps (62
steps at dt 2 fs). A 0.004 nm perturbation gives 9.9 /ps, so the rate is not an
artifact of one size.

Consequence for a learned propagator: a step error of size `e` has a useful
horizon of `ln(0.1/e) / 8` ps, which is `ln(0.1/e) / 0.016` steps at dt 2 fs.
Even a perfect model with error 1e-6 nm holds for only 1.4 ps (720 steps). To
stay accurate over a 200 ps run, the error must be about `0.1 * exp(-8*200)`,
far below float64 precision. A pure learned propagator cannot serve as a long
trajectory generator. The sound uses are a learned step coupled to a correct
force, or an equilibrium sampler such as a Boltzmann generator, which does not
need time-accurate paths. The horizon bounds every propagator method, so it is
the first measurement any such method needs.

## Experiment 8 - MTS failure map (`mtsres.py`)

Question (idea 4): is the multiple-timestep instability a smooth force bias that a
learned correction can remove, or a resonance that no force correction fixes? Map
the conserved-energy drift against the timestep and the PME update factor on the
water box, 20 ps each.

| dt (fs) | factor 0 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| 2 | 0.0020 | 0.0027 | 0.0025 | 0.0028 | 0.0026 |
| 3 | 0.0025 | -0.0005 | 0.0033 | 0.0090 | **0.313** |
| 4 | 0.0036 | 0.0030 | 0.0108 | **2.99** | crash |
| 5 | 0.0026 | 0.0046 | **0.503** | crash | crash |

Drift is kJ/mol/ps/atom; the limit is 0.02. Every setting whose slow-force update
interval `factor * dt` is 12 fs or less passes. Every setting at 15 fs or more
fails or crashes. The boundary follows the product, not the factor alone.

Conclusion: the MTS failure is a stability limit on the PME update interval, near
15 fs. It is not a smooth force bias. A learned correction cannot extend the
interval past that limit without changing the splitting, because the instability
is not a force-magnitude error that a network can subtract. Idea 4 is closed on
this build.

## Experiment 9 - capacity test with neural networks (`torchfit.py`)

Experiment 4 reached R squared 0.814 with a linear model on the 384
pairwise-plus-angular features and failed the 0.95 gate. Experiment 9 separates
"too little capacity" from "too little information". It fits the same features
with neural networks of growing width.

| hidden width | force RMS | R squared |
|---|---|---|
| linear (exact) | 251.5 | 0.814 |
| 32 | 341.1 | 0.657 |
| 128 | 270.1 | 0.785 |
| 512 | 255.9 | 0.807 |

Force magnitude is 582.5 kJ/mol/nm. Extra capacity does not help: the widest
network reaches 0.807, below the exact linear fit at 0.814. The limit is the
information in the descriptor, not the model size. The descriptor is a function
of the instantaneous local positions, and the missing force is long-range PME
plus the velocity-dependent SETTLE constraint. Neither is present in the input.
Idea 1 stays closed even with a deep model on this hardware.

## Experiment 10 - learned conditional water proposal (`patchlearn.py`)

Goal: beat the random patch-move baseline of `patch.py` (fit `dU/kT = 1.872 n`)
with a conditional proposal over local water, scored by the exact GROMACS
potential and Metropolis.

Model: a ridge regression learns the water dipole direction from a local
descriptor (radial density, radial times `u.rhat`, and frame-directional sums)
gathered over atoms within 0.6 nm.

Result: the proposal loses on raw acceptance. Per-water cost is about 6.5 kT at
n=1, against the baseline 0.9 kT, and the joint slope is 11.6 kT/water against
1.87. An early small run showed 0.026 kT at n=1; that was seed noise. A local
frame built from three neighbour oxygens was unstable (R squared 0.13); the
dipole is only predictable in the fixed box frame (R squared 0.75, residual std
about 0.3).

Fair metric (`patchfair.py`). The raw comparison is unfair, because the learned
move is larger than the baseline move. At radius 0.4 nm (5 waters) the baseline
move reaches 0.017 to 0.055 nm RMS over rotation sigma 0.05 to 0.50 rad, while
the learned move reaches 0.096 nm. The learned `dU/kT` is 32 to 44 at that size.
Interpolating the baseline to 0.096 nm gives 110 to 170 kT, so the learned
proposal is 3 to 4 times more efficient per unit displacement. The model does
find lower-energy poses than a random move of the same size.

The limitation is sharp and not a noise effect. Scaling the sampling noise to
zero leaves the RMS at 0.096 nm, because the conditional mean itself is displaced
from the current pose by 20 to 60 degrees. The ridge conditional mean is not the
current configuration: given a frozen neighbourhood the mean sits away from the
correlated sample at hand. So the model cannot act as a local proposal at any
step size in this form.

Diagnosis: the conditional model is informative, but sampling from its mean makes
a large rotation where the baseline rotates 6 degrees. For n>1 the independent
per-water draws are not a joint sample, so the patch breaks its mutual
hydrogen-bond network.

What a winning proposal needs: a joint (correlated) proposal such as an
autoregressive or flow model trained to match the sampling context; a proposal
centred on the current state, not on the conditional mode; and an equivariant
model, since the box-frame trick works for only one fixed box.

## Experiment 11 - acceptance against displacement, three moves (`patchfrontier.py`)

Deliverable from Claude: one plot, acceptance against RMS displacement, three
curves, matched patch sizes, one exact scorer. The moves are a random jiggle
(`patch.py`), a rigid patch rotation (`phase2/surface_vs_volume.py`), and the
learned residual.

Fix first: `patchlearn.water_from` placed the hydrogens from the dipole alone,
which leaves the rotation about the dipole axis free. That added a spurious
rotation and explained part of Experiment 10's loss. The residual here rotates
the whole molecule by the minimal rotation that takes its current dipole to the
predicted one.

Setup: trained on the water-small reference trajectory, scored on
`phase2/npt_6.0.gro` (21087 atoms, 7029 waters) so the patch sizes match phase2.
Four radii give 19, 47, 195 and 1096 waters; 40 proposals per point.

Median frontier, `exp(-median dU/kT)`:

| waters | move | RMS nm | dU/kT | accept |
|---|---|---|---|---|
| 47 | rotate | 0.0063 | 2.44 | 8.7e-2 |
| 47 | learned | 0.0064 | -2.24 | 1.0 |
| 195 | rotate | 0.0064 | 6.99 | 9.3e-4 |
| 195 | learned | 0.0062 | -23.4 | 1.0 |
| 1096 | rotate | 0.0063 | 25.3 | 1.0e-11 |
| 1096 | learned | 0.0059 | -181.8 | 1.0 |

The learned curve sits above the rotation curve at matched RMS, and the gap
grows with patch size: about 4x at 19 waters, and over 70 orders of magnitude at
1096 waters. The random jiggle is below the rotation curve at every size. This is
the surface-scaling hope made concrete: the learned direction improves internal
structure, so its cost does not scale with patch volume.

Caveat, and it is important. The residual is a deterministic drift toward the
conditional mean, so it lowers the energy (`dU` negative) and is not a reversible
proposal. The raw `exp(-dU)` is not a valid Metropolis acceptance for it: the
reverse move is not the mirror, so a Hastings correction is required. The plot
therefore measures the information in the learned direction, not the acceptance
of a finished sampler. A valid stochastic form, and a detailed-balance check, are
the next step.

## Experiment 12 - valid Metropolis test of the learned drift (`patchbalance.py`)

Experiment 11's advantage was measured without the Hastings term, so it is not an
acceptance. This experiment puts the same learned direction into a valid
Metropolis-Hastings proposal and compares it with a symmetric random walk.

  rw     each patch water turns by a random rotation vector `N(0, sigma^2 I)`.
         The proposal is symmetric, so acceptance = `min(1, exp(-dU))`.
  drift  the rotation vector is `t*r + sigma*xi`, where `r` is the minimal
         rotation from the current dipole to the model conditional mean. The
         Hastings ratio `q(x|x')/q(x'|x)` is applied exactly.

Both moves are legal. Result, medians over 40 proposals, drift sigma 0.04:

| waters | move | RMS nm | dU/kT | log Hastings | true accept |
|---|---|---|---|---|---|
| 44 | rw | 0.0043 | 11.98 | 0.00 | 6.5e-6 |
| 44 | drift | 0.0044 | 11.35 | -4.05 | 4.7e-7 |
| 201 | rw | 0.0044 | 28.61 | 0.00 | 4.0e-13 |
| 201 | drift | 0.0046 | 19.04 | -78.22 | 7.1e-42 |

The learned drift does lower the energy a little more per unit displacement
(11.4 against 12.0 at 44 waters; 19-24 against 29 at 201 waters). But the
Hastings term is large and negative, and it grows with patch size, so the true
acceptance is far below the random walk at every matched size. As the drift step
grows past the noise, the penalty scales as `(t*r/sigma)^2` per water and kills
the proposal outright (log Hastings of -400 to -2200 in an earlier sweep).

Why: the conditional mean is a target, not a reversible drift. A proposal that
moves each water toward its conditional mean is biased downhill, and the Hastings
ratio charges back exactly that bias. A reversible proposal needs a score or
force estimate, not a mean displacement, and Experiments 3 and 4 already showed
that a position-only model cannot supply that force.

Conclusion: the Experiment 11 frontier measures the information in the learned
direction, and that information is real. It is not usable as a Metropolis
proposal in this form. Making it legal removes more than it adds.

## Verdict

- **Learned propagator (idea 1):** the compute is affordable, but the accuracy
  gate failed. The three-body test reached R squared 0.814 against a 0.95 gate,
  and oxygen-only forces reached only 0.849. The residual is long-range PME force
  and a velocity-dependent constraint force, neither of which is a function of
  the instantaneous local positions. A local position-only force model is closed.
  A model trained end-to-end to the next positions, not to forces, could still
  absorb these parts statistically, but that is a larger model off this CPU.
  The Lyapunov probe (experiment 7) then closes the route: the trajectory
  decorrelates in about 0.1 ps, so any per-step error larger than about 1e-70
  ruins a 200 ps path. A learned step cannot replace the true force over a long
  run. It can only sit inside a corrector, or sample equilibrium.
- **Learned MTS splitting (idea 4):** closed. The instability is a resonance at a
  slow-force update interval near 15 fs (experiment 8), not a smooth force bias,
  so a learned correction cannot remove it without changing the splitting.
- **Coarse-grain plus backmap (idea 3):** the only route with a predictable large
  gain, but three probes now show the cost. It needs about 16:1 for 10x
  (experiment 2). One-site water cannot hold the first peak and the second shell
  at once (experiment 6), so orientational or multi-site terms are required, and
  they cut the saving. The protein perturbs water only inside 0.6 nm (experiment
  5), so a distance-adaptive map is mostly uniform bulk and gains little. A
  coarse-grain water model with the required quality is a research project, not a
  small change to this tool.
- **Boltzmann generator (idea 2):** best for the equilibrium checks only. Torch is
  now installed in `ml/.venv311`, so this is no longer blocked by the stack, but
  training still wants a GPU, and the tool's drift check measures dynamics, which
  a generator does not satisfy.

## Next small experiment (if idea 3 is pursued)

The three probes close the cheap part of idea 3. A useful next step is larger: a
many-site or orientational coarse-grain water model, fitted to the O-O, O-H and
H-H distributions together, then a measured speed and RDF check against this
reference. Treat it as a research track with its own budget, not a small probe.

## Files

- `ceiling.py` - experiments 1 and 2, writes `ceiling.txt`.
- `forcefit.py` - experiment 3, caches `forcefit_data.npz`.
- `threebody.py` - experiment 4, writes `threebody.txt`.
- `cgprobe.py` - bulk-water order and spatial correlation, writes `cgprobe.txt`.
- `cgvillin.py` - water order against distance to the protein, writes
  `cgvillin.txt`.
- `cg_ibi.py` - one-site water IBI, writes `cg_ibi_rdf.txt`.
- `lyapunov.py` - trajectory divergence rate, writes `lyapunov.txt`.
- `patch.py` - exact-GROMACS water-patch scorer and random baseline, writes
  `patch/patch.txt`.
- `patchlearn.py` - learned conditional water proposal, writes `patchlearn.txt`.
- `patchfair.py` - fair metric for the two water proposals, writes
  `patchfair.txt`.
- `patchbalance.py` - valid Metropolis comparison of rw against learned drift,
  writes `patchbalance.txt` and `patchbalance.png`.
- `patchfrontier.py` - three-move acceptance frontier, writes `patchfrontier.txt`
  and `patchfrontier.png`; `patchfrontier_plot.py` aggregates and plots it.
- `mtsres.py` - MTS stability map, writes `mtsres.txt`.
- `torchfit.py` - neural-network capacity test, writes `torchfit.txt`. Runs under
  `ml/.venv311` (python3.11 plus torch 2.14.0).
- `forcedump.mdp` - mdp for the 1 ps coordinate plus force dump.
- `dump.trr`, `dump.txt` - rerun trajectory and its text dump (large, do not
  commit).

## Experiment 13 - three-body coarse-grain water (`cg3.py`)

Goal: a one-site coarse-grain (CG) water, one particle per molecule, would cut
the particle count 3x. The earlier one-site IBI (`cg_ibi.py`) diverged. `cg3.py`
adds a three-body angular term `u3(cos theta)` localised to the first
coordination shell by `w(r) = exp(-((r-0.28)/0.05)^2)`, and updates `u2(r)` and
`u3(cos)` jointly by IBI against the atomistic O-O RDF and the envelope-weighted
`P(cos theta)`. The three-body force is a stated derivation and was verified
against a finite-difference gradient (max error 5e-8 against a force of 286; net
force 8.5e-14).

Result (n 256, box 1.97 nm, 25 iterations): the angular distribution converges
(`max|P-Pt|` 0.019 -> 0.002), but the O-O RDF plateaus at `max|g-gt|` ~0.45
(rms 0.11) and then drifts up. Final comparison:

| r (nm) | model g | target g |
|---|---|---|
| 0.26 | 0.78 | 0.30 |
| 0.28 | 1.82 | 2.42 |
| 0.30 | 1.86 | 2.00 |
| 0.34 | 1.05 | 0.93 |
| 0.40 | 1.01 | 0.97 |

The model first peak sits at r=0.30 nm with g=1.82; the target sits at 0.28 nm
with g=2.42. The second shell and beyond match (g 0.93-1.02).

Conclusion: the 3-body term fixes the angular order and the outer structure, but
a one-site model still cannot reproduce the sharp first peak, which is 25 percent
too low and shifted by 0.02 nm. The one-site CG water route is closed.
Reproducing the first peak needs a multi-site model (for example two sites), and
that erodes the particle-count saving that was the point.

Files: `cg3.py` (three-body IBI, writes `cg3_state.npy` and `cg3_final_g.npy`);
`bg/` (Boltzmann generator prototype, see `bg/README.md`).
