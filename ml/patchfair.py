"""Fair comparison of the random and learned water-patch proposals (phase0/ml).

Milestone 2b. The earlier comparison was unfair: the learned proposal rotates a
water by about 20 degrees while the random baseline rotates about 6 degrees, so
the baseline wins on a smaller move. This script reports both methods against the
mean root-mean-square displacement of the patch, and compares them at matched
displacement.

The metric is the exact GROMACS potential of the full box (patch.energies), the
same as patch.py and patchlearn.py. Acceptance is reported as exp(-median dU/kT),
the phase2/RESULTS.md convention.

Writes patchfair.txt. Scratch in patchfair/.
"""

import argparse
import math
import os

import numpy as np

import patch
import patchlearn

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "patchfair")


def displacement(coords, new, waters, box):
    diff = patch.min_image(new - coords, box)
    idx = np.concatenate([np.array(w) for w in waters])
    return float(np.sqrt(np.mean(np.sum(diff[idx] ** 2, axis=1))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", default=os.path.join(PROJECT, "phase0", "systems", "npt_3.0.gro"))
    parser.add_argument("--top", default=os.path.join(PROJECT, "phase0", "systems", "topol_3.0.top"))
    parser.add_argument("--radius", type=float, default=0.4)
    parser.add_argument("--proposals", type=int, default=80)
    parser.add_argument("--sigmas", default="0.05,0.10,0.15,0.20,0.30,0.40,0.50")
    parser.add_argument("--alphas", default="0.10,0.20,0.35,0.50,1.00")
    parser.add_argument("--sigma-t", type=float, default=0.01)
    parser.add_argument("--max-frames", type=int, default=40)
    parser.add_argument("--seed", type=int, default=21)
    args = parser.parse_args()

    os.makedirs(WORK, exist_ok=True)
    patch.WORK = WORK
    patchlearn.WORK = WORK
    patch.prepare(args.gro, args.top)

    names, coords, box = patch.read_gro(args.gro)
    waters = patch.water_index(names)
    oxygen = [w[0] for w in waters]
    print(f"atoms {len(names)}  waters {len(waters)}  box {box}")

    frames = patchlearn.ref_frames(names, box)
    X, Y = patchlearn.training_set(frames, box, oxygen, waters, len(names), args.max_frames)
    model = patchlearn.fit(X, Y, 1e-3)
    print(f"training examples {len(X)}  dipole R2 {np.round(model['r2'], 3)}")

    rng = np.random.default_rng(args.seed)
    radius = args.radius
    inner = rng.uniform(radius + 0.3, min(box) - radius - 0.3, size=3)
    dist = np.linalg.norm(patch.min_image(coords[np.array(oxygen)] - inner, box), axis=1)
    patch_oxy = [oxygen[k] for k in np.where(dist <= radius)[0]]
    n0 = len(patch_oxy)
    print(f"patch waters {n0}")

    frames_out, labels, rms = [coords.copy()], [("base", 0.0)], [0.0]
    for sigma in [float(s) for s in args.sigmas.split(",")]:
        for _ in range(args.proposals):
            new, n = patch.make_proposal(coords, waters, patch_oxy, rng, args.sigma_t, sigma)
            frames_out.append(new)
            labels.append((f"random sigma={sigma:.2f}", sigma))
            rms.append(displacement(coords, new, [w for w in waters if any(o == w[0] for o in patch_oxy)], box))
    for alpha in [float(a) for a in args.alphas.split(",")]:
        for _ in range(args.proposals):
            new, n = patchlearn.sample_patch(coords, box, oxygen, waters, patch_oxy, model,
                                             rng, len(names), args.sigma_t, scale=alpha)
            frames_out.append(new)
            labels.append((f"learned alpha={alpha:.2f}", alpha))
            rms.append(displacement(coords, new, [w for w in waters if any(o == w[0] for o in patch_oxy)], box))

    path = os.path.join(WORK, "frames.gro")
    with open(path, "w") as handle:
        for i, frame in enumerate(frames_out):
            patch.write_gro_frame(handle, f"frame {i}", names, frame, box)

    pot = patch.energies(path)[:, 0]
    base = pot[0]
    print(f"base potential {base:.3f} kJ/mol  repeat spread {pot[0] - base:.4f}")

    lines = [f"system {os.path.basename(args.gro)}  patch waters {n0}  radius {radius}",
             f"proposals per method {args.proposals}  kT {patch.KT:.4f} kJ/mol",
             "",
             f"{'method':>20} {'rms_nm':>8} {'dU/kT med':>10} {'A_med':>10} {'A_mean':>10}"]
    idx = 1
    rows = []
    for sigma in [float(s) for s in args.sigmas.split(",")]:
        seg = slice(idx, idx + args.proposals)
        du = np.array([pot[idx + j] - base for j in range(args.proposals)])
        r = float(np.mean(rms[idx:idx + args.proposals]))
        med = float(np.median(du / patch.KT))
        acc = float(np.mean(np.minimum(1.0, np.exp(-du / patch.KT))))
        rows.append((f"random {sigma:.2f}", r, med, acc))
        lines.append(f"{f'random sigma={sigma:.2f}':>20} {r:8.4f} {med:10.3f} "
                     f"{math.exp(-max(0.0, med)):10.3e} {acc:10.3e}")
        idx += args.proposals
    learned = []
    for alpha in [float(a) for a in args.alphas.split(",")]:
        du = np.array([pot[idx + j] - base for j in range(args.proposals)])
        r = float(np.mean(rms[idx:idx + args.proposals]))
        med = float(np.median(du / patch.KT))
        acc = float(np.mean(np.minimum(1.0, np.exp(-du / patch.KT))))
        learned.append((r, med, acc))
        lines.append(f"{f'learned alpha={alpha:.2f}':>20} {r:8.4f} {med:10.3f} "
                     f"{math.exp(-max(0.0, med)):10.3e} {acc:10.3e}")
        idx += args.proposals

    rr = np.array([row[1] for row in rows])
    mm = np.array([row[2] for row in rows])
    lines += ["", "matched-displacement comparison (random dU/kT interpolated at the learned rms):"]
    for r, med, acc in learned:
        if rr.min() <= r <= rr.max():
            pred = float(np.interp(r, rr, mm))
            lines.append(f"  learned rms {r:.4f} nm  dU/kT {med:.3f}  random {pred:.3f}  "
                         f"ratio {med / max(pred, 1e-9):.2f}")
        else:
            lines.append(f"  learned rms {r:.4f} nm outside random range "
                         f"[{rr.min():.4f}, {rr.max():.4f}]")

    out = os.path.join(HERE, "patchfair.txt")
    with open(out, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
