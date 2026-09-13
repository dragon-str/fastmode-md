"""Valid Metropolis test of the learned water direction (phase0/ml).

The residual move in patchfrontier.py drifts each water toward the model
conditional mean and keeps the move only when the energy drops.  That move is not
reversible, so its raw exp(-dU) is not a Metropolis acceptance.  This script puts
the same learned direction into a proper Metropolis-Hastings proposal and compares
it against a symmetric random walk at matched RMS displacement.

Both moves are valid, so the comparison is fair:

  rw      each patch water turns by a random rotation vector N(0, sigma^2 I).
          The proposal is symmetric, so acceptance = min(1, exp(-dU)).
  drift   the rotation vector is t * r + sigma * xi, where r is the minimal
          rotation that takes the current dipole to the model conditional mean.
          The proposal is not symmetric, so the Hastings ratio q(x|x')/q(x'|x)
          is applied.

Scores forward frames with the exact GROMACS potential (patch.energies).  Writes
patchbalance.txt and patchbalance.png.  Scratch in patchbalance/.
"""

import argparse
import math
import os

import numpy as np

import patch
import patchlearn as pl
from patchfrontier import rodrigues, rms_disp

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
TRAIN = os.path.join(HERE, "patchbalance", "train")
SCORE = os.path.join(HERE, "patchbalance", "score")

RW_SIGMAS = [0.02, 0.04, 0.06, 0.09, 0.12, 0.16]
DRIFT_T = [0.02, 0.04, 0.07, 0.10, 0.14, 0.18]
DRIFT_SIGMA = 0.04


def rotvec(u0, u):
    axis = np.cross(u0, u)
    s = float(np.linalg.norm(axis))
    c = float(u0 @ u)
    if s < 1e-12:
        if c >= 0.0:
            return np.zeros(3)
        perp = np.array([1.0, 0.0, 0.0]) if abs(u0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        k = np.cross(u0, perp)
        k /= np.linalg.norm(k)
        return k * math.pi
    return axis / s * math.atan2(s, c)


def predict(coords, box, oxygen, U, w):
    i = oxygen.index(w[0])
    sel, delta, dist = pl.context_of(coords, oxygen, i, box, w)
    if len(sel) < 3:
        return None
    feat = pl.features(sel, delta, dist, U)
    mu = (feat - model_ref["xm"]) / model_ref["xs"] @ model_ref["w"]
    n = np.linalg.norm(mu)
    if n < 1e-9:
        return None
    return mu / n


def targets(coords, box, oxygen, waters, U):
    out = {}
    for w in waters:
        mu = predict(coords, box, oxygen, U, w)
        if mu is not None:
            out[w[0]] = rotvec(U[w[0]], mu)
    return out


def apply(coords, w, rho):
    ang = float(np.linalg.norm(rho))
    if ang < 1e-12:
        return
    k = rho / ang
    o = coords[w[0]]
    coords[w[1]] = o + rodrigues(coords[w[1]] - o, k, ang)
    coords[w[2]] = o + rodrigues(coords[w[2]] - o, k, ang)


def proposal(base, box, oxygen, waters, U, rng, t, sigma):
    out = base.copy()
    r = targets(base, box, oxygen, waters, U)
    rr = {}
    logq_fwd = 0.0
    logq_rev = 0.0
    applied = {}
    for w in waters:
        if w[0] not in r:
            continue
        xi = rng.normal(size=3)
        rho = t * r[w[0]] + sigma * xi
        applied[w[0]] = rho
        apply(out, w, rho)
        logq_fwd += -0.5 * float(xi @ xi) - 1.5 * math.log(2 * math.pi * sigma * sigma)
    U1 = pl.build_dipoles(out, waters, len(base))
    r2 = targets(out, box, oxygen, waters, U1)
    for w in waters:
        if w[0] not in r2 or w[0] not in applied:
            continue
        rho = applied[w[0]]
        resid = rho / sigma + (t / sigma) * r2[w[0]]
        logq_rev += -0.5 * float(resid @ resid) - 1.5 * math.log(2 * math.pi * sigma * sigma)
    return np.mod(out, box), logq_fwd, logq_rev


def emit_plot(rows, radii, n_by_radius, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    moves = [("rw", "o", "tab:gray"), ("drift", "^", "tab:red")]
    fig, axes = plt.subplots(1, len(radii), figsize=(3.6 * len(radii), 4.0),
                             sharey=True)
    if len(radii) == 1:
        axes = [axes]
    for ax, radius in zip(axes, radii):
        for move, marker, colour in moves:
            pts = sorted([(r["rms"], r["accept"]) for r in rows
                          if r["move"] == move and r["radius"] == radius])
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [max(p[1], 1e-12) for p in pts]
            ax.plot(xs, ys, marker=marker, color=colour, label=move, lw=1.5)
        ax.set_yscale("log")
        ax.set_xlabel("RMS displacement (nm)")
        ax.set_title(f"{n_by_radius[radius]:.0f} waters")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("true MH acceptance")
    axes[0].legend()
    fig.suptitle("Valid Metropolis: random walk vs learned drift (exact scorer)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", default=os.path.join(PROJECT, "phase2", "npt_6.0.gro"))
    parser.add_argument("--top", default=os.path.join(PROJECT, "phase2", "topol_6.0.top"))
    parser.add_argument("--train-gro", default=os.path.join(PROJECT, "phase0", "systems", "npt_3.0.gro"))
    parser.add_argument("--train-top", default=os.path.join(PROJECT, "phase0", "systems", "topol_3.0.top"))
    parser.add_argument("--radii", default="0.7,1.14,2.0")
    parser.add_argument("--proposals", type=int, default=40)
    parser.add_argument("--max-frames", type=int, default=40)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    os.makedirs(TRAIN, exist_ok=True)
    os.makedirs(SCORE, exist_ok=True)

    patch.WORK = TRAIN
    patch.prepare(args.train_gro, args.train_top)
    pl.WORK = TRAIN
    tnames, tcoords, tbox = patch.read_gro(args.train_gro)
    twaters = patch.water_index(tnames)
    toxygen = [w[0] for w in twaters]
    frames = pl.ref_frames(tnames, tbox)
    X, Y = pl.training_set(frames, tbox, toxygen, twaters, len(tnames), args.max_frames)
    global model_ref
    model_ref = pl.fit(X, Y, 1e-3)
    print(f"train examples {len(X)}  R2 {np.round(model_ref['r2'], 3)}")

    patch.WORK = SCORE
    patch.prepare(args.gro, args.top)
    names, coords, box = patch.read_gro(args.gro)
    waters = patch.water_index(names)
    oxygen = [w[0] for w in waters]
    natoms = len(names)
    print(f"score system {os.path.basename(args.gro)}  atoms {natoms}  waters {len(waters)}")

    rng = np.random.default_rng(args.seed)
    radii = [float(r) for r in args.radii.split(",")]
    if args.smoke:
        radii = radii[:1]
    props = 3 if args.smoke else args.proposals
    rw_sigmas = RW_SIGMAS[:2] if args.smoke else RW_SIGMAS
    drift_ts = DRIFT_T[:2] if args.smoke else DRIFT_T

    O = np.array(oxygen)
    selected = {}
    for radius in radii:
        centre = rng.uniform(radius + 0.05, min(box) - radius - 0.05, size=3)
        dist = np.linalg.norm(patch.min_image(coords[O] - centre, box), axis=1)
        keep = np.where(dist <= radius)[0]
        selected[radius] = (centre, [waters[k] for k in keep])
        print(f"radius {radius:.2f}  waters {len(keep)}")

    frames_out = [coords.copy()]
    rows = []

    for radius in radii:
        _, pw = selected[radius]
        patch_atoms = [w[0] for w in pw]
        for sigma in rw_sigmas:
            for _ in range(props):
                U = pl.build_dipoles(coords, waters, natoms)
                new, lf, lr = proposal(coords, box, oxygen, pw, U, rng, 0.0, sigma)
                rows.append({"move": "rw", "radius": radius, "step": sigma,
                             "rms": rms_disp(coords, new, patch_atoms, box),
                             "log_hast": lr - lf})
                frames_out.append(new)
        for t in drift_ts:
            for _ in range(props):
                U = pl.build_dipoles(coords, waters, natoms)
                new, lf, lr = proposal(coords, box, oxygen, pw, U, rng, t, DRIFT_SIGMA)
                rows.append({"move": "drift", "radius": radius, "step": t,
                             "rms": rms_disp(coords, new, patch_atoms, box),
                             "log_hast": lr - lf})
                frames_out.append(new)

    path = os.path.join(SCORE, "frames.gro")
    with open(path, "w") as handle:
        for i, frame in enumerate(frames_out):
            patch.write_gro_frame(handle, f"frame {i}", names, frame, box)
    print(f"frames {len(frames_out)}")

    pot = patch.energies(path)[:, 0]
    base = pot[0]
    print(f"base potential {base:.3f} kJ/mol  spread {pot[0] - base:.4f}")

    for r, p in zip(rows, pot[1:]):
        r["du"] = float((p - base) / patch.KT)
        log_alpha = -max(0.0, r["du"]) + r["log_hast"]
        r["accept"] = math.exp(min(0.0, log_alpha))

    n_by_radius = {radius: len(selected[radius][1]) for radius in radii}
    lines = [f"score system {os.path.basename(args.gro)}  atoms {natoms}",
             f"train system {os.path.basename(args.train_gro)}  R2 {np.round(model_ref['r2'], 3).tolist()}",
             f"drift sigma {DRIFT_SIGMA}  proposals per point {props}  seed {args.seed}"
             f"  kT {patch.KT:.4f}",
             ""]
    for radius in radii:
        lines.append(f"radius {radius:.2f}  waters {n_by_radius[radius]}")
        lines.append(f"{'move':>6} {'step':>7} {'RMS nm':>9} {'dU/kT':>9} {'logH':>8} {'accept':>10}")
        for move in ("rw", "drift"):
            pts = [r for r in rows if r["move"] == move and r["radius"] == radius]
            pts.sort(key=lambda r: (r["step"], r["rms"]))
            for r in pts:
                lines.append(f"{move:>6} {r['step']:7.3f} {r['rms']:9.4f}"
                             f" {r['du']:9.3f} {r['log_hast']:8.3f} {r['accept']:10.3e}")
            lines.append("")

    out = os.path.join(HERE, "patchbalance.txt")
    with open(out, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")

    emit_plot(rows, radii, n_by_radius, os.path.join(HERE, "patchbalance.png"))


if __name__ == "__main__":
    main()
