"""Acceptance versus displacement for three local water moves (phase0/ml).

Answers Claude's deliverable: one plot, acceptance against RMS displacement, for
three moves scored by the same exact GROMACS potential at matched patch sizes.

  jiggle    every water gets its own random translation and rotation
  rotate    the whole patch turns rigidly about its own centre
  learned   each water's dipole is nudged toward the conditional mean by t

The learned move is a residual form of patchlearn.py: it interpolates between the
current dipole and the predicted one, so t sets the step size and traces a
frontier. Training is on the water-small reference trajectory; scoring is on the
large box so the patch sizes reach about 1090 waters, matching phase2.

Writes patchfrontier.txt and patchfrontier.png. Scratch in patchfrontier/.
"""

import argparse
import math
import os

import numpy as np

import patch
import patchlearn as pl

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
TRAIN = os.path.join(HERE, "patchfrontier", "train")
SCORE = os.path.join(HERE, "patchfrontier", "score")

JIGGLE_SCALES = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
ROTATE_SIGMAS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12]
RESIDUAL_T = [0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0]


def rodrigues(rel, axis, angle):
    return (rel * math.cos(angle) + np.cross(axis, rel) * math.sin(angle)
            + np.outer(rel @ axis, axis) * (1 - math.cos(angle)))


def jiggle(coords, patch_atoms, rng, sigma, rot_sigma):
    out = coords.copy()
    for b in patch_atoms:
        grp = out[b:b + 3]
        com = grp.mean(axis=0)
        ax = rng.normal(size=3)
        ax /= np.linalg.norm(ax)
        out[b:b + 3] = (rodrigues(grp - com, ax, rng.normal(scale=rot_sigma))
                        + com + rng.normal(scale=sigma, size=3))
    return out


def rotate(coords, box, patch_atoms, centre, rng, sigma, radius):
    out = coords.copy()
    ax = rng.normal(size=3)
    ax /= np.linalg.norm(ax)
    angle = sigma / max(radius, 1e-3)
    idx = np.concatenate([np.arange(b, b + 3) for b in patch_atoms])
    rel = patch.min_image(out[idx] - centre, box)
    out[idx] = rodrigues(rel, ax, angle) + centre
    return out


def residual(coords, box, waters, oxygen, patch_waters, U, model, rng, t):
    out = coords.copy()
    order = list(patch_waters)
    rng.shuffle(order)
    for w in order:
        i = oxygen.index(w[0])
        sel, delta, dist = pl.context_of(out, oxygen, i, box, w)
        if len(sel) < 3:
            continue
        feat = pl.features(sel, delta, dist, U)
        mu = (feat - model["xm"]) / model["xs"] @ model["w"]
        n = np.linalg.norm(mu)
        if n < 1e-9:
            continue
        mu = mu / n
        u0 = U[w[0]]
        cand = u0 + t * (mu - u0)
        nc = np.linalg.norm(cand)
        u = cand / nc if nc > 1e-9 else u0
        axis = np.cross(u0, u)
        s = float(np.linalg.norm(axis))
        c = float(u0 @ u)
        if s < 1e-12:
            if c >= 0.0:
                continue
            perp = np.array([1.0, 0.0, 0.0]) if abs(u0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            k = np.cross(u0, perp)
            k /= np.linalg.norm(k)
            theta = math.pi
        else:
            k = axis / s
            theta = math.atan2(s, c)
        o = out[w[0]]
        out[w[1]] = o + rodrigues(out[w[1]] - o, k, theta)
        out[w[2]] = o + rodrigues(out[w[2]] - o, k, theta)
        U[w[0]] = u
    return np.mod(out, box)


def rms_disp(base, new, patch_atoms, box):
    idx = np.concatenate([np.arange(b, b + 3) for b in patch_atoms])
    d = patch.min_image(new[idx] - base[idx], box)
    return float(np.sqrt((d ** 2).sum(axis=1).mean()))


def make_plot(rows, radii, n_by_radius, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    moves = [("jiggle", "o", "tab:gray"), ("rotate", "s", "tab:blue"),
             ("learned", "^", "tab:red")]
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
    axes[0].set_ylabel("acceptance  exp(-median dU/kT)")
    axes[0].legend()
    fig.suptitle("Local water move efficiency: exact GROMACS scorer, matched patches")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", default=os.path.join(PROJECT, "phase2", "npt_6.0.gro"))
    parser.add_argument("--top", default=os.path.join(PROJECT, "phase2", "topol_6.0.top"))
    parser.add_argument("--train-gro", default=os.path.join(PROJECT, "phase0", "systems", "npt_3.0.gro"))
    parser.add_argument("--train-top", default=os.path.join(PROJECT, "phase0", "systems", "topol_3.0.top"))
    parser.add_argument("--radii", default="0.5,0.7,1.14,2.0")
    parser.add_argument("--proposals", type=int, default=50)
    parser.add_argument("--max-frames", type=int, default=40)
    parser.add_argument("--seed", type=int, default=21)
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
    model = pl.fit(X, Y, 1e-3)
    print(f"train examples {len(X)}  R2 {np.round(model['r2'], 3)}")

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
    props = 4 if args.smoke else args.proposals
    jiggle_scales = JIGGLE_SCALES[:2] if args.smoke else JIGGLE_SCALES
    rotate_sigmas = ROTATE_SIGMAS[:2] if args.smoke else ROTATE_SIGMAS
    residual_ts = RESIDUAL_T[:2] if args.smoke else RESIDUAL_T

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

    def emit(move, radius, step, new, patch_waters):
        patch_atoms = [w[0] for w in patch_waters]
        rows.append({"move": move, "radius": radius, "step": step,
                     "rms": rms_disp(coords, new, patch_atoms, box)})
        frames_out.append(new)

    for radius in radii:
        centre, pw = selected[radius]
        patch_atoms = [w[0] for w in pw]
        for s in jiggle_scales:
            for _ in range(props):
                emit("jiggle", radius, s,
                     jiggle(coords, patch_atoms, rng, 0.01 * s, 0.10 * s), pw)
        for sig in rotate_sigmas:
            for _ in range(props):
                emit("rotate", radius, sig,
                     rotate(coords, box, patch_atoms, centre, rng, sig, radius), pw)
        for t in residual_ts:
            for _ in range(props):
                U = pl.build_dipoles(coords, waters, natoms)
                emit("learned", radius, t,
                     residual(coords, box, waters, oxygen, pw, U, model, rng, t), pw)

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
        r["accept"] = math.exp(-max(0.0, r["du"]))

    n_by_radius = {radius: len(selected[radius][1]) for radius in radii}
    lines = [f"score system {os.path.basename(args.gro)}  atoms {natoms}",
             f"train system {os.path.basename(args.train_gro)}  R2 {np.round(model['r2'], 3).tolist()}",
             f"proposals per point {props}  seed {args.seed}  kT {patch.KT:.4f}",
             ""]
    for radius in radii:
        lines.append(f"radius {radius:.2f}  waters {n_by_radius[radius]}")
        lines.append(f"{'move':>8} {'step':>7} {'RMS nm':>9} {'dU/kT':>9} {'accept':>10}")
        for move in ("jiggle", "rotate", "learned"):
            pts = [r for r in rows if r["move"] == move and r["radius"] == radius]
            pts.sort(key=lambda r: r["rms"])
            step_key = {"jiggle": "scale", "rotate": "sigma", "learned": "t"}[move]
            for r in pts:
                lines.append(f"{move:>8} {r['step']:7.3f} {r['rms']:9.4f}"
                             f" {r['du']:9.3f} {r['accept']:10.3e}")
            lines.append("")
            _ = step_key

    out = os.path.join(HERE, "patchfrontier.txt")
    with open(out, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")

    make_plot(rows, radii, n_by_radius, os.path.join(HERE, "patchfrontier.png"))


if __name__ == "__main__":
    main()
