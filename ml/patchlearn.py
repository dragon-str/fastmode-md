"""Learned local water proposal (phase0/ml).

Milestone 2: replace the random patch move of patch.py with a conditional
proposal. A ridge regression learns p(dipole | local context) from the reference
trajectory. A proposal samples each patch water from that conditional, in random
order, so already sampled patch waters enter the context of the later ones.

OUTCOME: the proposal loses to the random baseline. Single-water cost is about
6.5 kT against 0.9 kT, and the joint slope is 11.6 kT/water against 1.87. The
conditional model is informative (dipole R2 0.75 in the box frame), but sampling
from it makes a larger rotation than the baseline, and independent draws are not
a joint sample, so a multi-water patch breaks its hydrogen-bond network. See
patchlearn.txt for the diagnosis. A winning proposal needs a joint model.

Writes patchlearn.txt. Scratch in patchlearn/.
"""

import argparse
import math
import os

import numpy as np

import patch

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "patchlearn")
RC = 0.6
N_RBF = 16
RBF_W = 0.05
OH = 0.09572
HALF_ANGLE = math.radians(104.52 / 2.0)


def rbf_centers():
    return np.linspace(0.0, RC, N_RBF)


CENTERS = rbf_centers()


def frame_from(p1, p2, p3):
    e1 = p2 - p1
    n1 = np.linalg.norm(e1)
    if n1 < 1e-9:
        return np.eye(3)
    e1 = e1 / n1
    e2 = p3 - p1
    e2 = e2 - np.dot(e2, e1) * e1
    n2 = np.linalg.norm(e2)
    if n2 < 1e-9:
        e2 = np.array([-e1[1], e1[0], 0.0])
        n2 = np.linalg.norm(e2)
        if n2 < 1e-9:
            e2 = np.array([0.0, -e1[2], e1[1]])
            n2 = np.linalg.norm(e2)
    e2 = e2 / n2
    e3 = np.cross(e1, e2)
    return np.stack([e1, e2, e3])


def context_of(coords, oxygen, i, box, w):
    oi = coords[oxygen[i]]
    delta = patch.min_image(coords - oi, box)
    dist = np.linalg.norm(delta, axis=1)
    sel = np.where(dist <= RC)[0]
    sel = sel[~np.isin(sel, w)]
    return sel, delta, dist


def choose_frame(coords, sel, dist, U):
    order = sel[np.argsort(dist[sel])]
    anchors = [a for a in order if np.linalg.norm(U[a]) > 0.5]
    if len(anchors) < 3:
        anchors = list(order)
    p1, p2, p3 = coords[anchors[0]], coords[anchors[1]], coords[anchors[2]]
    return frame_from(p1, p2, p3)


def features(sel, delta, dist, U):
    if len(sel) == 0:
        return np.zeros(5 * N_RBF)
    r = dist[sel]
    g = np.exp(-((r[:, None] - CENTERS[None, :]) / RBF_W) ** 2)
    rhat = delta[sel] / np.maximum(r[:, None], 1e-9)
    cu = np.einsum("mc,mc->m", U[sel], rhat)
    directional = np.einsum("mk,mc->kc", g, rhat).ravel()
    polarized = np.einsum("mk,mc->kc", g, U[sel]).ravel()
    return np.concatenate([g.sum(0), (g * cu[:, None]).sum(0),
                           (g * (cu ** 2)[:, None]).sum(0), directional, polarized])


def build_dipoles(coords, waters, natoms):
    U = np.zeros((natoms, 3))
    for w in waters:
        o, h1, h2 = coords[w[0]], coords[w[1]], coords[w[2]]
        u = (h1 - o) + (h2 - o)
        n = np.linalg.norm(u)
        if n > 1e-9:
            U[w[0]] = u / n
    return U


def training_set(frames, box, oxygen, waters, natoms, max_frames):
    X, Y = [], []
    for coords in frames[:max_frames]:
        U = build_dipoles(coords, waters, natoms)
        for i in range(len(oxygen)):
            w = waters[i]
            sel, delta, dist = context_of(coords, oxygen, i, box, w)
            if len(sel) < 3:
                continue
            feat = features(sel, delta, dist, U)
            X.append(feat)
            Y.append(U[oxygen[i]])
    return np.array(X), np.array(Y)


def water_from(o, u):
    axis = np.cross([0.0, 0.0, 1.0], u)
    s = np.linalg.norm(axis)
    if s < 1e-9:
        if u[2] >= 0:
            rot = np.eye(3)
        else:
            rot = np.diag([1.0, -1.0, -1.0])
    else:
        k = axis / s
        theta = math.acos(max(-1.0, min(1.0, float(u[2]))))
        cross = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        rot = np.eye(3) + math.sin(theta) * cross + (1 - math.cos(theta)) * (cross @ cross)
    canon = np.array([[math.sin(HALF_ANGLE), 0.0, math.cos(HALF_ANGLE)],
                      [-math.sin(HALF_ANGLE), 0.0, math.cos(HALF_ANGLE)]]) * OH
    return o + canon @ rot.T


def fit(X, Y, lam):
    xm, xs = X.mean(0), X.std(0) + 1e-9
    Xs = (X - xm) / xs
    A = Xs.T @ Xs + lam * np.eye(Xs.shape[1])
    w = np.linalg.solve(A, Xs.T @ Y)
    pred = Xs @ w
    resid = Y - pred
    return {"xm": xm, "xs": xs, "w": w, "std": resid.std(0), "r2": 1 - resid.var(0) / Y.var(0)}


def sample_patch(coords, box, oxygen, waters, patch_oxy, model, rng, natoms, sigma_t, scale=1.0):
    out = coords.copy()
    U = build_dipoles(out, waters, natoms)
    patch = [w for w in waters if any(o == w[0] for o in patch_oxy)]
    rng.shuffle(patch)
    for w in patch:
        i = oxygen.index(w[0])
        sel, delta, dist = context_of(out, oxygen, i, box, w)
        if len(sel) < 3:
            continue
        feat = features(sel, delta, dist, U)
        z = (feat - model["xm"]) / model["xs"]
        mu = z @ model["w"]
        u = mu + scale * rng.normal(scale=model["std"])
        n = np.linalg.norm(u)
        u = u / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])
        o_new = out[w[0]] + rng.normal(scale=sigma_t, size=3)
        hw = water_from(o_new, u)
        out[w[0]] = o_new
        out[w[1]] = hw[0]
        out[w[2]] = hw[1]
        U[w[0]] = u
    return np.mod(out, box), len(patch)


def ref_frames(names, box):
    tgro = os.path.join(WORK, "ref_frames.gro")
    xtc = os.path.join(PROJECT, "phase0", "out", "water_small", "runs", "reference",
                       "dt2fs_hmroff_mtsoff_nstlist10_tol0.005", "run.xtc")
    if not os.path.exists(tgro):
        res = patch.run([patch.GMX, "trjconv", "-f", xtc, "-s", "rerun.tpr",
                         "-skip", "40", "-pbc", "mol", "-o", "ref_frames.gro"],
                        WORK, stdin="0\n")
        if res.returncode != 0:
            raise SystemExit(res.stderr[-2000:])
    with open(tgro) as handle:
        lines = handle.readlines()
    frames, i = [], 0
    while i < len(lines):
        natoms = int(lines[i + 1])
        chunk = lines[i + 2:i + 2 + natoms]
        coords = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in chunk])
        frames.append(coords)
        i += 2 + natoms + 1
    return frames


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", default=os.path.join(PROJECT, "phase0", "systems", "npt_3.0.gro"))
    parser.add_argument("--top", default=os.path.join(PROJECT, "phase0", "systems", "topol_3.0.top"))
    parser.add_argument("--max-frames", type=int, default=40)
    parser.add_argument("--lam", type=float, default=1e-3)
    parser.add_argument("--radii", default="0.3,0.4,0.5,0.6")
    parser.add_argument("--proposals", type=int, default=40)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--sigma-t", type=float, default=0.005)
    args = parser.parse_args()

    os.makedirs(WORK, exist_ok=True)
    patch.WORK = WORK
    patch.prepare(args.gro, args.top)

    names, coords, box = patch.read_gro(args.gro)
    waters = patch.water_index(names)
    oxygen = [w[0] for w in waters]
    print(f"atoms {len(names)}  waters {len(waters)}")

    frames = ref_frames(names, box)
    print(f"reference frames {len(frames)}")
    X, Y = training_set(frames, box, oxygen, waters, len(names), args.max_frames)
    print(f"training examples {len(X)}  features {X.shape[1]}")
    model = fit(X, Y, args.lam)
    print(f"train R2 per output: {np.round(model['r2'], 3)}")
    print(f"residual std: {np.round(model['std'], 4)}")

    rng = np.random.default_rng(args.seed)
    radii = [float(r) for r in args.radii.split(",")]
    frames_out, labels = [coords.copy()], [(0.0, 0.0)]
    for radius in radii:
        inner = rng.uniform(radius + 0.3, min(box) - radius - 0.3, size=3)
        dist = np.linalg.norm(patch.min_image(coords[np.array(oxygen)] - inner, box), axis=1)
        patch_oxy = [oxygen[k] for k in np.where(dist <= radius)[0]]
        for _ in range(args.proposals):
            new, n = sample_patch(coords, box, oxygen, waters, patch_oxy, model, rng,
                                  len(names), args.sigma_t)
            frames_out.append(new)
            labels.append((radius, n))

    path = os.path.join(WORK, "frames.gro")
    with open(path, "w") as handle:
        for i, frame in enumerate(frames_out):
            patch.write_gro_frame(handle, f"frame {i}", names, frame, box)

    pot = patch.energies(path)[:, 0]
    base = pot[0]
    print(f"base potential {base:.3f} kJ/mol  repeat spread {pot[0] - base:.4f}")

    lines = [f"system {os.path.basename(args.gro)}  waters {len(waters)}",
             f"training frames {args.max_frames}  examples {len(X)}  lambda {args.lam}",
             f"train R2 {np.round(model['r2'], 3).tolist()}",
             "",
             f"{'radius':>7} {'n_waters':>8} {'dU/kT med':>10} {'A_med':>10} {'A_mean':>10}"]
    ns, dus = [], []
    idx = 1
    for radius in radii:
        du = np.array([pot[idx + j] - base for j in range(args.proposals)])
        idx += args.proposals
        n = float(np.mean([labels[idx - args.proposals + j][1] for j in range(args.proposals)]))
        med = float(np.median(du / patch.KT))
        acc = float(np.mean(np.minimum(1.0, np.exp(-du / patch.KT))))
        ns.append(n)
        dus.append(med)
        lines.append(f"{radius:7.2f} {n:8.1f} {med:10.3f} {math.exp(-max(0.0, med)):10.3e} {acc:10.3e}")
    slope, intercept = np.polyfit(ns, dus, 1)
    lines += ["", f"fit dU/kT = {slope:.3f} * n {intercept:+.3f}",
              "baseline patch.py fit dU/kT = 1.872 * n - 0.276"]
    out = os.path.join(HERE, "patchlearn.txt")
    with open(out, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
