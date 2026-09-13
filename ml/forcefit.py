#!/usr/bin/env python3
"""Force-fit probe.

Can a cheap, rotationally equivariant model reproduce the GROMACS forces on the
water box? If yes, a learned force term is plausible. If no, a learned
propagator starts from a bad force field.

The model is a sum over neighbours of a learned radial coefficient times the
unit vector:
    F_i = sum_j c(r_ij) rhat_ij
The coefficient c(r) is a radial basis function (RBF) expansion. This form is
equivariant by construction, so the fit needs no rotation augmentation.

Data comes from `gmx dump -f dump.trr`, which writes coordinates and forces as
text. One frame in `--stride` frames is used, so the samples decorrelate.

Usage:
  python3 forcefit.py [--frames 40] [--stride 20] [--cutoff 0.5] [--rbf 16]
"""
import argparse
import os
import re
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
GMX = os.path.join(PROJECT, "build", "bin", "gmx")
TRR = os.path.join(HERE, "dump.trr")
DUMP = os.path.join(HERE, "dump.txt")
CACHE = os.path.join(HERE, "forcefit_data.npz")


def parse_dump(path, want_frames, stride):
    natoms = None
    frame = -1
    keep = False
    xs = []
    fs = []
    cur_x = None
    cur_f = None
    section = None
    arr = None
    fh = open(path)
    for line in fh:
        s = line.strip()
        if s.startswith("dump.trr frame"):
            if keep and cur_x is not None and cur_f is not None:
                xs.append(cur_x)
                fs.append(cur_f)
                if len(xs) >= want_frames:
                    break
            frame = int(re.search(r"frame (\d+):", s).group(1))
            keep = (frame % stride == 0)
            cur_x = cur_f = None
            section = None
            arr = None
            continue
        if s.startswith("natoms="):
            natoms = int(re.search(r"\d+", s.split("=", 1)[1]).group(0))
            continue
        if re.match(r"[xf] \(", s):
            section = s[0]
            arr = []
            continue
        if section is not None and s.startswith(section + "["):
            vals = s.split("{", 1)[1].split("}", 1)[0].split(",")
            arr.append([float(v) for v in vals])
            if len(arr) == natoms:
                a = np.array(arr)
                if section == "x":
                    cur_x = a
                else:
                    cur_f = a
                section = None
                arr = None
            continue
    fh.close()
    if keep and cur_x is not None and cur_f is not None and len(xs) < want_frames:
        xs.append(cur_x)
        fs.append(cur_f)
    return np.array(xs), np.array(fs)


def build_neighbours(x, box):
    d = x[:, None, :] - x[None, :, :]
    d -= box * np.round(d / box)
    r = np.sqrt((d ** 2).sum(-1))
    np.fill_diagonal(r, np.inf)
    return d, r


def design(x_frames, box, cutoff, centers, width, natoms):
    nf = len(centers) * 3
    X = np.zeros((len(x_frames) * natoms, nf))
    for fi, x in enumerate(x_frames):
        d, r = build_neighbours(x, box)
        mask = r < cutoff
        rr = r[mask]
        u = d[mask] / rr[:, None]
        src = np.repeat(np.arange(natoms), mask.sum(1))
        phi = np.exp(-((rr[:, None] - centers[None, :]) / width) ** 2)
        off = fi * natoms
        for k in range(len(centers)):
            contrib = np.zeros((natoms, 3))
            np.add.at(contrib, src, phi[:, k][:, None] * u)
            X[off:off + natoms, k * 3:k * 3 + 3] = contrib
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=40)
    ap.add_argument("--stride", type=int, default=20)
    ap.add_argument("--cutoff", type=float, default=0.5)
    ap.add_argument("--rbf", type=int, default=16)
    args = ap.parse_args()

    if os.path.exists(CACHE):
        z = np.load(CACHE)
        xs, fs = z["x"], z["f"]
    else:
        if not os.path.exists(DUMP):
            with open(DUMP, "w") as out:
                subprocess.run([GMX, "dump", "-f", TRR], stdout=out,
                               stderr=subprocess.DEVNULL, check=True)
        xs, fs = parse_dump(DUMP, args.frames, args.stride)
        np.savez(CACHE, x=xs, f=fs)
    natoms = xs.shape[1]
    print(f"frames={len(xs)} atoms={natoms} box=2.9982 nm cutoff={args.cutoff}")

    box = 2.9982
    centers = np.linspace(0.12, args.cutoff, args.rbf)
    width = centers[1] - centers[0]
    ntr = max(2, len(xs) - 5)
    X = design(xs[:ntr], box, args.cutoff, centers, width, natoms)
    Y = fs[:ntr].reshape(-1, 3)
    A = X.T @ X + 1e-6 * np.eye(X.shape[1])
    W = np.linalg.solve(A, X.T @ Y)
    Xte = design(xs[ntr:], box, args.cutoff, centers, width, natoms)
    pred = (Xte @ W).reshape(fs[ntr:].shape)
    err = pred - fs[ntr:]
    rms = float(np.sqrt((err ** 2).mean()))
    mag = float(np.sqrt((fs[ntr:] ** 2).mean()))
    print(f"train frames={ntr} test frames={len(xs)-ntr}")
    print(f"force RMS magnitude   = {mag:8.1f} kJ/mol/nm")
    print(f"model force RMS error = {rms:8.1f} kJ/mol/nm")
    print(f"fraction of variance  = {1 - (rms/mag)**2:.3f}")
    print("note: force RMS error does not map linearly to trajectory error; "
          "small errors still grow through Lyapunov instability")


if __name__ == "__main__":
    main()
