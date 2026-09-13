#!/usr/bin/env python3
"""Three-body descriptor probe.

`forcefit.py` shows a pairwise radial model plateaus at R squared 0.81 on the
water box. The missing part is many-body: angle terms and the SETTLE constraint
force. This script asks whether a cheap angular descriptor closes the gap.

Model (both channels are rotationally equivariant by construction):

    pairwise:  F_i = sum_j c(r_ij) u_ij
    angular:   F_i = sum_j c(r_ij, q_ij) u_ij

where q_ij is a local anisotropy scalar for the neighbour j of atom i:

    q_ij = sum_k cos(theta_jik) * g(r_ik) / n_i

with g a soft cutoff weight and n_i the neighbour count. q measures whether the
local shell is ordered, so the coefficient can respond to structure, not only
distance. This is the cheapest way to add three-body information without a full
symmetry-function machinery.

Usage:
  python3 threebody.py [--frames 20] [--stride 20] [--cutoff 0.5] [--rbf 16]
                       [--qbf 8]
"""
import argparse
import os

import numpy as np

from forcefit import CACHE, build_neighbours


def load():
    z = np.load(CACHE)
    return z["x"], z["f"]


def q_descriptor(x, box, cutoff):
    d, r = build_neighbours(x, box)
    mask = r < cutoff
    u = np.zeros_like(d)
    np.divide(d, r[:, :, None], out=u, where=(r[:, :, None] > 1e-12))
    g = np.exp(-((r / cutoff) ** 2))
    natoms = x.shape[0]
    q = np.zeros((natoms, natoms))
    for i in range(natoms):
        nbr = np.where(mask[i])[0]
        if len(nbr) == 0:
            continue
        uu = u[i, nbr]
        cos = uu @ uu.T
        w = g[i, nbr]
        q[i, nbr] = (cos @ w) / len(nbr)
    return mask, q, u


def rbf_1d(v, centers, width):
    return np.exp(-((v[:, None] - centers[None, :]) / width) ** 2)


def design_pair(x_frames, box, cutoff, centers, width, natoms):
    X = np.zeros((len(x_frames) * natoms, len(centers) * 3))
    for fi, x in enumerate(x_frames):
        d, r = build_neighbours(x, box)
        mask = r < cutoff
        rr = r[mask]
        u = d[mask] / rr[:, None]
        src = np.repeat(np.arange(natoms), mask.sum(1))
        phi = rbf_1d(rr, centers, width)
        off = fi * natoms
        for k in range(len(centers)):
            c = np.zeros((natoms, 3))
            np.add.at(c, src, phi[:, k][:, None] * u)
            X[off:off + natoms, k * 3:k * 3 + 3] = c
    return X


def design_angular(x_frames, box, cutoff, centers, width, qcenters, qwidth,
                   natoms):
    nf = len(centers) * len(qcenters) * 3
    X = np.zeros((len(x_frames) * natoms, nf))
    for fi, x in enumerate(x_frames):
        mask, q, u = q_descriptor(x, box, cutoff)
        d, r = build_neighbours(x, box)
        rr = r[mask]
        qq = q[mask]
        uf = u[mask]
        src = np.repeat(np.arange(natoms), mask.sum(1))
        pr = rbf_1d(rr, centers, width)
        pq = rbf_1d(qq, qcenters, qwidth)
        off = fi * natoms
        col = 0
        for a in range(len(centers)):
            for b in range(len(qcenters)):
                c = np.zeros((natoms, 3))
                np.add.at(c, src, (pr[:, a] * pq[:, b])[:, None] * uf)
                X[off:off + natoms, col:col + 3] = c
                col += 3
    return X


def fit(Xtr, Ytr, Xte, Fte):
    A = Xtr.T @ Xtr + 1e-6 * np.eye(Xtr.shape[1])
    W = np.linalg.solve(A, Xtr.T @ Ytr)
    pred = (Xte @ W).reshape(Fte.shape)
    err = pred - Fte
    rms = float(np.sqrt((err ** 2).mean()))
    mag = float(np.sqrt((Fte ** 2).mean()))
    return W, rms, mag, 1 - (rms / mag) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--stride", type=int, default=20)
    ap.add_argument("--cutoff", type=float, default=0.5)
    ap.add_argument("--rbf", type=int, default=16)
    ap.add_argument("--qbf", type=int, default=8)
    args = ap.parse_args()

    xs, fs = load()
    if len(xs) > args.frames:
        idx = np.linspace(0, len(xs) - 1, args.frames).round().astype(int)
        xs, fs = xs[idx], fs[idx]
    natoms = xs.shape[1]
    box = 2.9982
    ntr = max(2, len(xs) - 5)
    xtr, xts = xs[:ntr], xs[ntr:]
    ftr, fts = fs[:ntr], fs[ntr:]

    centers = np.linspace(0.12, args.cutoff, args.rbf)
    width = centers[1] - centers[0]
    qcenters = np.linspace(-1.0, 1.0, args.qbf)
    qwidth = qcenters[1] - qcenters[0]

    print(f"frames={len(xs)} train={ntr} test={len(xs)-ntr} atoms={natoms} "
          f"cutoff={args.cutoff}")
    Xp = design_pair(xtr, box, args.cutoff, centers, width, natoms)
    _, rms_p, mag, r2_p = fit(Xp, ftr.reshape(-1, 3),
                              design_pair(xts, box, args.cutoff, centers, width,
                                          natoms),
                              fts)
    print(f"pairwise          rms={rms_p:7.1f}  R2={r2_p:5.3f}  "
          f"features={Xp.shape[1]}")

    Xa = design_angular(xtr, box, args.cutoff, centers, width, qcenters,
                        qwidth, natoms)
    _, rms_a, _, r2_a = fit(Xa, ftr.reshape(-1, 3),
                            design_angular(xts, box, args.cutoff, centers, width,
                                           qcenters, qwidth, natoms),
                            fts)
    print(f"pairwise+angular  rms={rms_a:7.1f}  R2={r2_a:5.3f}  "
          f"features={Xa.shape[1]}")
    print(f"force magnitude   rms={mag:7.1f} kJ/mol/nm")
    print(f"R2 gate 0.95: {'PASSED' if r2_a >= 0.95 else 'FAILED'}")


if __name__ == "__main__":
    main()
