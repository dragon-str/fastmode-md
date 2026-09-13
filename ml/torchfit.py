#!/usr/bin/env python3
"""Capacity test on the angular descriptor.

Experiment 4 fit the 384-feature pairwise-plus-angular design matrix with a
linear model and reached R squared 0.814, far below the 0.95 gate. Two reasons
are possible: the model has too little capacity, or the descriptor lacks the
information (long-range PME force and the velocity-dependent SETTLE constraint
force). This script separates them. It fits the SAME features with neural
networks of growing width. If capacity were the limit, the test R squared would
rise with width. If it stays flat, the limit is information.

Features are cached in `torchfit_feats.npz`.  Needs the python3.11 venv with
torch at `ml/.venv311`.

Usage:
  .venv311/bin/python torchfit.py [--frames 20] [--epochs 80]
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from forcefit import CACHE
from threebody import design_angular, load

FEATS = os.path.join(HERE, "torchfit_feats.npz")
WIDTHS = [0, 32, 128, 512]


def build_features(frames, cutoff, rbf, qbf):
    xs, fs = load()
    if len(xs) > frames:
        idx = np.linspace(0, len(xs) - 1, frames).round().astype(int)
        xs, fs = xs[idx], fs[idx]
    centers = np.linspace(0.12, cutoff, rbf)
    width = centers[1] - centers[0]
    qcenters = np.linspace(-1.0, 1.0, qbf)
    qwidth = qcenters[1] - qcenters[0]
    X = design_angular(xs, 2.9982, cutoff, centers, width, qcenters, qwidth,
                       xs.shape[1])
    np.savez(FEATS, X=X.astype(np.float32), f=fs.astype(np.float32))
    return X, fs


def load_features(frames, cutoff, rbf, qbf):
    if os.path.exists(FEATS):
        z = np.load(FEATS)
        return z["X"], z["f"]
    return build_features(frames, cutoff, rbf, qbf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch", type=int, default=8192)
    ap.add_argument("--cutoff", type=float, default=0.5)
    ap.add_argument("--rbf", type=int, default=16)
    ap.add_argument("--qbf", type=int, default=8)
    args = ap.parse_args()

    import torch

    torch.manual_seed(0)
    X, fs = load_features(args.frames, args.cutoff, args.rbf, args.qbf)
    nf = X.shape[0] // fs.shape[0]
    ntr = max(2, fs.shape[0] - 5) * nf
    Xtr = torch.tensor(X[:ntr], dtype=torch.float32)
    Ytr = torch.tensor(fs.reshape(-1, 3)[:ntr], dtype=torch.float32)
    Xte = torch.tensor(X[ntr:], dtype=torch.float32)
    Fte = fs.reshape(-1, 3)[ntr:]
    mag = float(np.sqrt((Fte ** 2).mean()))
    mean, std = Xtr.mean(0), Xtr.std(0).clamp_min(1e-6)
    Xtr, Xte = (Xtr - mean) / std, (Xte - mean) / std

    rows = []
    for width in WIDTHS:
        if width == 0:
            A = torch.cat([Xtr, torch.ones(Xtr.shape[0], 1)], 1)
            sol = torch.linalg.lstsq(A, Ytr).solution
            with torch.no_grad():
                pred = torch.cat(
                    [Xte, torch.ones(Xte.shape[0], 1)], 1) @ sol
            pred = pred.numpy()
            rms = float(np.sqrt(((pred - Fte) ** 2).mean()))
            rows.append((width, rms, 1 - (rms / mag) ** 2))
            print(f"width {width:4d}  rms {rms:7.1f}  R2 {1-(rms/mag)**2:5.3f}")
            continue
        model = torch.nn.Sequential(
                torch.nn.Linear(Xtr.shape[1], width), torch.nn.ReLU(),
                torch.nn.Linear(width, width), torch.nn.ReLU(),
                torch.nn.Linear(width, 3))
        opt = torch.optim.Adam(model.parameters(), lr=1e-3)
        lossf = torch.nn.MSELoss()
        n = Xtr.shape[0]
        for ep in range(args.epochs):
            perm = torch.randperm(n)
            for i in range(0, n, args.batch):
                j = perm[i:i + args.batch]
                opt.zero_grad()
                loss = lossf(model(Xtr[j]), Ytr[j])
                loss.backward()
                opt.step()
        with torch.no_grad():
            pred = model(Xte).numpy()
        rms = float(np.sqrt(((pred - Fte) ** 2).mean()))
        r2 = 1 - (rms / mag) ** 2
        rows.append((width, rms, r2))
        print(f"width {width:4d}  rms {rms:7.1f}  R2 {r2:5.3f}")

    write(rows, mag, nf)


def write(rows, mag, nf):
    path = os.path.join(HERE, "torchfit.txt")
    with open(path, "w") as fh:
        fh.write("Capacity test on the pairwise-plus-angular descriptor\n")
        fh.write("water small, 20 frames (15 train, 5 test), 384 features\n")
        fh.write(f"force magnitude {mag:.1f} kJ/mol/nm; features {nf}\n\n")
        fh.write(f"{'hidden':>6} {'rms':>9} {'R2':>7}\n")
        for width, rms, r2 in rows:
            name = "linear" if width == 0 else str(width)
            fh.write(f"{name:>6} {rms:9.1f} {r2:7.3f}\n")
        fh.write("\nR2 gate 0.95. A flat row means capacity is not the limit.\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
