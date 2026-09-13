#!/usr/bin/env python3
"""Coarse-graining probes on the existing reference trajectory.

Intelligent coarse-graining means coarsen the degrees of freedom that cost the
least physics. That is only possible when the local environment is not uniform:
if every water is equivalent, only a uniform map helps. This script measures the
two things that decide whether a non-uniform map has any basis.

  1. Heterogeneity. Distribution of a per-water local order parameter (the
     tetrahedral order q) and the coordination number. A broad distribution means
     the local environments differ, so a per-environment map can specialise.
  2. Spatial correlation. The correlation of the q deviation against distance.
     A correlation length above zero means the environments form coherent
     regions, so coarsening can be done by region, not by random selection.

The trajectory is the fastmode water-small 2 fs reference (200 ps, 2652 atoms).
Only the oxygen sites are used, because they carry the structure.

Usage:
  python3 cgprobe.py [--xtc <path>] [--tpr <path>] [--skip 40] [--rc 0.35]
"""
import argparse
import os
import subprocess
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
GMX = os.path.join(PROJECT, "build", "bin", "gmx")
GMXLIB = os.path.join(PROJECT, "gromacs", "share", "top")
REF = os.path.join(PROJECT, "phase0", "out", "water_small", "runs", "reference",
                   "dt2fs_hmroff_mtsoff_nstlist10_tol0.005")


def extract_gro(xtc, tpr, skip, out):
    env = dict(os.environ, GMXLIB=GMXLIB)
    subprocess.run([GMX, "trjconv", "-f", xtc, "-s", tpr, "-o", out,
                    "-skip", str(skip), "-pbc", "mol"],
                   input="System\n", text=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=True)


def parse_gro_frames(path):
    frames = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        title = lines[i]
        if i + 1 >= len(lines):
            break
        try:
            n = int(lines[i + 1].strip())
        except ValueError:
            i += 1
            continue
        atoms = lines[i + 2:i + 2 + n]
        names = [a[10:15].strip() for a in atoms]
        xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])]
                        for a in atoms])
        box = np.array([float(v) for v in lines[i + 2 + n].split()[:3]])
        frames.append((names, xyz, box))
        i += 3 + n
    return frames


def min_image(d, box):
    d -= box * np.round(d / box)
    return d


def coordination_and_q(xyz, box, rc):
    o = xyz[::3]
    d = o[:, None, :] - o[None, :, :]
    d = min_image(d, box)
    r = np.sqrt((d ** 2).sum(-1))
    np.fill_diagonal(r, np.inf)
    nbr = r < rc
    cn = nbr.sum(1)
    q = np.full(len(o), np.nan)
    for i in range(len(o)):
        js = np.where(nbr[i])[0]
        if len(js) < 4:
            continue
        order = js[np.argsort(r[i, js])[:4]]
        u = d[i, order] / r[i, order][:, None]
        c = u @ u.T
        iu = np.triu_indices(4, 1)
        q[i] = 1.0 - 3.0 / 8.0 * np.sum((c[iu] + 1.0 / 3.0) ** 2)
    return r, o, cn, q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xtc", default=os.path.join(REF, "run.xtc"))
    ap.add_argument("--tpr", default=os.path.join(REF, "run.tpr"))
    ap.add_argument("--skip", type=int, default=40)
    ap.add_argument("--rc", type=float, default=0.35)
    ap.add_argument("--rmax", type=float, default=1.5)
    ap.add_argument("--nbins", type=int, default=30)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="cgprobe_")
    gro = os.path.join(tmp, "ref_sub.gro")
    print(f"extracting every {args.skip}th frame ...")
    extract_gro(args.xtc, args.tpr, args.skip, gro)
    frames = parse_gro_frames(gro)
    print(f"frames={len(frames)} atoms={len(frames[0][0])} "
          f"oxygens={len(frames[0][0]) // 3} rc={args.rc} nm")

    qs, cns = [], []
    corr_num = np.zeros(args.nbins)
    corr_den = np.zeros(args.nbins)
    dq2_sum = 0.0
    lag = (args.rmax / args.nbins, )
    edges = np.linspace(0, args.rmax, args.nbins + 1)
    for names, xyz, box in frames:
        r, o, cn, q = coordination_and_q(xyz, box, args.rc)
        qs.append(q)
        cns.append(cn)
        good = np.isfinite(q)
        dq = np.zeros(len(q))
        dq[good] = q[good] - q[good].mean()
        dq2_sum += (dq[good] ** 2).sum()
        dd = min_image(o[:, None, :] - o[None, :, :], box)
        rr = np.sqrt((dd ** 2).sum(-1))
        prod = dq[:, None] * dq[None, :]
        iu = np.triu_indices(len(o), 1)
        hist, _ = np.histogram(rr[iu], bins=edges,
                               weights=prod[iu])
        cnt, _ = np.histogram(rr[iu], bins=edges)
        corr_num += hist
        corr_den += cnt
    q_all = np.concatenate([q[np.isfinite(q)] for q in qs])
    cn_all = np.concatenate(cns)

    print("\n[1] HETEROGENEITY")
    print(f"  tetrahedral q : mean {q_all.mean():.3f}  std {q_all.std():.3f}  "
          f"min {q_all.min():.3f}  max {q_all.max():.3f}")
    print(f"  coordination  : mean {cn_all.mean():.2f}  std {cn_all.std():.2f}  "
          f"min {cn_all.min()}  max {cn_all.max()}")
    frac_ordered = np.mean(q_all > 0.8)
    print(f"  fraction q>0.8 (ordered)     = {frac_ordered:.3f}")
    print(f"  fraction q<0.5 (disordered)  = {np.mean(q_all < 0.5):.3f}")
    print("  -> a broad distribution means per-environment maps can specialise")

    print("\n[2] SPATIAL CORRELATION of the q deviation")
    g = np.divide(corr_num, corr_den, out=np.zeros_like(corr_num),
                  where=corr_den > 0)
    n_o = len(frames[0][0]) // 3
    norm = dq2_sum / max(1.0, len(frames) * n_o)
    g = g / norm if norm else g
    print("  r (nm)   C(r)")
    for b in range(args.nbins):
        print(f"  {(edges[b] + edges[b + 1]) / 2:5.2f}    {g[b]:+.3f}")
    print("  -> C(r) above zero at large r means coherent regions exist")


if __name__ == "__main__":
    main()
