#!/usr/bin/env python3
"""Adaptive-resolution landscape on the villin trajectory.

The bulk-water probe showed that tetrahedral order has a correlation length of
about one molecule, so pure water has no coherent regions to partition. Around a
solute the order is organised: it decays with distance from the protein. This
script measures that decay, and the fraction of water that lies beyond a given
distance. Together they set the ceiling of a distance-adaptive coarse-grain map:
coarsen the far water, keep the shell atomistic.

Output: the mean and width of the tetrahedral order q against the distance to the
nearest protein atom, and the water fraction in distance shells.

Usage:
  python3 cgvillin.py [--skip 40]
"""
import argparse
import os
import tempfile

import numpy as np

from cgprobe import (GMX, GMXLIB, HERE, PROJECT, extract_gro, min_image,
                     parse_gro_frames)

REF = os.path.join(PROJECT, "phase0", "out", "villin", "runs", "reference",
                   "dt2fs_hmroff_mtsoff_nstlist10_tol0.005")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xtc", default=os.path.join(REF, "run.xtc"))
    ap.add_argument("--tpr", default=os.path.join(REF, "run.tpr"))
    ap.add_argument("--skip", type=int, default=40)
    ap.add_argument("--rc", type=float, default=0.35)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="cgvillin_")
    gro = os.path.join(tmp, "ref_sub.gro")
    print(f"extracting every {args.skip}th frame ...")
    extract_gro(args.xtc, args.tpr, args.skip, gro)
    frames = parse_gro_frames(gro)
    print(f"frames={len(frames)} atoms={len(frames[0][0])}")

    names, xyz0, _ = frames[0]
    is_water_atom = np.array([n in ("OW", "HW1", "HW2") for n in names])
    is_water = np.array([n == "OW" for n in names])
    is_protein = ~is_water_atom
    nprot = int(is_protein.sum())
    now = int(is_water.sum())
    print(f"protein atoms={nprot} water oxygens={now}")

    bins = np.arange(0.0, 2.61, 0.2)
    q_sum = np.zeros(len(bins) - 1)
    q_sq = np.zeros(len(bins) - 1)
    cnt = np.zeros(len(bins) - 1)
    dn_mean = 0.0
    nlist = []

    for names, xyz, box in frames:
        o = xyz[is_water]
        p = xyz[is_protein]
        d = o[:, None, :] - p[None, :, :]
        d = min_image(d, box)
        r = np.sqrt((d ** 2).sum(-1))
        dn = r.min(1)

        oo = o[:, None, :] - o[None, :, :]
        oo = min_image(oo, box)
        ro = np.sqrt((oo ** 2).sum(-1))
        np.fill_diagonal(ro, np.inf)
        nbr = ro < args.rc
        ro_m = np.where(nbr, ro, np.inf)
        order = np.argpartition(ro_m, 4, axis=1)[:, :4]
        u = np.take_along_axis(oo, order[:, :, None], axis=1)
        u = u / np.take_along_axis(ro_m, order, axis=1)[:, :, None]
        c = np.einsum("ijk,ilk->ijl", u, u)
        sq = (c + 1.0 / 3.0) ** 2
        iu = np.triu_indices(4, 1)
        q = 1.0 - 3.0 / 8.0 * sq[:, iu[0], iu[1]].sum(1)
        q[nbr.sum(1) < 4] = np.nan
        good = np.isfinite(q)
        idx = np.digitize(dn[good], bins) - 1
        for b in range(len(bins) - 1):
            m = idx == b
            q_sum[b] += q[good][m].sum()
            q_sq[b] += (q[good][m] ** 2).sum()
            cnt[b] += m.sum()
        nlist.append(dn)

    print("\nq against distance to nearest protein atom")
    print("  d (nm)    count     q mean   q std")
    for b in range(len(bins) - 1):
        if cnt[b] == 0:
            continue
        mu = q_sum[b] / cnt[b]
        var = max(0.0, q_sq[b] / cnt[b] - mu ** 2)
        print(f"  {bins[b]:.1f}-{bins[b+1]:.1f}  {int(cnt[b]):7d}   "
              f"{mu:+.3f}   {np.sqrt(var):.3f}")

    dn = np.concatenate(nlist)
    print("\nwater fraction beyond a distance")
    print("  d (nm)   fraction beyond")
    for th in (0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0):
        print(f"  {th:4.1f}     {np.mean(dn > th):.3f}")
    print(f"\nmean nearest-protein distance = {dn.mean():.3f} nm")
    print("-> the coarsenable fraction is the fraction beyond the shell radius")


if __name__ == "__main__":
    main()
