#!/usr/bin/env python3
"""One-site coarse-grain water with a learned 3-body (angular) term.

Pairwise IBI (`cg_ibi.py`) diverged: a one-site isotropic pair potential cannot
place the O-O first peak at 0.28 nm and the second shell at once. This probe adds
three-body terms and asks whether that closure can reproduce the atomistic O-O
radial distribution.

Potential:
    U = sum_{i<j} u2(r_ij)
      + sum_i sum_{j<k} w(r_ij) w(r_ik) u3(cos theta_jik)
with a Gaussian envelope w(r) = exp(-((r-0.28)/0.05)^2) that localises the
angular term to the first coordination shell. Both u2(r) and u3(cos theta) are
tables, updated by iterative Boltzmann inversion against targets measured from
the atomistic trajectory.

Usage: python3 cg3.py --n 256 --iters 8 --steps 1500
"""
import argparse
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("bgmod", os.path.join(HERE, "bg", "bg.py"))
bg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bg)
BUILD = os.path.join(HERE, "bg", "build")
KB = 0.0083144621
T = 300.0
KT = KB * T
MASS = 18.015
WC, WSIG = 0.28, 0.05


def envelope(r):
    return np.exp(-((r - WC) / WSIG) ** 2)


def denvelope(r):
    w = envelope(r)
    return np.where(np.isfinite(r), w * (-2.0 * (r - WC) / WSIG ** 2), 0.0)


def load_targets(path, rmax, dr, cgrid, stride=5):
    frames = bg.read_frames(path)[::stride]
    box = frames[0][1]
    L = box[0]
    o = np.array([c[0::3, :] for c, b in frames])
    nf = len(o)
    n = o.shape[1]
    bins = np.arange(0.0, rmax + dr, dr)
    gh = np.zeros(len(bins) - 1)
    ch = np.zeros(len(cgrid) - 1)
    for f in range(nf):
        d = o[f][:, None, :] - o[f][None, :, :]
        d -= L * np.round(d / L)
        dist = np.sqrt((d ** 2).sum(-1))
        np.fill_diagonal(dist, np.inf)
        rr = dist.ravel()
        m = (rr > 0) & (rr < rmax)
        gh += np.histogram(rr[m], bins=bins)[0]
        within = (dist > 0) & (dist < 0.45)
        e = d / np.where(dist > 0, dist, 1.0)[..., None]
        for i in range(n):
            nb = np.nonzero(within[i])[0]
            if len(nb) < 2:
                continue
            jj, kk = np.triu_indices(len(nb), 1)
            j = nb[jj]
            k = nb[kk]
            c = (e[i, j] * e[i, k]).sum(1)
            wgt = envelope(dist[i, j]) * envelope(dist[i, k])
            ch += np.histogram(c, bins=cgrid, weights=wgt)[0]
    r = 0.5 * (bins[:-1] + bins[1:])
    rho = n / L ** 3
    shell = 4 * np.pi * r ** 2 * dr * rho
    g = gh / (nf * n * shell)
    c = 0.5 * (cgrid[:-1] + cgrid[1:])
    P = ch / max(ch.sum(), 1e-30)
    return r, g, c, P, L, rho


def pair_force_table(r, u, dr):
    return -np.gradient(u, dr)


def run_md(pos, L, rgrid, f2, dr, cedges, ccent, u3, du3, dt, steps, collect_every,
           gamma, rng, collect=True):
    n = len(pos)
    vel = rng.normal(0.0, np.sqrt(KT / MASS), size=(n, 3))
    rmax = rgrid[-1]
    gh = np.zeros(len(rgrid))
    ch = np.zeros(len(cedges) - 1)
    ncfg = 0
    bins = np.concatenate([[0.0], rgrid + 0.5 * dr])
    for step in range(steps):
        d = pos[:, None, :] - pos[None, :, :]
        d -= L * np.round(d / L)
        dist = np.sqrt((d ** 2).sum(-1))
        np.fill_diagonal(dist, np.inf)
        rr = dist.ravel()
        rn = np.where(np.isfinite(rr), rr, 2 * rmax)
        idx = np.clip((rn / dr).astype(np.int64), 0, len(f2) - 1)
        rinv = np.where(np.isfinite(rr) & (rr > 0), 1.0 / np.where(rr > 0, rr, 1.0), 0.0)
        coef = (f2[idx] * rinv).reshape(n, n)
        coef = np.where(dist < rmax, coef, 0.0)
        force = np.einsum("ij,ijk->ik", coef, d)
        e = d / np.where(dist > 0, dist, 1.0)[..., None]
        w = envelope(dist)
        dw = denvelope(dist)
        within = (dist > 0) & (dist < 0.45)
        for i in range(n):
            nb = np.nonzero(within[i])[0]
            if len(nb) < 2:
                continue
            jj, kk = np.triu_indices(len(nb), 1)
            j = nb[jj]
            k = nb[kk]
            rij = dist[i, j]
            rik = dist[i, k]
            eij = e[i, j]
            eik = e[i, k]
            c = (eij * eik).sum(1)
            u3c = np.interp(c, ccent, u3)
            du3c = np.interp(c, ccent, du3)
            wij = w[i, j]
            wik = w[i, k]
            dwij = dw[i, j]
            dwik = dw[i, k]
            gc_j = (eik - c[:, None] * eij) / rij[:, None]
            gc_k = (eij - c[:, None] * eik) / rik[:, None]
            gj = dwij[:, None] * wik[:, None] * u3c[:, None] * eij \
                + wij[:, None] * wik[:, None] * du3c[:, None] * gc_j
            gk = dwik[:, None] * wij[:, None] * u3c[:, None] * eik \
                + wij[:, None] * wik[:, None] * du3c[:, None] * gc_k
            gi = -dwij[:, None] * wik[:, None] * u3c[:, None] * eij \
                - dwik[:, None] * wij[:, None] * u3c[:, None] * eik \
                - wij[:, None] * wik[:, None] * du3c[:, None] * (gc_j + gc_k)
            np.add.at(force, i, gi.sum(0))
            np.add.at(force, j, gj)
            np.add.at(force, k, gk)
            if collect and step % collect_every == 0:
                ch += np.histogram(c, bins=cedges, weights=wij * wik)[0]
        if collect and step % collect_every == 0:
            m = (rr > 0) & (rr < rmax) & np.isfinite(rr)
            gh += np.histogram(rr[m], bins=bins)[0]
            ncfg += 1
        vel += force / MASS * dt
        vel += -gamma * vel * dt + np.sqrt(2 * gamma * KT / MASS * dt) \
            * rng.normal(size=(n, 3))
        pos = pos + vel * dt
        pos -= L * np.floor(pos / L)
    return gh, ch, ncfg


def rdf_from_hist(gh, n, ncfg, L, rgrid, dr):
    rho = n / L ** 3
    shell = 4 * np.pi * rgrid ** 2 * dr * rho
    return gh / (ncfg * n * shell)


def smooth(u, width=2):
    k = np.ones(2 * width + 1) / (2 * width + 1)
    return np.convolve(u, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--rmax", type=float, default=0.9)
    ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--iters", type=int, default=8)
    ap.add_argument("--steps", type=int, default=1500)
    ap.add_argument("--dt", type=float, default=0.005)
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--damp", type=float, default=0.3)
    args = ap.parse_args()

    cgrid = np.linspace(-1.0, 1.0, 61)
    r, g_t, c, P_t, Lbox, rho = load_targets(
        os.path.join(BUILD, "traj.gro"), args.rmax, args.dr, cgrid)
    print(f"target from traj: {len(r)} r points, box {Lbox:.3f}, rho {rho:.3f}, "
          f"{len(c)} angle bins")
    L = (args.n / rho) ** (1.0 / 3.0)
    rgrid = np.arange(args.dr, args.rmax + args.dr, args.dr)
    bins = np.concatenate([[0.0], rgrid + 0.5 * args.dr])
    g_grid = np.maximum(np.interp(rgrid, r, g_t, left=0.0, right=1.0), 1e-6)
    P_grid = np.maximum(np.interp(c, c, P_t), 1e-12)

    u2 = -KT * np.log(g_grid)
    u2 = np.minimum(u2, 60 * KT)
    u2 -= u2[-1]
    u3 = np.zeros_like(c)

    rng = np.random.default_rng(args.seed)
    nside = int(np.ceil(args.n ** (1.0 / 3.0)))
    grid = np.arange(nside) * (L / nside)
    xx, yy, zz = np.meshgrid(grid, grid, grid, indexing="ij")
    pos = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], 1)[:args.n]
    pos += rng.normal(0, 0.02, size=pos.shape)

    print(f"\nIBI (pair + 3-body), n {args.n}, box {L:.3f}, dt {args.dt}")
    print("  iter  max|g-gt|  rms|g-gt|  max|P-Pt|")
    for it in range(args.iters):
        f2 = pair_force_table(rgrid, u2, args.dr)
        du3 = np.gradient(u3, c)
        collect_every = max(1, args.steps // 500)
        gh, ch, ncfg = run_md(pos, L, rgrid, f2, args.dr, cgrid, c, u3, du3,
                              args.dt, args.steps, collect_every, args.gamma, rng)
        g = rdf_from_hist(gh, args.n, ncfg, L, rgrid, args.dr)
        P = ch / max(ch.sum(), 1e-30)
        gt = np.interp(rgrid, r, g_t, left=0.0, right=1.0)
        dev = g - gt
        print(f"  {it:4d}  {np.max(np.abs(dev)):9.4f}  "
              f"{np.sqrt(np.mean(dev ** 2)):9.4f}  {np.max(np.abs(P - P_grid)):9.4f}")
        if it < args.iters - 1:
            corr2 = KT * np.log(np.maximum(g, 1e-6) / np.maximum(gt, 1e-6))
            corr2 = np.clip(corr2, -2 * KT, 2 * KT)
            corr2[rgrid < 0.18] = 0.0
            corr2 = args.damp * smooth(corr2, 3)
            u2 = np.minimum(u2 + corr2, 60 * KT)
            u2 -= u2[-1]
            corr3 = KT * np.log(np.maximum(P, 1e-12) / np.maximum(P_grid, 1e-12))
            corr3 = np.clip(corr3, -2 * KT, 2 * KT)
            corr3 = args.damp * smooth(corr3, 2)
            u3 = u3 + corr3
            m = (c > -0.99) & (c < 0.99)
            np.save(os.path.join(HERE, "cg3_state.npy"),
                    {"r": rgrid, "u2": u2, "c": c, "u3": u3}, allow_pickle=True)

    # timestep stability of the learned potential
    print("\nTimestep stability")
    f2 = pair_force_table(rgrid, u2, args.dr)
    du3 = np.gradient(u3, c)
    for dt in (0.005, 0.010, 0.020, 0.030):
        gh, ch, ncfg = run_md(pos.copy(), L, rgrid, f2, args.dr, cgrid, c, u3, du3,
                              dt, 1000, 10, args.gamma,
                              np.random.default_rng(11))
        g = rdf_from_hist(gh, args.n, ncfg, L, rgrid, args.dr)
        gt = np.interp(rgrid, r, g_t, left=0.0, right=1.0)
        dev = g - gt
        print(f"  dt {dt:5.3f}  max|g-gt| {np.max(np.abs(dev)):8.4f}  "
              f"rms {np.sqrt(np.mean(dev ** 2)):8.4f}")


if __name__ == "__main__":
    main()
