#!/usr/bin/env python3
"""One-site water by iterative Boltzmann inversion (idea 3, coarse-graining).

The force probe closed the position-only atomistic surrogate: long-range PME and
the constraint force are not local functions of the coordinates. Coarse-graining
does not copy forces. It invents an effective pair potential whose purpose is to
reproduce the target structure. That is a different and much weaker demand.

This is a standalone numpy probe, not a GROMACS change. It answers two questions
before any GROMACS table work:

  1. Can one site per water reproduce the atomistic O-O radial distribution?
     Method: iterative Boltzmann inversion, u_new = u_old + kT ln(g/g_target).
  2. What timestep does the resulting soft potential tolerate?
     A soft coarse pair potential has no bond vibration, so dt should rise.

One site per water removes the orientational degrees of freedom and cuts the site
count by 3. The O positions are kept, so the O-O structure is the right target.

Usage:
  python3 cg_ibi.py --n 512 --iters 20 --steps 3000 --dt 0.005
"""
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
RDF = os.path.join(PROJECT, "phase0", "out", "water_small", "runs", "reference",
                   "dt2fs_hmroff_mtsoff_nstlist10_tol0.005", "run_rdf.xvg")
KB = 0.0083144621
T = 300.0
KT = KB * T
MASS = 18.015


def load_target(path, rmax):
    r, g = [], []
    with open(path) as fh:
        for line in fh:
            if line[0] in "@#":
                continue
            p = line.split()
            if len(p) < 2:
                continue
            r.append(float(p[0]))
            g.append(float(p[1]))
    r = np.array(r)
    g = np.array(g)
    m = r <= rmax
    return r[m], g[m]


def interpolate(x, y, xi):
    return np.interp(xi, x, y, left=y[0], right=np.nan)


def force_table(r, u, dr):
    f = -np.gradient(u, dr)
    return f


def run_md(pos, L, u_tab, f_tab, dr, dt, steps, target_g, bins, collect_every,
           gamma, rng, collect=True):
    n = len(pos)
    mass = MASS
    vel = rng.normal(0.0, np.sqrt(KT / mass), size=(n, 3))
    rmax = bins[-1]
    hist = np.zeros(len(bins) - 1)
    ncfg = 0
    for step in range(steps):
        d = pos[:, None, :] - pos[None, :, :]
        d -= L * np.round(d / L)
        dist = np.sqrt((d ** 2).sum(-1))
        np.fill_diagonal(dist, np.inf)
        r = dist.ravel()
        rn = np.where(np.isfinite(r), r, 2 * rmax)
        idx = np.clip((rn / dr).astype(np.int64), 0, len(u_tab) - 1)
        fr = f_tab[idx]
        rinv = np.where(np.isfinite(r) & (r > 0), 1.0 / np.where(r > 0, r, 1.0), 0.0)
        coef = (fr * rinv).reshape(n, n)
        coef = np.where(np.isfinite(dist) & (dist < rmax), coef, 0.0)
        force = np.einsum("ij,ijk->ik", coef, d)
        if collect and step % collect_every == 0:
            m = (r > 0) & (r < rmax) & np.isfinite(r)
            hist += np.histogram(r[m], bins=bins)[0]
            ncfg += 1
        vel += force / mass * dt
        vel += -gamma * vel * dt + np.sqrt(2 * gamma * KT / mass * dt) \
            * rng.normal(size=(n, 3))
        pos = pos + vel * dt
        pos -= L * np.floor(pos / L)
    return hist, ncfg, vel


def rdf_from_hist(hist, n, ncfg, L, bins):
    dr = bins[1] - bins[0]
    r = 0.5 * (bins[:-1] + bins[1:])
    rho = n / L ** 3
    shell = 4 * np.pi * r ** 2 * dr * rho
    g = hist / (ncfg * n * shell)
    return r, g


def smooth(u, width=2):
    k = np.ones(2 * width + 1) / (2 * width + 1)
    return np.convolve(u, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=512)
    ap.add_argument("--rmax", type=float, default=1.0)
    ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--iters", type=int, default=25)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--dt", type=float, default=0.005)
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--damp", type=float, default=0.5)
    ap.add_argument("--core", type=float, default=0.18)
    args = ap.parse_args()

    r_t, g_t = load_target(RDF, args.rmax)
    rho = 884 / 2.99359 ** 3
    L = (args.n / rho) ** (1.0 / 3.0)
    print(f"target: {len(r_t)} points, rmax {args.rmax} nm")
    print(f"density {rho:.4f} /nm^3   n {args.n}   box {L:.4f} nm")

    bins = np.arange(0.0, args.rmax + args.dr, args.dr)
    r_grid = 0.5 * (bins[:-1] + bins[1:])
    g_on_grid = np.interp(r_grid, r_t, g_t, left=0.0, right=1.0)
    g_on_grid = np.maximum(g_on_grid, 1e-6)

    u = -KT * np.log(g_on_grid)
    u = np.minimum(u, 60 * KT)
    u -= u[-1]
    rng = np.random.default_rng(args.seed)
    nside = int(np.ceil(args.n ** (1.0 / 3.0)))
    grid = np.arange(nside) * (L / nside)
    xx, yy, zz = np.meshgrid(grid, grid, grid, indexing="ij")
    pos = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], 1)[:args.n]
    pos += rng.normal(0, 0.02, size=pos.shape)

    print(f"\nIBI iterations (dt {args.dt} ps, {args.steps} steps each)")
    print("  iter   max|g-g_t|   rms|g-g_t|   first_peak   first_min")
    best = None
    for it in range(args.iters):
        f_tab = force_table(r_grid, u, args.dr)
        collect_every = max(1, args.steps // 2000)
        hist, ncfg, vel = run_md(pos, L, u, f_tab, args.dr, args.dt, args.steps,
                                 g_on_grid, bins, collect_every, args.gamma, rng)
        r, g = rdf_from_hist(hist, args.n, ncfg, L, bins)
        gt = np.interp(r, r_t, g_t, left=0.0, right=1.0)
        dev = g - gt
        peak = r[np.argmax(g)]
        first_min = r[np.argmin(np.where((r > peak) & (r < peak + 0.15), g, 1e9))]
        print(f"  {it:4d}   {np.max(np.abs(dev)):9.4f}   "
              f"{np.sqrt(np.mean(dev ** 2)):9.4f}   {peak:8.3f}   {first_min:8.3f}")
        if it < args.iters - 1:
            corr = KT * np.log(np.maximum(g, 1e-6) / np.maximum(gt, 1e-6))
            corr = np.clip(corr, -2 * KT, 2 * KT)
            corr[r < args.core] = 0.0
            corr[r > r_t[-1]] = 0.0
            corr = args.damp * smooth(corr, 3)
            u = u + corr
            u = np.minimum(u, 60 * KT)
            u -= u[-1]
        best = (r.copy(), g.copy())
    np.savetxt(os.path.join(HERE, "cg_ibi_rdf.txt"),
               np.column_stack([r, g, gt]), header="r g_fit g_target")

    print("\nTimestep stability of the final potential")
    print("  dt (ps)   max|g-g_t|   rms|g-g_t|")
    f_tab = force_table(r_grid, u, args.dr)
    for dt in (0.005, 0.010, 0.020, 0.030):
        try:
            hist, ncfg, _ = run_md(pos.copy(), L, u, f_tab, args.dr, dt,
                                   2000, g_on_grid, bins, 10, args.gamma,
                                   np.random.default_rng(11))
            r, g = rdf_from_hist(hist, args.n, ncfg, L, bins)
            gt = np.interp(r, r_t, g_t, left=0.0, right=1.0)
            dev = g - gt
            print(f"  {dt:6.3f}   {np.max(np.abs(dev)):9.4f}   "
                  f"{np.sqrt(np.mean(dev ** 2)):9.4f}")
        except Exception as exc:
            print(f"  {dt:6.3f}   unstable: {exc}")
    print("\n-> one site per water is 3x fewer sites; a stable large dt adds more")


if __name__ == "__main__":
    main()
