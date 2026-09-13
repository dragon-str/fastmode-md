#!/usr/bin/env python3
"""Ceiling check for a learned force or propagator on a CPU.

Experiment 1 measures the wall-clock cost of one GROMACS step. It then measures
the cost of a NumPy neural-network forward pass over the same atoms. The ratio
tells us how many steps a learned propagator must cover before it can pay for
itself.

Experiment 2 fits the measured step cost against particle count. It uses that
fit to estimate the speedup from a coarse-grained model with fewer particles.

Usage:
  python3 ceiling.py [--tpr PATH] [--steps 300] [--ntomp 4]
"""
import argparse
import glob
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE0 = os.path.dirname(HERE)
PROJECT = os.path.dirname(PHASE0)
GMX = os.path.join(PROJECT, "build", "bin", "gmx")

WATER_SMALL = os.path.join(
    PHASE0, "out", "water_small", "runs", "runs",
    "dt5fs_hmroff_mtsoff_nstlist10_tol0.005", "run.tpr")
VILLIN = os.path.join(
    PHASE0, "out", "villin", "runs", "runs",
    "dt4fs_hmroff_mtsoff_nstlist10_tol0.005", "run.tpr")


def time_mdrun(tpr, steps, ntomp, tag):
    out = os.path.join(HERE, tag)
    cmd = [GMX, "mdrun", "-s", tpr, "-deffnm", out, "-nsteps", str(steps),
           "-ntmpi", "1", "-ntomp", str(ntomp), "-nb", "cpu", "-pin", "off"]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        sys.exit("mdrun failed: " + proc.stderr[-2000:])
    log = out + ".log"
    text = open(log).read()
    m = re.search(r"Performance:\s+([0-9.]+)", text)
    nsday = float(m.group(1)) if m else float("nan")
    dt = None
    m = re.search(r"dt\s+=\s+([0-9.]+)", text)
    if m:
        dt = float(m.group(1))
    return {"wall_s": wall, "ns_per_day": nsday, "dt_ps": dt,
            "t_step_s": wall / steps}


def numpy_net_timing(n_atoms, in_dim, hidden, out_dim, iters=30):
    rng = np.random.default_rng(0)
    x = rng.standard_normal((n_atoms, in_dim)).astype(np.float32)
    W1 = rng.standard_normal((in_dim, hidden)).astype(np.float32) * 0.05
    b1 = np.zeros(hidden, np.float32)
    W2 = rng.standard_normal((hidden, hidden)).astype(np.float32) * 0.05
    b2 = np.zeros(hidden, np.float32)
    W3 = rng.standard_normal((hidden, out_dim)).astype(np.float32) * 0.05
    b3 = np.zeros(out_dim, np.float32)

    def forward():
        h = np.tanh(x @ W1 + b1)
        h = np.tanh(h @ W2 + b2)
        return h @ W3 + b3

    forward()
    ts = []
    for _ in range(iters):
        t0 = time.perf_counter()
        forward()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def experiment1(tpr, steps, ntomp):
    print("=== Experiment 1: step cost vs net cost ===")
    md = time_mdrun(tpr, steps, ntomp, "bench_md")
    text = open(os.path.join(HERE, "bench_md.log")).read()
    m = re.search(r"There are:\s+(\d+) Atoms", text)
    n_atoms = int(m.group(1)) if m else None
    if not n_atoms:
        sys.exit("could not find the atom count in bench_md.log")
    md["t_step_s"] = (md["dt_ps"] / 1000.0) * 86400.0 / md["ns_per_day"]
    print(f"mdrun: {n_atoms} atoms, {steps} steps, "
          f"{md['wall_s']:.2f} s wall, {md['ns_per_day']:.1f} ns/day, "
          f"{md['t_step_s']*1e3:.4f} ms/step")
    rows = []
    for hidden in (16, 32, 64, 128, 256):
        t = numpy_net_timing(n_atoms, 32, hidden, 3)
        need = t / md["t_step_s"]
        rows.append((hidden, t, need))
        print(f"  net in=32 hidden={hidden:<3d} out=3 : "
              f"{t*1e3:9.3f} ms/call  break-even K = {need:12.1f} steps")
    return md, n_atoms, rows


def load_system_stats():
    stats = []
    for path in sorted(glob.glob(os.path.join(PHASE0, "out", "*", "results.json"))):
        d = json.load(open(path))
        b = d.get("best")
        if not b:
            continue
        dt = b["spec"]["dt"]
        nsday = b["ns_per_day"]
        t_step = (dt / 1000.0) * 86400.0 / nsday
        stats.append((d["system"], d["atoms"], dt, nsday, t_step))
    return stats


def experiment2():
    print("\n=== Experiment 2: step cost vs particle count ===")
    stats = load_system_stats()
    name = {"water_small": "water small", "water_medium": "water medium",
            "water_large": "water large", "villin": "villin"}
    N = np.array([s[1] for s in stats], float)
    T = np.array([s[4] for s in stats], float)
    for s in stats:
        print(f"  {name.get(s[0], s[0]):<14s} {s[1]:6d} atoms  "
              f"{s[2]*1000:.0f} fs  {s[3]:7.1f} ns/day  {s[4]*1e3:8.4f} ms/step")
    A = np.vstack([np.log(N), np.ones_like(N)]).T
    slope, intercept = np.linalg.lstsq(A, np.log(T), rcond=None)[0]
    print(f"  fit: t_step ~ N^{slope:.3f}  (cost per step grows with N^{slope:.2f})")
    print("  a 4:1 coarse-grain at fixed step count would give:")
    for ratio in (2, 4, 8):
        pred = math.exp(intercept) * (N[0] / ratio) ** slope
        print(f"    {ratio}:1 fewer atoms -> {T[0]/pred:.2f}x fewer ms/step "
              f"for the small box")
    return stats, slope


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpr", default=WATER_SMALL)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--ntomp", type=int, default=4)
    args = ap.parse_args()
    md, n_atoms, rows = experiment1(args.tpr, args.steps, args.ntomp)
    stats, slope = experiment2()
    with open(os.path.join(HERE, "ceiling.txt"), "w") as fh:
        fh.write(f"mdrun t_step = {md['t_step_s']*1e3:.4f} ms/step "
                 f"at {n_atoms} atoms\n")
        for hidden, t, need in rows:
            fh.write(f"net hidden={hidden} {t*1e3:.3f} ms/call "
                     f"break-even K={need:.1f}\n")
        fh.write(f"step-cost exponent = {slope:.3f}\n")
    print(f"\nwrote {os.path.join(HERE, 'ceiling.txt')}")


if __name__ == "__main__":
    main()
