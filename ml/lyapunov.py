#!/usr/bin/env python3
"""Measure the Lyapunov horizon of the water trajectory.

A learned propagator replaces several MD steps with one model call. Its error is
a perturbation. If that error grows exponentially, the useful horizon of any such
model is short, and the model must be near exact. This script measures the growth
rate directly: it runs two trajectories from the same start, with one atom
displaced by 0.001 nm, and tracks the RMS per-atom separation against time.

The negative result that matters: the time for the separation to reach 0.1 nm
(one bond length, structure lost). It bounds every learned-step idea, whatever
the model family.

Usage:
  python3 lyapunov.py [--ps 20] [--delta 0.001] [--ntomp 4]
"""
import argparse
import os
import re
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
GMX = os.path.join(PROJECT, "build", "bin", "gmx")
GMXLIB = os.path.join(PROJECT, "gromacs", "share", "top")
GRO = os.path.join(PROJECT, "phase0", "systems", "npt_3.0.gro")
REF_MDP = os.path.join(PROJECT, "phase0", "out", "water_small", "runs",
                       "reference", "dt2fs_hmroff_mtsoff_nstlist10_tol0.005",
                       "run.mdp")
WORK = os.path.join(HERE, "lyapunov")


def perturb(gro_in, gro_out, delta, mol=(0, 1, 2)):
    with open(gro_in) as fh:
        lines = fh.read().splitlines()
    for i in mol:
        line = lines[2 + i]
        x = float(line[20:28]) + delta
        lines[2 + i] = line[:20] + f"{x:8.3f}" + line[28:]
    with open(gro_out, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def make_mdp(ps, nstxout, dt=0.002):
    with open(REF_MDP) as fh:
        text = fh.read()
    nsteps = round(ps / dt)
    out = []
    for line in text.splitlines():
        key = line.split("=")[0].strip()
        if key == "nsteps":
            line = f"nsteps          = {nsteps}"
        elif key == "nstxout":
            line = f"nstxout         = {nstxout}"
        elif key == "nstxout-compressed":
            line = "nstxout-compressed = 0"
        elif key == "nstvout":
            line = "nstvout         = 0"
        elif key == "nstenergy":
            line = "nstenergy       = 100000"
        elif key == "nstlog":
            line = "nstlog          = 100000"
        out.append(line)
    return "\n".join(out) + "\n"


def run(cmd, cwd):
    env = dict(os.environ, GMXLIB=GMXLIB)
    r = subprocess.run(cmd, cwd=cwd, env=env, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise SystemExit(f"failed: {' '.join(cmd)}\n{r.stdout[-2000:]}")
    return r.stdout


def parse_trr(path, natoms):
    env = dict(os.environ, GMXLIB=GMXLIB)
    r = subprocess.run([GMX, "dump", "-f", path], env=env, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    pat = re.compile(r"^\s*x\[\s*(\d+)\]=\{\s*([-0-9.eE+]+),\s*([-0-9.eE+]+),"
                     r"\s*([-0-9.eE+]+)\}")
    frames = []
    cur = np.zeros((natoms, 3))
    count = 0
    for line in r.stdout.splitlines():
        if line.strip().startswith("x (") and "x3" in line:
            cur = np.zeros((natoms, 3))
            count = 0
            frames.append(cur)
            continue
        m = pat.search(line)
        if m and frames:
            i = int(m.group(1))
            cur[i] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
            count += 1
    return frames


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ps", type=float, default=20.0)
    ap.add_argument("--delta", type=float, default=0.001)
    ap.add_argument("--nstxout", type=int, default=20)
    ap.add_argument("--ntomp", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(WORK, exist_ok=True)
    natoms = int(open(GRO).read().splitlines()[1])
    perturb(GRO, os.path.join(WORK, "start_pert.gro"), args.delta)
    mdp = make_mdp(args.ps, args.nstxout)
    with open(os.path.join(WORK, "run.mdp"), "w") as fh:
        fh.write(mdp)

    for tag, start in (("ref", GRO), ("pert", os.path.join(WORK, "start_pert.gro"))):
        run([GMX, "grompp", "-f", "run.mdp", "-c", start, "-p",
             os.path.join(PROJECT, "phase0", "systems", "topol_3.0.top"),
             "-o", f"{tag}.tpr"], cwd=WORK)
        run([GMX, "mdrun", "-s", f"{tag}.tpr", "-deffnm", tag,
             "-ntmpi", "1", "-ntomp", str(args.ntomp), "-nb", "cpu",
             "-pin", "off", "-resetstep", "2000"], cwd=WORK)

    ref = parse_trr(os.path.join(WORK, "ref.trr"), natoms)
    pert = parse_trr(os.path.join(WORK, "pert.trr"), natoms)
    n = min(len(ref), len(pert))
    if n < 3:
        raise SystemExit(f"too few frames: ref {len(ref)} pert {len(pert)}")

    box = 2.99820
    times = np.array([args.nstxout * 0.002 * k for k in range(n)])
    rms = np.zeros(n)
    for k in range(n):
        d = ref[k] - pert[k]
        d -= box * np.round(d / box)
        rms[k] = np.sqrt(np.mean(np.sum(d * d, axis=1)))

    print("n_frames", n, "natoms", natoms)
    for k in range(0, n, max(1, n // 10)):
        print(f"  t={times[k]:6.1f} ps  rms={rms[k]:.3e} nm")

    good = (rms > 5e-4) & (rms < 0.05)
    lam = float("nan")
    if good.sum() >= 3:
        p = np.polyfit(times[good], np.log(rms[good]), 1)
        lam = p[0]
    lines = [
        "Lyapunov horizon probe",
        f"system: water small, {natoms} atoms, dt 2 fs, {args.ps} ps, "
        f"rigid water molecule translated by {args.delta} nm",
        "",
        "  time(ps)    rms(nm)",
    ]
    for k in range(n):
        lines.append(f"  {times[k]:7.1f}   {rms[k]:.4e}")
    lines.append("")
    if np.isfinite(lam) and lam > 0:
        lines.append(f"growth rate lambda = {lam:.3f} /ps  "
                     f"(Lyapunov time {1/lam:.3f} ps)")
        for thr in (0.03, 0.1):
            t = np.log(thr / args.delta) / lam
            lines.append(f"  rms reaches {thr} nm at t = {t:.2f} ps")
        lines.append("")
        lines.append(
            "Interpretation: a learned step whose error is a random position "
            "perturbation of size e has a useful horizon of "
            "ln(0.1/e)/lambda ps, before the structure is lost.")
    else:
        lines.append("no exponential growth observed in the measurement window")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(HERE, "lyapunov.txt"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
