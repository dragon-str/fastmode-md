#!/usr/bin/env python3
"""Compare protein observables between a 2 fs reference and fast-mode winners.

Checks that the all-bonds and HMR settings used to raise the timestep do not
change the protein ensemble. Reports radius of gyration, backbone RMSD and
per-residue RMSF, each against the reference run.

Usage:
    valprotein.py --ref <rundir> <rundir> [<rundir> ...]
"""

import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import fastmode as fm  # noqa: E402

GMX = os.environ.get("GMX", os.path.join(PROJECT, "build", "bin", "gmx"))
GMXLIB = os.environ.get("GMXLIB", os.path.join(PROJECT, "gromacs", "share", "top"))


def run_tool(rundir, args, stdin="", out=""):
    env = dict(os.environ)
    env["GMXLIB"] = GMXLIB
    p = subprocess.run(
        [GMX] + args,
        cwd=rundir,
        env=env,
        input=stdin,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if out and not os.path.exists(os.path.join(rundir, out)):
        raise RuntimeError(f"{args[0]} failed in {rundir}:\n{p.stdout[-800:]}")
    return p.stdout


def observables(rundir):
    run_tool(rundir, ["gyrate", "-s", "run.tpr", "-f", "run.xtc", "-o", "val_gyr.xvg"], "1\n", "val_gyr.xvg")
    run_tool(rundir, ["rms", "-s", "run.tpr", "-f", "run.xtc", "-o", "val_rms.xvg"], "4\n4\n", "val_rms.xvg")
    run_tool(rundir, ["rmsf", "-s", "run.tpr", "-f", "run.xtc", "-o", "val_rmsf.xvg", "-res"], "4\n", "val_rmsf.xvg")
    gyr, _ = fm.read_xvg(os.path.join(rundir, "val_gyr.xvg"))
    rms, _ = fm.read_xvg(os.path.join(rundir, "val_rms.xvg"))
    rmsf, _ = fm.read_xvg(os.path.join(rundir, "val_rmsf.xvg"))
    return gyr[1], rms[1], rmsf[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True)
    ap.add_argument("runs", nargs="+")
    args = ap.parse_args()

    ref_rg, ref_rms, ref_rmsf = observables(args.ref)
    print(f"reference {args.ref}")
    print(f"  Rg mean {ref_rg.mean():.4f} std {ref_rg.std():.4f} nm")
    print(f"  backbone RMSD mean {ref_rms.mean():.4f} std {ref_rms.std():.4f} nm")
    print()
    print(f"{'run':<52} {'Rg nm':>14} {'RMSD nm':>14} {'RMSF corr':>10} {'RMSF maxdev':>11}")
    for d in args.runs:
        rg, rms, rmsf = observables(d)
        n = min(len(rmsf), len(ref_rmsf))
        corr = float(np.corrcoef(rmsf[:n], ref_rmsf[:n])[0, 1])
        maxdev = float(np.max(np.abs(rmsf[:n] - ref_rmsf[:n])))
        name = os.path.basename(d)
        print(f"{name:<52} {rg.mean():.4f}+-{rg.std():.4f} {rms.mean():.4f}+-{rms.std():.4f} {corr:10.3f} {maxdev:11.4f}")


if __name__ == "__main__":
    main()
