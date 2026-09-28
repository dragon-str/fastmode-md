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
sys.path.insert(0, HERE)
import fastmode as fm  # noqa: E402

GMX = fm.GMX_DEFAULT
GMXLIB = fm.GMXLIB_DEFAULT


def run_tool(rundir, args, stdin="", out=""):
    env = dict(os.environ)
    if GMXLIB:
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


def observables(rundir, begin_ps=0.0):
    b = ["-b", f"{begin_ps:.4f}"]
    run_tool(rundir, ["gyrate", "-s", "run.tpr", "-f", "run.xtc", "-o", "val_gyr.xvg"] + b, "1\n", "val_gyr.xvg")
    run_tool(rundir, ["rms", "-s", "run.tpr", "-f", "run.xtc", "-o", "val_rms.xvg"] + b, "4\n4\n", "val_rms.xvg")
    run_tool(rundir, ["rmsf", "-s", "run.tpr", "-f", "run.xtc", "-o", "val_rmsf.xvg", "-res"] + b, "4\n", "val_rmsf.xvg")
    gyr, _ = fm.read_xvg(os.path.join(rundir, "val_gyr.xvg"))
    rms, _ = fm.read_xvg(os.path.join(rundir, "val_rms.xvg"))
    rmsf, _ = fm.read_xvg(os.path.join(rundir, "val_rmsf.xvg"))
    return gyr[1], rms[1], rmsf[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--begin-ps", type=float, default=0.0,
                    help="skip this much of each trajectory (the warmup); "
                         "results/valprotein.txt used 0")
    args = ap.parse_args()

    ref_rg, ref_rms, ref_rmsf = observables(args.ref, args.begin_ps)
    print(f"reference {args.ref}")
    print(f"  Rg mean {ref_rg.mean():.4f} std {ref_rg.std():.4f} nm")
    print(f"  backbone RMSD mean {ref_rms.mean():.4f} std {ref_rms.std():.4f} nm")
    print()
    print(f"{'run':<52} {'Rg nm':>14} {'RMSD nm':>14} {'RMSF corr':>10} {'RMSF maxdev':>11}")
    for d in args.runs:
        rg, rms, rmsf = observables(d, args.begin_ps)
        n = min(len(rmsf), len(ref_rmsf))
        corr = float(np.corrcoef(rmsf[:n], ref_rmsf[:n])[0, 1])
        maxdev = float(np.max(np.abs(rmsf[:n] - ref_rmsf[:n])))
        name = os.path.basename(d)
        print(f"{name:<52} {rg.mean():.4f}+-{rg.std():.4f} {rms.mean():.4f}+-{rms.std():.4f} {corr:10.3f} {maxdev:11.4f}")


if __name__ == "__main__":
    main()
