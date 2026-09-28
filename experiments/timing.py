#!/usr/bin/env python3
"""Interleaved timing replicates for settings already built by revalidate.py.

A shared machine changes absolute ns/day between batches (the same 2 fs villin
reference measured 46.8 and 69.6 ns/day in two earlier batches). This script
times every setting once per round, for several rounds, in the same order, so
each setting sees the same background load. The speedup of a setting is then
its ns/day divided by the reference's ns/day in the same round, and the report
gives the median and range of that ratio across rounds.

Each timing run re-uses the setting's run.tpr, runs a fixed number of steps
with mdrun's counters reset after the warmup, and reads the Performance line.

Usage:
    python experiments/timing.py --runs <revalidate out>/runs --gmx <gmx> \\
        [--rounds 5] [--nsteps 10000] [--ntomp 4]
Writes timing.json and timing.txt in the --out directory (default: --runs/..).
"""
import argparse
import glob
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import fastmode as fm  # noqa: E402

SETTINGS = [
    ("2 fs reference", "ref_s11"),
    ("dt 4 fs, h-bonds", "plain4_s11"),
    ("dt 6 fs, HMR f3, all-bonds, default LINCS", "rec6_s11"),
    ("dt 6 fs, HMR f3, all-bonds, LINCS order 8 iter 2", "rec6_lincs8_s11"),
    ("dt 7 fs, HMR f3, all-bonds, LINCS order 8 iter 2", "rec7_lincs8_s11"),
]


def tpr_for(runs, name):
    hits = glob.glob(os.path.join(runs, name, "*", "run.tpr"))
    if len(hits) != 1:
        sys.exit(f"expected one run.tpr under {name}, found {len(hits)}")
    return hits[0]


def time_once(gmx, tpr, workdir, nsteps, ntomp):
    os.makedirs(workdir, exist_ok=True)
    deffnm = os.path.join(workdir, "t")
    for ext in (".log", ".edr", ".xtc", ".gro", ".cpt"):
        if os.path.exists(deffnm + ext):
            os.remove(deffnm + ext)
    load = os.getloadavg()[0]
    subprocess.run([gmx, "mdrun", "-s", tpr, "-deffnm", deffnm,
                    "-nsteps", str(nsteps), "-resetstep", str(fm.WARMUP_STEPS),
                    "-ntmpi", "1", "-ntomp", str(ntomp), "-nb", "cpu",
                    "-pin", "off"],
                   env=fm.gmx_env(), capture_output=True, text=True)
    return fm.performance_ns_per_day(deffnm + ".log"), load


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--rounds", type=int, default=5)
    ap.add_argument("--nsteps", type=int, default=10000)
    ap.add_argument("--ntomp", type=int, default=4)
    args = ap.parse_args()

    out = os.path.abspath(args.out or os.path.dirname(os.path.abspath(args.runs)))
    tprs = [(label, name, tpr_for(args.runs, name)) for label, name in SETTINGS]
    table = {name: [] for _, name, _ in tprs}
    loads = []
    for r in range(args.rounds):
        for label, name, tpr in tprs:
            perf, load = time_once(args.gmx, tpr, os.path.join(out, "timing", name),
                                   args.nsteps, args.ntomp)
            table[name].append(perf)
            loads.append(load)
            print(f"round {r + 1} {name:<18} {perf} ns/day  (1-min load {load:.2f})",
                  flush=True)

    ref = table["ref_s11"]
    lines = ["Interleaved timing, villin HP-35, 4 OpenMP threads",
             f"{args.rounds} rounds x {len(tprs)} settings, {args.nsteps} steps each,"
             f" counters reset after {fm.WARMUP_STEPS} steps",
             f"1-minute load average over the batch: {min(loads):.2f} to {max(loads):.2f}",
             "",
             f"  {'setting':<50} {'median ns/day':>13} {'speedup median':>15} {'range':>14}"]
    summary = {}
    for label, name, _ in tprs:
        ratios = [p / q for p, q in zip(table[name], ref)]
        med = statistics.median(ratios)
        summary[name] = {"label": label, "ns_per_day": table[name], "ratios": ratios,
                         "median_ratio": med}
        lines.append(f"  {label:<50} {statistics.median(table[name]):13.1f} "
                     f"{med:14.2f}x {min(ratios):6.2f}-{max(ratios):.2f}x")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(out, "timing.txt"), "w") as fh:
        fh.write(text)
    with open(os.path.join(out, "timing.json"), "w") as fh:
        json.dump({"rounds": args.rounds, "nsteps": args.nsteps, "loads": loads,
                   "settings": summary}, fh, indent=1)
    print()
    print(text)


if __name__ == "__main__":
    main()
