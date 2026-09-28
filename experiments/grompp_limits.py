#!/usr/bin/env python3
"""Which bond does grompp name when it rejects a villin setting?

grompp warns when a bond's estimated oscillation period is shorter than five
timesteps, and fastmode treats that warning as a rejection. This script builds
each h-bonds setting exactly as fastmode.evaluate() does (same mdp builder, same
HMR transform), runs grompp on it, and records the bond grompp names, with the
residue and atom names read from the topology it used.

No simulation runs; each grompp call takes seconds.

Usage:
    python experiments/grompp_limits.py --gmx <gmx> [--out results/grompp_limits_2026-09.txt]
"""
import argparse
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
import fastmode as fm  # noqa: E402

VILLIN = os.path.join(REPO, "systems", "villin")
DTS = [0.004, 0.005]
FACTORS = [None, 2.5, 3.0, 3.5, 4.0, 4.5]   # None = HMR off

WARN = re.compile(r"between\s+atoms\s+(\d+)\s+(\S+)\s+and\s+(\d+)\s+(\S+)\s+has\s+an\s+"
                  r"estimated\s+oscillational\s+period\s+of\s+([0-9.eE+-]+)\s+ps")


def atom_table(top_path):
    """Protein atom number -> (residue number, residue name, atom name)."""
    text = open(top_path).read()
    section = text.split("[ atoms ]")[1].split("[")[0]
    table = {}
    for line in section.splitlines():
        parts = line.split(";")[0].split()
        if len(parts) >= 5 and parts[0].isdigit():
            table[int(parts[0])] = (parts[2], parts[3], parts[4])
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gmx", default=fm.GMX_DEFAULT)
    ap.add_argument("--out", default=os.path.join(REPO, "results", "grompp_limits_2026-09.txt"))
    args = ap.parse_args()

    conf = os.path.join(VILLIN, "eq.gro")
    base_top = os.path.join(VILLIN, "topol.top")
    lines = ["grompp bond check, villin HP-35, constraints = h-bonds",
             "(a setting is rejected when a bond period is below 5 x dt)", "",
             f"  {'dt':>4} {'HMR':>6}  {'grompp':<9} named bond (residue atom - residue atom), period"]
    with tempfile.TemporaryDirectory() as work:
        for dt in DTS:
            for factor in FACTORS:
                spec = dict(fm.REF_SPEC, dt=dt, hmr=factor is not None)
                if factor is not None:
                    spec["hmr_factor"] = factor
                run = os.path.join(work, fm.setting_label(spec))
                os.makedirs(run)
                mdp, _ = fm.build_mdp(spec, 1000, 1.0)
                with open(os.path.join(run, "run.mdp"), "w") as fh:
                    fh.write(mdp)
                top = fm.derive_top(spec, base_top, run)
                if top == base_top:
                    for name, src in fm.local_includes(base_top):
                        os.symlink(src, os.path.join(run, name))
                    top = os.path.join(run, "topol.top")
                    os.symlink(base_top, top)
                g = fm.Gmx(args.gmx, run, 1)
                res = g.grompp(os.path.join(run, "run.mdp"), conf, top, "run")
                text = res.stdout + res.stderr
                hmr = "off" if factor is None else f"{factor:g}"
                m = WARN.search(text)
                if res.returncode == 0:
                    lines.append(f"  {dt * 1000:4.0f} {hmr:>6}  {'accepts':<9}")
                elif m:
                    atoms = atom_table(top)
                    a1, a2 = int(m.group(1)), int(m.group(3))
                    r1, r2 = atoms[a1], atoms[a2]
                    period_fs = float(m.group(5)) * 1000
                    lines.append(f"  {dt * 1000:4.0f} {hmr:>6}  {'rejects':<9} "
                                 f"{r1[1]}{r1[0]} {r1[2]} - {r2[1]}{r2[0]} {r2[2]}, "
                                 f"{period_fs:.0f} fs")
                else:
                    lines.append(f"  {dt * 1000:4.0f} {hmr:>6}  {'rejects':<9} "
                                 f"(other error: {fm.last_error(text)})")
                print(lines[-1], flush=True)
    text = "\n".join(lines) + "\n"
    with open(args.out, "w") as fh:
        fh.write(text)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
