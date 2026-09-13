"""Soft-boundary spherical solvent shell.

Same idea as spherical.py, but the outer water layer is not frozen. It is a
separate moleculetype (BWA) with harmonic position restraints, so it is
thermostatted and can fluctuate. This removes the rigid 0 K wall.

Build: protein + mobile SOL + restrained BWA + ions.
Run: grompp/mdrun with reaction-field or PME.
"""
import argparse
import os
import re
import shutil
import subprocess

import numpy as np

import spherical as sp

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
GMX = os.environ.get("GMX", os.path.join(PROJECT, "build", "bin", "gmx"))
GMXLIB = os.environ.get("GMXLIB", os.path.join(PROJECT, "gromacs", "share", "top"))

BOUNDARY_ITP = """[ moleculetype ]
; name  nrexcl
BWA     2

[ atoms ]
;   nr  type  resnr residue  atom   cgnr    charge     mass
     1   OW       1    BWA     OW      1   -0.834   16.00000
     2   HW       1    BWA    HW1      1    0.417    1.00800
     3   HW       1    BWA    HW2      1    0.417    1.00800

[ settles ]
; OW    funct   doh     dhh
1       1       0.09572 0.15139

[ exclusions ]
1       2       3
2       1       3
3       1       2

[ position_restraints ]
; atom  type      fx      fy      fz
1       1      {k:g}     {k:g}     {k:g}
2       1      {k:g}     {k:g}     {k:g}
3       1      {k:g}     {k:g}     {k:g}
"""


def run(cmd, cwd, stdin=None):
    env = dict(os.environ)
    env["GMXLIB"] = GMXLIB
    return subprocess.run(cmd, cwd=cwd, input=stdin, text=True,
                          capture_output=True, env=env)


def water_molecules(resname, name):
    mols = []
    i = 0
    while i < len(name):
        if resname[i] in ("SOL", "HOH", "WAT") and name[i] == "OW":
            mols.append((i, i + 1, i + 2))
            i += 3
        else:
            i += 1
    return mols


def build(args):
    resname, name, resid, coords, box = sp.read_gro(args.gro)
    resname = np.asarray(resname)
    name = np.asarray(name)

    protein = [i for i in range(len(name)) if resname[i] not in ("SOL", "CL", "NA")]
    ions = [i for i in range(len(name)) if resname[i] in ("CL", "NA")]
    waters = water_molecules(list(resname), list(name))

    prot = coords[protein]
    info = []
    for (o, h1, h2) in waters:
        d = coords[o] - prot
        d -= box * np.round(d / box)
        info.append(float(np.sqrt((d ** 2).sum(1)).min()))

    mobile, boundary = [], []
    for (o, h1, h2), d in zip(waters, info):
        if d <= args.shell - args.freeze:
            mobile.append((o, h1, h2))
        elif d <= args.shell:
            boundary.append((o, h1, h2))

    order = list(protein)
    for m in mobile:
        order += list(m)
    for m in boundary:
        order += list(m)
    order += list(ions)

    out = coords[order].copy()
    rn = [resname[i] for i in order]
    nm = [name[i] for i in order]
    rd = [resid[i] for i in order]

    newbox = box
    if args.box:
        newbox = np.array([args.box] * 3)
        out = (out + newbox / 2 - out.mean(0)) % newbox

    os.makedirs(args.out, exist_ok=True)
    sp.write_gro(os.path.join(args.out, "soft.gro"), "soft shell",
                 rn, nm, rd, out, newbox)

    with open(os.path.join(args.out, "boundary.itp"), "w") as fh:
        fh.write(BOUNDARY_ITP.format(k=args.k))

    src = os.path.join(PROJECT, "phase2b", "topol.top")
    if args.top:
        src = args.top
    top = open(src).read()
    top = top.replace('#include "amber99sb.ff/tip3p.itp"',
                      '#include "amber99sb.ff/tip3p.itp"\n#include "boundary.itp"')
    lines = top.splitlines()
    for i, ln in enumerate(lines):
        if re.match(r"^SOL\s+\d+", ln):
            lines[i] = f"SOL          {len(mobile)}"
            lines.insert(i + 1, f"BWA          {len(boundary)}")
            break
    top = "\n".join(lines) + "\n"
    with open(os.path.join(args.out, "topol.top"), "w") as fh:
        fh.write(top)
    for f in ("posre.itp",):
        s = os.path.join(PROJECT, "phase2b", f)
        if os.path.exists(s):
            shutil.copy(s, os.path.join(args.out, f))

    def write_group(fh, label, atoms):
        fh.write(f"[ {label} ]\n")
        for k in range(0, len(atoms), 15):
            fh.write(" ".join(str(a) for a in atoms[k:k + 15]) + "\n")

    with open(os.path.join(args.out, "index.ndx"), "w") as fh:
        prot_ix = list(range(1, len(protein) + 1))
        sol_ix = list(range(len(protein) + 1, len(protein) + 3 * len(mobile) + 1))
        bwa_ix = list(range(sol_ix[-1] + 1, sol_ix[-1] + 3 * len(boundary) + 1))
        ion_ix = list(range(bwa_ix[-1] + 1, bwa_ix[-1] + len(ions) + 1))
        write_group(fh, "Protein", prot_ix)
        write_group(fh, "SOL", sol_ix)
        write_group(fh, "BWA", bwa_ix)
        write_group(fh, "Ion", ion_ix)
        write_group(fh, "Boundary", bwa_ix)
        write_group(fh, "System", list(range(1, len(order) + 1)))

    print(f"mobile {len(mobile)} boundary {len(boundary)} atoms {len(order)} box {newbox[0]:.2f}")
    print("build done")


def mdp_text(dt, constraints, nsteps, coulomb, rcoulomb, gen_seed):
    text = sp.mdp(dt, constraints, nsteps, coulomb, rcoulomb=rcoulomb, gen_seed=gen_seed)
    text = text.replace("tc-grps = Protein SOL Ion", "tc-grps = Protein SOL BWA Ion")
    text = text.replace("tau-t = 0.1 0.1 0.1", "tau-t = 0.1 0.1 0.1 0.1")
    text = text.replace("ref-t = 300 300 300", "ref-t = 300 300 300 300")
    text = re.sub(r"^freezegrps.*\n", "", text, flags=re.M)
    text = re.sub(r"^freezedim.*\n", "", text, flags=re.M)
    return text


def run_one(args, tag, dt, constraints, nsteps, coulomb, rcoulomb, gen_seed=None):
    d = os.path.join(args.out, "runs", tag)
    os.makedirs(d, exist_ok=True)
    for f in ("soft.gro", "topol.top", "index.ndx", "boundary.itp", "posre.itp"):
        p = os.path.join(args.out, f)
        if os.path.exists(p):
            shutil.copy(p, os.path.join(d, f))
    md = mdp_text(dt, constraints, nsteps, coulomb, rcoulomb, gen_seed)
    with open(os.path.join(d, "run.mdp"), "w") as fh:
        fh.write(md)
    g = run([GMX, "grompp", "-f", "run.mdp", "-c", "soft.gro", "-r", "soft.gro",
             "-p", "topol.top", "-n", "index.ndx", "-o", "run.tpr",
             "-po", "mdout.mdp"], d)
    if g.returncode != 0:
        print(f"{tag}: grompp failed")
        print((g.stderr or "").strip().splitlines()[-3:])
        return False
    m = run([GMX, "mdrun", "-s", "run.tpr", "-deffnm", "run", "-ntmpi", "1",
             "-ntomp", str(args.ntomp), "-nb", "cpu", "-pin", "off",
             "-resetstep", "2000"], d)
    log = os.path.join(d, "run.log")
    perf = None
    if os.path.exists(log):
        for ln in open(log):
            if ln.startswith("Performance:"):
                perf = float(ln.split()[1])
    if perf:
        print(f"{tag}: {perf:.1f} ns/day")
    else:
        err = (m.stderr or "").strip().splitlines()
        print(f"{tag}: run failed; {err[-1] if err else 'no Performance line'}")
    return perf is not None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--gro", default=os.path.join(PROJECT, "phase2b", "eq.gro"))
    b.add_argument("--top", default=None)
    b.add_argument("--out", default=os.path.join(HERE, "softshell"))
    b.add_argument("--shell", type=float, default=1.2)
    b.add_argument("--freeze", type=float, default=0.3)
    b.add_argument("--box", type=float, default=8.0)
    b.add_argument("--k", type=float, default=1000.0)
    b.set_defaults(func=build)

    r = sub.add_parser("run")
    r.add_argument("--out", default=os.path.join(HERE, "softshell"))
    r.add_argument("--coulomb", default="reaction-field")
    r.add_argument("--rcoulomb", type=float, default=1.2)
    r.add_argument("--ntomp", type=int, default=4)
    r.add_argument("--only", default=None)
    r.set_defaults(func=None)

    args = ap.parse_args()
    if args.cmd == "build":
        build(args)
        return
    plans = [("dt4_allbonds", 0.004, "all-bonds", 50000),
             ("dt5_allbonds", 0.005, "all-bonds", 40000),
             ("dt6_allbonds", 0.006, "all-bonds", 33334)]
    only = args.only.split(",") if args.only else None
    for tag, dt, cons, nsteps in plans:
        if only and tag not in only:
            continue
        run_one(args, tag, dt, cons, nsteps, args.coulomb, args.rcoulomb)


if __name__ == "__main__":
    main()
