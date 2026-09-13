#!/usr/bin/env python3
"""Spherical solvent shell test.

Build a non-periodic droplet around the villin protein: keep the protein and the
water within R of it, freeze the outer shell to hold the droplet, and drop PME
for a reaction field. Compare protein observables and speed against the periodic
2 fs reference.

Rationale: the water-order probe showed the protein perturbs water only within
about 0.6 nm, so the far water is bulk. Removing it should save particles. The
same system is then run at a large timestep with constraints = all-bonds to test
whether the two gains stack.

Usage:
  python3 spherical.py build --gro ../phase2b/eq.gro --top ../phase2b/topol.top \
      --out spherical --shell 1.0 --freeze 0.3
  python3 spherical.py run --out spherical
"""
import argparse
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
GMX = os.path.join(PROJECT, "build", "bin", "gmx")
GMXLIB = os.path.join(PROJECT, "gromacs", "share", "top")
KT = 2.4943
NTOMP = 4


def read_gro(path):
    with open(path) as handle:
        lines = handle.readlines()
    n = int(lines[1])
    resname, name, resid, coords = [], [], [], []
    for line in lines[2:2 + n]:
        resname.append(line[5:10].strip())
        name.append(line[10:15].strip())
        resid.append(line[15:20].strip())
        coords.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
    box = [float(x) for x in lines[2 + n].split()[:3]]
    return resname, name, resid, np.array(coords, float), np.array(box, float)


def write_gro(path, title, resname, name, resid, coords, box, vel=None):
    with open(path, "w") as out:
        out.write(title + "\n")
        out.write(f"{len(name):5d}\n")
        for i, (rn, nm, ri, (x, y, z)) in enumerate(zip(resname, name, resid, coords), start=1):
            line = f"{ri:>5}{rn:<5}{nm:>5}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
            out.write(line + "\n")
        out.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


def read_top(path):
    with open(path) as handle:
        return handle.read()


def mdp(dt, constraints, nsteps, coulomb="PME", rcoulomb=1.0, gen_seed=None):
    if coulomb == "PME":
        elec = ["coulombtype = PME", "fourierspacing = 0.12"]
    else:
        elec = ["coulombtype = reaction-field", "epsilon-rf = 78.5"]
    lines = [
        f"dt = {dt}", f"nsteps = {nsteps}",
        "nstxout = 0", "nstvout = 0",
        f"nstxout-compressed = {max(1, round(0.05 / dt))}",
        "nstlog = 1000", "nstenergy = 1000",
        "integrator = md", "constraints = " + constraints,
        "constraint-algorithm = lincs",
        "cutoff-scheme = Verlet", "verlet-buffer-tolerance = 0.005",
    ] + elec + [
        f"rcoulomb = {rcoulomb}", f"rvdw = {rcoulomb}",
        "comm-mode = linear",
        "tcoupl = v-rescale", "tc-grps = Protein SOL Ion", "tau-t = 0.1 0.1 0.1", "ref-t = 300 300 300",
        "pcoupl = no",
        "freezegrps = Boundary", "freezedim = Y Y Y",
    ]
    if gen_seed is None:
        lines.append("gen-vel = no")
    else:
        lines += ["gen-vel = yes", "gen-temp = 300", f"gen-seed = {gen_seed}"]
    return "\n".join(lines) + "\n"


def run(cmd, cwd, stdin=None):
    env = dict(os.environ, GMXLIB=GMXLIB)
    return subprocess.run(cmd, cwd=cwd, env=env, input=stdin,
                          capture_output=True, text=True)


def build(args):
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    resname, name, resid, coords, box = read_gro(args.gro)
    protein = np.array([i for i, rn in enumerate(resname) if rn not in ("SOL", "CL", "NA")])
    water_o = [i for i, (rn, nm) in enumerate(zip(resname, name))
               if rn == "SOL" and nm == "OW"]
    water_mols = [(i, i + 1, i + 2) for i in water_o]
    ions = [i for i, rn in enumerate(resname) if rn in ("CL", "NA")]
    pc = coords[protein]
    d = np.linalg.norm(coords[water_o][:, None, :] - pc[None, :, :], axis=2).min(axis=1)
    keep = [m for m, dd in zip(water_mols, d) if dd <= args.shell]
    frozen = [m for m, dd in zip(water_mols, d) if args.shell - args.freeze < dd <= args.shell]
    mobile = [m for m in keep if m not in frozen]
    print(f"protein {len(protein)} atoms, water {len(water_mols)} molecules, ions {len(ions)}")
    print(f"shell <= {args.shell} nm: keep {len(keep)} waters, freeze {len(frozen)}, mobile {len(mobile)}")
    print(f"new atom count: {len(protein) + 3 * len(keep) + len(ions)}")
    keep_atoms = protein.tolist() + ions
    for m in keep:
        keep_atoms += list(m)
    keep_atoms.sort()
    newres = [resname[i] for i in keep_atoms]
    newname = [name[i] for i in keep_atoms]
    newid = [resid[i] for i in keep_atoms]
    newco = coords[keep_atoms]
    if args.box > 0:
        box = np.array([args.box, args.box, args.box], float)
        center = newco.mean(axis=0)
        newco = (newco + (box / 2.0 - center)) % box
    write_gro(os.path.join(out, "shell.gro"), "spherical shell", newres, newname, newid, newco, box)
    shutil.copy(os.path.join(os.path.dirname(os.path.abspath(args.top)), "posre.itp"),
                os.path.join(out, "posre.itp")) if os.path.exists(
        os.path.join(os.path.dirname(os.path.abspath(args.top)), "posre.itp")) else None
    top = read_top(args.top)
    top = re.sub(r"(?m)^SOL\s+\d+", f"SOL          {len(keep)}", top)
    with open(os.path.join(out, "topol.top"), "w") as handle:
        handle.write(top)
    frozen_atoms = []
    for j, m in enumerate(keep):
        if m in frozen:
            base = keep_atoms.index(m[0]) + 1
            frozen_atoms += [base, base + 1, base + 2]
    newnum = {orig: pos + 1 for pos, orig in enumerate(keep_atoms)}
    natoms = len(keep_atoms)
    water_atoms = [newnum[a] for m in keep for a in m]
    protein_atoms = [newnum[a] for a in protein]
    ion_atoms = [newnum[a] for a in ions]

    def write_group(handle, label, atoms):
        handle.write(f"[ {label} ]\n")
        for k in range(0, len(atoms), 15):
            handle.write(" ".join(str(a) for a in atoms[k:k + 15]) + "\n")

    with open(os.path.join(out, "index.ndx"), "w") as handle:
        write_group(handle, "System", list(range(1, natoms + 1)))
        write_group(handle, "Protein", protein_atoms)
        write_group(handle, "SOL", water_atoms)
        write_group(handle, "Ion", ion_atoms)
        write_group(handle, "Boundary", frozen_atoms)
    print(f"frozen atoms {len(frozen_atoms)}")
    with open(os.path.join(out, "shell.mdp"), "w") as handle:
        handle.write(mdp(args.dt, args.constraints, args.nsteps, args.coulomb))
    print("build done")


def run_plan(args):
    out = os.path.abspath(args.out)
    plans = [("dt2_hbonds", 0.002, "h-bonds", 100000),
             ("dt4_allbonds", 0.004, "all-bonds", 50000),
             ("dt5_allbonds", 0.005, "all-bonds", 40000),
             ("dt6_allbonds", 0.006, "all-bonds", 33334)]
    for tag, dt, cons, nsteps in plans:
        if args.only and tag not in args.only.split(","):
            continue
        wd = os.path.join(out, "runs", tag)
        os.makedirs(wd, exist_ok=True)
        for f in ("shell.gro", "topol.top", "posre.itp", "index.ndx"):
            src = os.path.join(out, f)
            if os.path.exists(src):
                shutil.copy(src, wd)
        with open(os.path.join(wd, "run.mdp"), "w") as handle:
            handle.write(mdp(dt, cons, nsteps, args.coulomb))
        gp = run([GMX, "grompp", "-f", "run.mdp", "-c", "shell.gro", "-p", "topol.top",
                  "-n", "index.ndx", "-o", "run.tpr", "-po", "mdout.mdp"], wd)
        if gp.returncode != 0:
            print(f"{tag}: grompp failed")
            print(gp.stderr[-1500:])
            continue
        md = run([GMX, "mdrun", "-s", "run.tpr", "-deffnm", "run", "-ntmpi", "1",
                  "-ntomp", str(NTOMP), "-nb", "cpu", "-pin", "off", "-resetstep", "2000"], wd)
        if md.returncode != 0:
            print(f"{tag}: mdrun failed")
            print((md.stderr + md.stdout)[-1500:])
            continue
        log = os.path.join(wd, "run.log")
        perf = math.nan
        for line in open(log):
            if line.strip().startswith("Performance:"):
                perf = float(line.split()[1])
        print(f"{tag}: {perf:.1f} ns/day")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--gro", default=os.path.join(PROJECT, "phase2b", "eq.gro"))
    b.add_argument("--top", default=os.path.join(PROJECT, "phase2b", "topol.top"))
    b.add_argument("--out", default=os.path.join(HERE, "spherical"))
    b.add_argument("--shell", type=float, default=1.0)
    b.add_argument("--freeze", type=float, default=0.3)
    b.add_argument("--box", type=float, default=0.0)
    b.add_argument("--coulomb", default="PME")
    b.add_argument("--dt", type=float, default=0.002)
    b.add_argument("--constraints", default="h-bonds")
    b.add_argument("--nsteps", type=int, default=20000)
    b.set_defaults(func=build)
    r = sub.add_parser("run")
    r.add_argument("--out", default=os.path.join(HERE, "spherical"))
    r.add_argument("--only", default=None)
    r.add_argument("--coulomb", default="PME")
    r.set_defaults(func=run_plan)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
