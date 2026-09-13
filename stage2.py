#!/usr/bin/env python3
"""Stage 2: virtual-site system build and timestep sweep.

Builds the villin HP-35 domain with hydrogen virtual sites (pdb2gmx -vsite h
-heavyh), solvates and ionises it, minimises it, equilibrates it for 300 ps,
then sweeps dt over {4, 5, 6, 7} fs against a 2 fs reference of the SAME
virtual-site system.  The four validation checks are those of fastmode.py.

Both the reference and the candidates use constraints = all-bonds.  With
h-bonds alone, grompp refuses dt >= 5 fs: the carboxylate CG-OD1 bond has an
oscillational period of 2.2e-02 ps, and no hydrogen transform can lengthen it.

The rebuilt system is not the phase2b system.  Its numbers compare only with
its own 2 fs reference, never with phase2b.

Usage:
    python3 stage2.py build           # pdb2gmx, editconf, solvate, genion, em, eq
    python3 stage2.py sweep [--ref-ps 50]
    python3 stage2.py all
"""

import argparse
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace

import fastmode as fm

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
BUILD = os.path.join(HERE, "vsites", "build")
PDB = os.path.join(PROJECT, "phase2b", "prot.pdb")
EM_MDP = os.path.join(PROJECT, "phase2b", "em.mdp")
IONS_MDP = os.path.join(PROJECT, "phase2b", "ions.mdp")
DTS = [0.004, 0.005, 0.006, 0.007]

EQ_MDP = """integrator      = md
dt              = 0.002
nsteps          = 150000
nstlist         = 20
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
fourierspacing  = 0.12
tcoupl          = v-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.5 0.5
ref-t           = 300 300
pcoupl          = C-rescale
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
constraints     = h-bonds
nstenergy       = 5000
nstlog          = 5000
nstxout-compressed = 5000
gen-vel         = yes
gen-temp        = 300
"""


def run(cmd, **kw):
    print("+ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def build(gmx, ntomp):
    os.makedirs(BUILD, exist_ok=True)
    env = dict(os.environ, GMXLIB=fm.GMXLIB_DEFAULT)
    md = [gmx, "mdrun", "-ntmpi", "1", "-ntomp", str(ntomp),
          "-nb", "cpu", "-pin", "off", "-resetstep", "2000"]
    run([gmx, "pdb2gmx", "-f", PDB, "-vsite", "h", "-heavyh",
         "-ff", "amber99sb", "-water", "tip3p", "-ignh",
         "-o", "prot.gro", "-p", "topol.top", "-i", "posre.itp"],
        cwd=BUILD, env=env, stdin=subprocess.DEVNULL)
    run([gmx, "editconf", "-f", "prot.gro", "-o", "box.gro",
         "-d", "1.2", "-bt", "cubic", "-c"], cwd=BUILD, env=env)
    run([gmx, "solvate", "-cp", "box.gro", "-cs", "spc216.gro",
         "-p", "topol.top", "-o", "solv.gro"], cwd=BUILD, env=env)
    shutil.copy(IONS_MDP, os.path.join(BUILD, "ions.mdp"))
    run([gmx, "grompp", "-f", "ions.mdp", "-c", "solv.gro",
         "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "1"],
        cwd=BUILD, env=env)
    run([gmx, "genion", "-s", "ions.tpr", "-o", "solv_ions.gro",
         "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-neutral"],
        cwd=BUILD, env=env, input=b"SOL\n")
    shutil.copy(EM_MDP, os.path.join(BUILD, "em.mdp"))
    with open(os.path.join(BUILD, "eq.mdp"), "w") as fh:
        fh.write(EQ_MDP)
    run([gmx, "grompp", "-f", "em.mdp", "-c", "solv_ions.gro",
         "-p", "topol.top", "-o", "em.tpr"], cwd=BUILD, env=env)
    run(md + ["-deffnm", "em"], cwd=BUILD, env=env)
    run([gmx, "grompp", "-f", "eq.mdp", "-c", "em.gro",
         "-p", "topol.top", "-o", "eq.tpr"], cwd=BUILD, env=env)
    run(md + ["-deffnm", "eq"], cwd=BUILD, env=env)


def prepare_inputs(out):
    inputs = os.path.join(out, "runs", "inputs")
    os.makedirs(inputs, exist_ok=True)
    top = os.path.join(BUILD, "topol.top")
    shutil.copy(top, os.path.join(inputs, os.path.basename(top)))
    for name, src in fm.local_includes(top):
        shutil.copy(src, os.path.join(inputs, name))
    return os.path.join(inputs, "topol.top")


def sweep(args):
    out = os.path.abspath(args.out)
    top = prepare_inputs(out)
    conf = os.path.join(BUILD, "eq.gro")
    if not os.path.exists(conf):
        sys.exit(f"missing {conf}; run `stage2.py build` first")
    natoms = fm.read_gro_natoms(conf)
    timed = fm.machine_quiet(args.ntomp)
    workdir = os.path.join(out, "runs")
    os.makedirs(workdir, exist_ok=True)
    gmx = fm.Gmx(args.gmx, workdir, args.ntomp)
    print(f"stage2: {natoms} atoms, other load "
          f"{fm.external_load(args.ntomp):.2f}, timed={timed}", flush=True)

    ref_spec = {"dt": 0.002, "hmr": False, "mts": 0, "nstlist": 10,
                "tol": 0.005, "constraints": "all-bonds"}
    ref_dir = os.path.join(workdir, "reference")
    os.makedirs(ref_dir, exist_ok=True)
    print("  running 2 fs reference ...", flush=True)
    ref = fm.evaluate(gmx, ref_spec, conf, top, args.ref_ps, ref_dir,
                      timed)
    if ref.get("rejected") or ref.get("drift") is None:
        sys.exit(f"reference run failed: {ref.get('reason', ref)}")
    ref["pass"] = True
    ref["rdf_dev"] = 0.0
    runs = [ref]

    for dt in DTS:
        spec = dict(ref_spec, dt=dt)
        print(f"  dt {dt * 1000:.0f} fs ...", flush=True)
        rec = fm.evaluate(gmx, spec, conf, top, args.ref_ps,
                          os.path.join(workdir, "dt_sweep"), timed)
        if rec.get("rejected"):
            print(f"    rejected: {rec.get('reason', '')}", flush=True)
        else:
            fm.validate(rec, ref)
            print(f"    {'PASS' if rec['pass'] else 'fail'}"
                  + ("" if rec["pass"]
                     else f" ({'; '.join(rec.get('reasons', []))})"),
                  flush=True)
        runs.append(rec)

    passing = [r for r in runs if r.get("pass") and r is not ref]
    best = (max(passing, key=lambda r: r.get("ns_per_day") or 0.0)
            if passing else None)
    ns = SimpleNamespace(gro=conf, top=top, ref_ps=args.ref_ps,
                         ntomp=args.ntomp, strategy="vsites dt sweep")
    text = fm.write_report(os.path.join(out, "report.txt"), ns, natoms, ref,
                           runs, best)
    print(text, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=["build", "sweep", "all"])
    ap.add_argument("--gmx", default=fm.GMX_DEFAULT)
    ap.add_argument("--ntomp", type=int, default=4)
    ap.add_argument("--ref-ps", type=float, default=200.0)
    ap.add_argument("--out", default=os.path.join(HERE, "vsites"))
    args = ap.parse_args()
    if args.command in ("build", "all"):
        build(args.gmx, args.ntomp)
    if args.command in ("sweep", "all"):
        sweep(args)


if __name__ == "__main__":
    main()
