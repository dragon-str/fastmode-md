#!/usr/bin/env python3
"""Validate the 12.9x spherical-shell reaction-field setting against a periodic
2 fs reference, with independent replicas.

Runs N replicas of each configuration with fresh velocities (different gen-seed),
computes protein observables (Rg, backbone RMSD, per-residue RMSF) and the water
O-O RDF, and reports the mean and spread across replicas. This gives the 12.9x
setting a statistical footing that a single 200 ps run cannot.

Usage:
    validate10x.py [--seeds 11,22,33,44,55] [--ns 0.3] [--configs name,name]
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import fastmode as fm  # noqa: E402
import spherical as sp  # noqa: E402
import valprotein as vp  # noqa: E402

GMX = os.environ.get("GMX", os.path.join(PROJECT, "build", "bin", "gmx"))
GMXLIB = os.environ.get("GMXLIB", os.path.join(PROJECT, "gromacs", "share", "top"))

EQ = os.path.join(PROJECT, "phase2b")
V3 = os.path.join(HERE, "spherical_v3")

CONFIGS = {
    "periodic_dt2_hbonds": dict(
        kind="periodic", gro=os.path.join(EQ, "eq.gro"), top=os.path.join(EQ, "topol.top"),
        dt=0.002, constraints="h-bonds", coulomb="PME"),
    "shellPME_dt4_allbonds": dict(
        kind="shell", gro=os.path.join(V3, "shell.gro"), top=os.path.join(V3, "topol.top"),
        index=os.path.join(V3, "index.ndx"), dt=0.004, constraints="all-bonds", coulomb="PME"),
    "shellRF_dt4_allbonds": dict(
        kind="shell", gro=os.path.join(V3, "shell.gro"), top=os.path.join(V3, "topol.top"),
        index=os.path.join(V3, "index.ndx"), dt=0.004, constraints="all-bonds", coulomb="reaction-field",
        rcoulomb=1.5),
    "shellRF12_dt4_allbonds": dict(
        kind="shell", gro=os.path.join(V3, "shell.gro"), top=os.path.join(V3, "topol.top"),
        index=os.path.join(V3, "index.ndx"), dt=0.004, constraints="all-bonds", coulomb="reaction-field",
        rcoulomb=1.2),
}

# The shell RDF is normalized by the box volume, not the droplet volume, so a
# shell-vs-periodic RDF comparison is meaningless. Compare a shell variant with
# the same-geometry control instead.
RDF_REF = {"shellRF_dt4_allbonds": "shellPME_dt4_allbonds"}


def run(cmd, cwd, stdin=None):
    env = dict(os.environ, GMXLIB=GMXLIB)
    return subprocess.run(cmd, cwd=cwd, env=env, input=stdin, text=True, capture_output=True)


def setup(d, cfg, seed, nsteps):
    os.makedirs(d, exist_ok=True)
    shutil.copy(cfg["gro"], os.path.join(d, "conf.gro"))
    shutil.copy(cfg["top"], os.path.join(d, "topol.top"))
    if cfg["kind"] == "shell":
        shutil.copy(cfg["index"], os.path.join(d, "index.ndx"))
    for name, src in fm.local_includes(cfg["top"]):
        shutil.copy(src, os.path.join(d, name))
    if cfg["kind"] == "periodic":
        text = fm.build_mdp(dict(fm.REF_SPEC), nsteps, nsteps * cfg["dt"])[0]
        text = text.replace("gen-vel = no",
                            f"gen-vel = yes\ngen-temp = 300\ngen-seed = {seed}")
    else:
        text = sp.mdp(cfg["dt"], cfg["constraints"], nsteps, cfg["coulomb"],
                      rcoulomb=cfg.get("rcoulomb", 1.0), gen_seed=seed)
    with open(os.path.join(d, "run.mdp"), "w") as fh:
        fh.write(text)


def finished(d):
    log = os.path.join(d, "run.log")
    if not os.path.exists(log):
        return False
    with open(log) as fh:
        return "Finished mdrun" in fh.read()


def run_replica(d, cfg, nsteps):
    if finished(d):
        return True
    grompp = [GMX, "grompp", "-f", "run.mdp", "-c", "conf.gro", "-p", "topol.top", "-o", "run.tpr"]
    if cfg["kind"] == "shell":
        grompp += ["-n", "index.ndx"]
    p = run(grompp, d)
    if not os.path.exists(os.path.join(d, "run.tpr")):
        print(f"grompp failed in {d}:\n{p.stdout[-600:]}")
        return False
    p = run([GMX, "mdrun", "-s", "run.tpr", "-o", "run.trr", "-x", "run.xtc",
             "-c", "run.gro", "-e", "run.edr", "-g", "run.log",
             "-ntmpi", "1", "-ntomp", "4", "-nb", "cpu", "-pin", "off",
             "-resetstep", "2000"], d)
    if not finished(d):
        print(f"mdrun failed in {d}:\n{p.stdout[-600:]}")
        return False
    trr = os.path.join(d, "run.trr")
    if os.path.exists(trr):
        os.remove(trr)
    return True


def parse_index(path):
    groups = {}
    name = None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("["):
            name = s.strip("[]").strip()
            groups[name] = []
        elif name:
            groups[name] += [int(x) for x in s.split()]
    return groups


def add_groups(d):
    """Append MobileO and BoundaryO groups (O atoms only) to a shell index."""
    idx = os.path.join(d, "index.ndx")
    groups = parse_index(idx)
    if "MobileO" in groups:
        return
    _, names, _, _, _ = sp.read_gro(os.path.join(d, "conf.gro"))
    ow = [i + 1 for i, n in enumerate(names) if n == "OW"]
    bound = set(groups.get("Boundary", []))
    mob = [i for i in ow if i not in bound]
    bou = [i for i in ow if i in bound]

    def fmt(nm, at):
        out = f"[ {nm} ]\n"
        for j in range(0, len(at), 15):
            out += " ".join(str(x) for x in at[j:j + 15]) + "\n"
        return out

    with open(idx, "a") as fh:
        fh.write(fmt("MobileO", mob) + fmt("BoundaryO", bou))


def rdf(d, warmup, group="name OW"):
    cmd = [GMX, "rdf", "-s", "run.tpr", "-f", "run.xtc", "-o", "rdf.xvg",
           "-ref", group, "-sel", group, "-rmax", "1.0", "-bin", "0.01",
           "-b", f"{warmup:.1f}"]
    if os.path.exists(os.path.join(d, "index.ndx")):
        cmd += ["-n", "index.ndx"]
    run(cmd, d)
    cols, _ = fm.read_xvg(os.path.join(d, "rdf.xvg"))
    return cols


def com_contact(d):
    """Protein COM range and minimum protein-boundary O distance over the run."""
    run([GMX, "traj", "-s", "run.tpr", "-f", "run.xtc", "-n", "index.ndx",
         "-com", "-ox", "com.xvg"], d, stdin="1\n")
    cols, _ = fm.read_xvg(os.path.join(d, "com.xvg"))
    com = np.array(cols[1:4]).T if len(cols) >= 4 else np.array([cols[1]]).T
    com_range = float(np.max(np.linalg.norm(com - com[0], axis=1)))
    run([GMX, "mindist", "-s", "run.tpr", "-f", "run.xtc", "-n", "index.ndx",
         "-od", "mindist.xvg"], d, stdin="1\n6\n")
    cols, _ = fm.read_xvg(os.path.join(d, "mindist.xvg"))
    dist = np.array(cols[1])
    return com_range, float(np.min(dist)), float(np.mean(dist < 0.35))



def water_order(d, dt, ns, nframes=10):
    xtc = "run.xtc"
    tpr = "run.tpr"
    out = "frames.gro"
    skip = max(1, int(round(ns * 1000.0 / 0.05 / nframes)))
    run([GMX, "trjconv", "-f", xtc, "-s", tpr, "-o", out,
         "-skip", str(skip), "-pbc", "mol"], d, stdin="System\n")
    lines = open(os.path.join(d, out)).read().splitlines()
    i = 0
    qs, nns = [], []
    while i < len(lines):
        try:
            n = int(lines[i + 1].strip())
        except (ValueError, IndexError):
            i += 1
            continue
        atoms = lines[i + 2:i + 2 + n]
        names = [a[10:15].strip() for a in atoms]
        xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])]
                        for a in atoms])
        box = np.array([float(v) for v in lines[i + 2 + n].split()[:3]])
        o = xyz[[j for j, nm in enumerate(names) if nm == "OW"]]
        if len(o) > 8:
            dm = o[:, None, :] - o[None, :, :]
            dm -= box * np.round(dm / box)
            r = np.sqrt((dm ** 2).sum(-1))
            np.fill_diagonal(r, np.inf)
            nns.append(np.min(r, axis=1))
            order = np.argsort(r, axis=1)[:, :4]
            u = dm[np.arange(len(o))[:, None], order] / r[np.arange(len(o))[:, None], order][:, :, None]
            c = u @ u.transpose(0, 2, 1)
            iu = np.triu_indices(4, 1)
            q = 1.0 - 3.0 / 8.0 * np.sum((c[:, iu[0], iu[1]] + 1.0 / 3.0) ** 2, axis=1)
            qs.append(q)
        i += 3 + n
    if not qs:
        return None
    return float(np.mean(np.concatenate(qs))), float(np.std(np.concatenate(qs))), \
        float(np.mean(np.concatenate(nns)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", default="11,22,33,44,55")
    ap.add_argument("--ns", type=float, default=0.3)
    ap.add_argument("--configs", default=",".join(CONFIGS))
    ap.add_argument("--out", default=os.path.join(HERE, "validate10x"))
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    names = [s for s in args.configs.split(",") if s]
    res = {}

    os.makedirs(args.out, exist_ok=True)
    state = dict(system="villin spherical-shell validation", seeds=seeds, ns=args.ns,
                 configs=names, out=args.out, started=time.time(),
                 done={}, current=None)

    def save_state():
        tmp = os.path.join(args.out, "validate.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(state, fh, indent=1)
        os.replace(tmp, os.path.join(args.out, "validate.json"))

    save_state()

    for name in names:
        cfg = CONFIGS[name]
        nsteps = int(round(args.ns * 1000.0 / cfg["dt"]))
        rg, rms, rmsf, rdfg, rdfr, drift = [], [], [], [], [], []
        comrng, mind, frac = [], [], []
        wo_q, wo_s, wo_nn = [], [], []
        for seed in seeds:
            d = os.path.join(args.out, name, f"seed_{seed}")
            print(f"[{name}] seed {seed}", flush=True)
            state["current"] = dict(config=name, seed=seed, phase="setup",
                                    started=time.time())
            save_state()
            setup(d, cfg, seed, nsteps)
            if cfg["kind"] == "shell":
                add_groups(d)
            state["current"]["phase"] = "mdrun"
            save_state()
            if not run_replica(d, cfg, nsteps):
                print(f"[{name}] seed {seed} FAILED", flush=True)
                state["done"][f"{name}/seed_{seed}"] = dict(ok=False)
                state["current"] = None
                save_state()
                continue
            state["current"]["phase"] = "analysis"
            save_state()
            g, r, rmsf_i = vp.observables(d)
            rg.append(float(np.mean(g)))
            rms.append(float(np.mean(r)))
            rmsf.append(rmsf_i)
            drift.append(fm.conserved_drift(os.path.join(d, "run.log")))
            group = "MobileO" if cfg["kind"] == "shell" else "name OW"
            cols = rdf(d, min(50.0, args.ns * 250.0), group)
            rdfr.append(cols[0])
            rdfg.append(cols[1])
            if cfg["kind"] == "shell":
                cr, mc, fc = com_contact(d)
                comrng.append(cr)
                mind.append(mc)
                frac.append(fc)
            wo = water_order(d, cfg["dt"], args.ns)
            if wo is not None:
                wo_q.append(wo[0])
                wo_s.append(wo[1])
                wo_nn.append(wo[2])
            state["done"][f"{name}/seed_{seed}"] = dict(
                ok=True, rg=float(np.mean(g)), rms=float(np.mean(r)),
                drift=drift[-1], q=(wo[0] if wo is not None else None))
            state["current"] = None
            save_state()
        res[name] = dict(rg=rg, rms=rms, rmsf=rmsf, rdfg=rdfg, rdfr=rdfr, drift=drift,
                         comrng=comrng, mind=mind, frac=frac,
                         wo_q=wo_q, wo_s=wo_s, wo_nn=wo_nn)
        print(f"[{name}] {len(rg)} replicas", flush=True)

    def peak(r, g):
        m = (r >= 0.22) & (r <= 0.35)
        i = int(np.argmax(g[m]))
        return float(r[m][i]), float(g[m][i])

    ref = res.get("periodic_dt2_hbonds")
    ref_peak = None
    if ref and ref["rdfg"]:
        rr = np.mean([a[:100] for a in ref["rdfr"]], axis=0)
        gg = np.mean([a[:100] for a in ref["rdfg"]], axis=0)
        ref_peak = peak(rr, gg)

    lines = []
    lines.append("validate10x: independent-seed replicas, %.2f ns each, seeds %s"
                 % (args.ns, args.seeds))
    lines.append("")
    lines.append(f"{'config':<26} {'n':>2} {'Rg nm':>16} {'RMSD nm':>16} {'RMSF corr':>10} "
                 f"{'RMSF maxdev':>12} {'O-O r':>7} {'q':>7} {'q std':>7} {'nn nm':>7} "
                 f"{'COM rng':>8} {'mind':>7} {'frac':>7} {'drift':>10}")
    for name in names:
        r = res[name]
        if not r["rg"]:
            lines.append(f"{name:<26}  0  (no replicas)")
            continue
        rg_m, rg_s = np.mean(r["rg"]), np.std(r["rg"])
        rm_m, rm_s = np.mean(r["rms"]), np.std(r["rms"])
        if ref and name != "periodic_dt2_hbonds" and r["rmsf"] and ref["rmsf"]:
            n = min([len(a) for a in r["rmsf"]] + [len(a) for a in ref["rmsf"]])
            m = np.mean([a[:n] for a in r["rmsf"]], axis=0)
            mref = np.mean([a[:n] for a in ref["rmsf"]], axis=0)
            corr = float(np.corrcoef(m, mref)[0, 1])
            maxdev = float(np.max(np.abs(m - mref)))
        else:
            corr, maxdev = 1.0, 0.0
        if r["rdfg"]:
            rr = np.mean([a[:100] for a in r["rdfr"]], axis=0)
            gg = np.mean([a[:100] for a in r["rdfg"]], axis=0)
            pr, pg = peak(rr, gg)
        else:
            pr, pg = float("nan"), float("nan")
        dr = np.mean([d for d in r["drift"] if d is not None]) if r["drift"] else float("nan")
        cr = np.mean(r["comrng"]) if r["comrng"] else float("nan")
        mi = np.mean(r["mind"]) if r["mind"] else float("nan")
        fc = np.mean(r["frac"]) if r["frac"] else float("nan")
        qm = np.mean(r["wo_q"]) if r["wo_q"] else float("nan")
        qsd = np.mean(r["wo_s"]) if r["wo_s"] else float("nan")
        qnn = np.mean(r["wo_nn"]) if r["wo_nn"] else float("nan")
        lines.append(f"{name:<26} {len(r['rg']):>2} {rg_m:>8.4f}+-{rg_s:<6.4f} {rm_m:>8.4f}+-{rm_s:<6.4f} "
                     f"{corr:>10.3f} {maxdev:>12.4f} {pr:>7.3f} {qm:>7.3f} {qsd:>7.3f} {qnn:>7.4f} "
                     f"{cr:>8.4f} {mi:>7.4f} {fc:>7.3f} {dr:>10.5f}")
    if ref_peak:
        lines.append("")
        lines.append("periodic O-O first peak: r %.3f nm, g %.3f" % ref_peak)
        lines.append("Shell RDF uses mobile O only (MobileO group); its g(r) is not "
                     "box-normalized, so compare the peak, not the plateau.")

    text = "\n".join(lines) + "\n"
    out = os.path.join(HERE, "validate10x.txt") if args.out == os.path.join(HERE, "validate10x") else os.path.join(args.out, "validate10x.txt")
    with open(out, "w") as fh:
        fh.write(text)
    state["current"] = None
    state["finished"] = time.time()
    state["summary"] = {name: dict(
        n=len(res[name]["rg"]),
        rg=float(np.mean(res[name]["rg"])) if res[name]["rg"] else None,
        rms=float(np.mean(res[name]["rms"])) if res[name]["rms"] else None,
        q=float(np.mean(res[name]["wo_q"])) if res[name]["wo_q"] else None,
        drift=float(np.mean([d for d in res[name]["drift"] if d is not None]))
        if any(d is not None for d in res[name]["drift"]) else None)
        for name in names}
    save_state()
    print(text)


if __name__ == "__main__":
    main()
