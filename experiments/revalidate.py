#!/usr/bin/env python3
"""Re-validation experiments for the villin fast-mode settings (September 2026).

An independent review of v0.1.0 found that the protein results rested on single
runs, that the protein-observable check tested a different setting from the one
recommended, and that no run-to-run noise floor existed for the protein
observables. This script closes those gaps. It reuses fastmode.evaluate() and
fastmode.validate() unchanged, so every run is built exactly as the audit
builds it.

It runs, in sequence and on one machine:

  1. Three independent 2 fs references (seeds 11, 12, 13). Their spread is the
     noise floor for every check, including the protein observables.
  2. Three seeds each of plain dt 4 fs (h-bonds, no HMR) and of the
     recommended setting (dt 6 fs, HMR factor 3, all-bonds). These validate the
     recommended setting itself and give timing replicates.
  3. Three seeds of the recommended setting with tighter LINCS (order 8,
     2 iterations), and single runs at dt 7 fs with both LINCS settings.

Every candidate group goes through fastmode.validate() against every reference.
Rerunning on the same --out directory reuses finished runs (no new MD) and
recomputes every observable, so the analysis can be corrected without new
simulations.

Protein observables (radius of gyration, backbone RMSD, per-residue RMSF) and a
production-only energy drift exclude the warmup steps.

Usage:
    python experiments/revalidate.py --gro eq.gro --top topol.top \\
        --out <workdir> --gmx <path/to/gmx> [--ntomp 4]
Writes <out>/revalidate.json and <out>/revalidate.txt.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import fastmode as fm  # noqa: E402

SEEDS = [11, 12, 13]
REF = dict(fm.REF_SPEC)
PLAIN4 = dict(fm.REF_SPEC, dt=0.004)
REC6 = dict(fm.REF_SPEC, dt=0.006, hmr=True, hmr_factor=3.0, constraints="all-bonds")
REC7 = dict(REC6, dt=0.007)
LINCS_TIGHT = ["lincs-order = 8", "lincs-iter = 2"]

_build_mdp = fm.build_mdp


def build_mdp_with_extra(spec, nsteps, ref_ps):
    """fastmode.build_mdp plus optional extra mdp lines carried in the spec."""
    text, nstxout = _build_mdp(spec, nsteps, ref_ps)
    extra = spec.get("extra_mdp") or []
    if extra:
        text += "\n".join(extra) + "\n"
    return text, nstxout


fm.build_mdp = build_mdp_with_extra


def tool(gmx, workdir, args, stdin):
    return subprocess.run([gmx.binary] + args, cwd=workdir, env=gmx.env,
                          input=stdin, text=True, capture_output=True)


def production_drift(gmx, workdir, warmup_ps, natoms):
    """Slope of the conserved energy over production only, per atom."""
    out = os.path.join(workdir, "prod_conserved.xvg")
    tool(gmx, workdir, ["energy", "-f", "run.edr", "-b", f"{warmup_ps:.4f}",
                        "-o", out], "Conserved-En.\n\n")
    cols, _ = fm.read_xvg(out)
    slope = np.polyfit(cols[0], cols[1], 1)[0]
    return float(slope / natoms)


def protein_observables(gmx, workdir, warmup_ps):
    b = f"{warmup_ps:.4f}"
    tool(gmx, workdir, ["gyrate", "-s", "run.tpr", "-f", "run.xtc", "-b", b,
                        "-o", "val_gyr.xvg"], "1\n")
    tool(gmx, workdir, ["rms", "-s", "run.tpr", "-f", "run.xtc", "-b", b,
                        "-o", "val_rms.xvg"], "4\n4\n")
    tool(gmx, workdir, ["rmsf", "-s", "run.tpr", "-f", "run.xtc", "-b", b,
                        "-o", "val_rmsf.xvg", "-res"], "4\n")
    gyr, _ = fm.read_xvg(os.path.join(workdir, "val_gyr.xvg"))
    rms, _ = fm.read_xvg(os.path.join(workdir, "val_rms.xvg"))
    rmsf, _ = fm.read_xvg(os.path.join(workdir, "val_rmsf.xvg"))
    return {"rg": float(np.mean(gyr[1])), "rmsd": float(np.mean(rms[1])),
            "rmsf": [float(x) for x in rmsf[1]]}


def run(gmx, spec, name, conf, top, ref_ps, out, natoms):
    strategy_dir = os.path.join(out, "runs", name)
    os.makedirs(strategy_dir, exist_ok=True)
    load = fm.external_load(gmx.ntomp)
    print(f"  {name} ... other load {load:.2f}", flush=True)
    rec = fm.evaluate(gmx, spec, conf, top, ref_ps, strategy_dir, load <= fm.LOAD_LIMIT)
    rec["name"] = name
    rec["load"] = load
    if rec.get("rejected") or rec.get("drift") is None:
        print(f"    failed: {rec.get('reason', '')}", flush=True)
        return rec
    warmup_ps = fm.warmup_time_ps(spec["dt"])
    rec["drift_production"] = production_drift(gmx, rec["workdir"], warmup_ps, natoms)
    rec["protein"] = protein_observables(gmx, rec["workdir"], warmup_ps)
    print(f"    {rec['ns_per_day']:.1f} ns/day  drift {rec['drift']:.4f} "
          f"(production {rec['drift_production']:.4f})  "
          f"T {rec['temperature']:.2f}  Rg {rec['protein']['rg']:.4f}", flush=True)
    return rec


def rmsf_corr(a, b):
    n = min(len(a), len(b))
    return float(np.corrcoef(a[:n], b[:n])[0, 1])


def rdf_dev(a, b):
    ga, gb = np.asarray(a["rdf_g"]), np.asarray(b["rdf_g"])
    n = min(len(ga), len(gb))
    return float(np.max(np.abs(ga[:n] - gb[:n])))


def summarize(recs, natoms):
    ok = {k: v for k, v in recs.items() if v.get("protein")}
    refs = [ok[f"ref_s{s}"] for s in SEEDS if f"ref_s{s}" in ok]
    lines = ["fastmode re-validation, villin HP-35 (%d atoms)" % natoms,
             "=" * 48, ""]

    def pairs(xs, ys, fn, same):
        vals = []
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                if same and j <= i:
                    continue
                vals.append(fn(x, y))
        return vals

    checks = [
        ("temperature (K)", lambda a, b: abs(a["temperature"] - b["temperature"])),
        ("density (%)", lambda a, b: 100 * abs(a["density"] - b["density"]) / b["density"]),
        ("O-O RDF max dev", rdf_dev),
        ("Rg (nm)", lambda a, b: abs(a["protein"]["rg"] - b["protein"]["rg"])),
        ("RMSF correlation", lambda a, b: rmsf_corr(a["protein"]["rmsf"], b["protein"]["rmsf"])),
    ]
    groups = [("plain dt 4 fs, h-bonds", "plain4"), ("dt 6 fs, HMR f3, all-bonds", "rec6"),
              ("dt 6 fs, ..., LINCS 8/2", "rec6_lincs8")]
    lines.append("Differences from the 2 fs references. The first column is the noise")
    lines.append("floor: every pair of independent 2 fs references. The others compare")
    lines.append("each candidate seed with each reference seed. min / mean / max over pairs.")
    lines.append("")
    head = f"  {'quantity':<18} {'2 fs vs 2 fs':>24}"
    for g, _ in groups:
        head += f" {g:>30}"
    lines.append(head)
    for label, fn in checks:
        row = f"  {label:<18}"
        floor = pairs(refs, refs, fn, True)
        row += f" {min(floor):8.4f}/{np.mean(floor):7.4f}/{max(floor):7.4f}"
        for _, key in groups:
            cands = [ok[f"{key}_s{s}"] for s in SEEDS if f"{key}_s{s}" in ok]
            vals = pairs(cands, refs, fn, False)
            row += f"   {min(vals):8.4f}/{np.mean(vals):7.4f}/{max(vals):7.4f}"
        lines.append(row)
    lines.append("")

    lines.append("Signed mean density difference from the references (%):")
    ref_rho = np.mean([r["density"] for r in refs])
    for g, key in groups:
        cands = [ok[f"{key}_s{s}"]["density"] for s in SEEDS if f"{key}_s{s}" in ok]
        lines.append(f"  {g:<30} {100 * (np.mean(cands) - ref_rho) / ref_rho:+.3f}")
    lines.append("")

    lines.append("Four-check verdicts, each candidate run against each reference seed:")
    names = sorted(n for n in ok if not n.startswith("ref_"))
    for name in names:
        c = ok[name]
        verdicts = []
        for r in refs:
            trial = dict(c)
            trial["reasons"] = []
            fm.validate(trial, {"temperature": r["temperature"],
                                "density": r["density"], "rdf_g": r["rdf_g"]})
            verdicts.append("pass" if trial["pass"] else "FAIL: " + "; ".join(trial["reasons"]))
        n_pass = sum(v == "pass" for v in verdicts)
        lines.append(f"  {name:<18} {n_pass}/{len(refs)}  " + " | ".join(verdicts))
    lines.append("")

    lines.append("Timing: not reported here. These runs were not timed under one")
    lines.append("controlled batch; see experiments/timing.py and results/timing_2026-09.txt.")
    lines.append("")

    lines.append("Energy drift, kJ/mol/ps per atom (whole run / production only):")
    for name in sorted(recs):
        r = recs[name]
        if r.get("drift") is None:
            lines.append(f"  {name:<24} failed: {r.get('reason', '')}")
            continue
        lines.append(f"  {name:<24} {r['drift']:9.5f} / {r['drift_production']:9.5f}")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--ntomp", type=int, default=4)
    ap.add_argument("--ref-ps", type=float, default=200.0)
    ap.add_argument("--summary-only", action="store_true",
                    help="rewrite revalidate.txt from revalidate.json, no runs")
    args = ap.parse_args()
    if args.summary_only:
        out = os.path.abspath(args.out)
        recs = json.load(open(os.path.join(out, "revalidate.json")))
        text = summarize(recs, fm.read_gro_natoms(os.path.abspath(args.gro)))
        with open(os.path.join(out, "revalidate.txt"), "w") as fh:
            fh.write(text)
        print(text)
        return

    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    conf, top = os.path.abspath(args.gro), os.path.abspath(args.top)
    natoms = fm.read_gro_natoms(conf)
    gmx = fm.Gmx(args.gmx, out, args.ntomp)
    print(f"revalidate: {natoms} atoms, {args.ref_ps:.0f} ps per run", flush=True)

    plan = []
    for s in SEEDS:
        plan += [(f"ref_s{s}", dict(REF, seed=s)),
                 (f"plain4_s{s}", dict(PLAIN4, seed=s)),
                 (f"rec6_s{s}", dict(REC6, seed=s))]
    s = SEEDS[0]
    plan += [(f"rec6_lincs8_s{s}", dict(REC6, seed=s, extra_mdp=LINCS_TIGHT)),
             (f"rec7_s{s}", dict(REC7, seed=s)),
             (f"rec7_lincs8_s{s}", dict(REC7, seed=s, extra_mdp=LINCS_TIGHT))]
    plan += [(f"rec6_lincs8_s{s}", dict(REC6, seed=s, extra_mdp=LINCS_TIGHT))
             for s in SEEDS[1:]]

    recs = {}
    path = os.path.join(out, "revalidate.json")
    for name, spec in plan:
        recs[name] = run(gmx, spec, name, conf, top, args.ref_ps, out, natoms)
        with open(path, "w") as fh:
            json.dump(recs, fh, indent=1, default=str)

    text = summarize(recs, natoms)
    with open(os.path.join(out, "revalidate.txt"), "w") as fh:
        fh.write(text)
    print()
    print(text)


if __name__ == "__main__":
    main()
