#!/usr/bin/env python3
"""search: brute-force a large setting space for the fastest safe GROMACS run.

fastmode sweeps a few axes by hand.  This tool widens the space and lets it run
for hours.  The space is the hydrogen mass factor, the constraint set, the
timestep and an optional global mass scale.  The objective is the largest
ns/day that still passes the four fastmode checks.

Two stages, because the checks need different run lengths:

  1. SCREEN.  A short run (default 20 ps) rejects crashes, grompp refusals and
     energy drift.  It cannot resolve temperature or the RDF, so it is a filter,
     not a verdict.
  2. CERTIFY.  The fastest screens by ns/day are re-run at production length
     (default 200 ps) and judged by the full fastmode validator.

State goes to `search.json` after every trial, so the dashboard reads live
progress and the run can resume after a stop.  A candidate already recorded in
either stage is skipped.

The mass factor is a real integration device (HMR at factor 4 is standard).  The
global mass scale is a pure frequency shift: it slows the dynamics, so a speedup
that comes from it is not the same physics.  It is off by default and recorded
apart when used.
"""

import argparse
import itertools
import json
import os
import shutil
import sys
import time

import fastmode as fm


def parse_list(text, cast=float):
    return [cast(x) for x in text.split(",") if x.strip() != ""]


def candidates(args, protein):
    dts = parse_list(args.dts, float)
    factors = parse_list(args.hmr_factors, float)
    scales = parse_list(args.mass_scales, float)
    for dt, factor, cons, scale in itertools.product(
            dts, factors, args.constraints, scales):
        hmr = factor > 0 and protein
        if not protein and (factor > 0 or scale != 1.0):
            continue
        spec = dict(fm.REF_SPEC)
        spec.update(dt=dt, hmr=hmr, constraints=cons, mass_scale=scale)
        if hmr:
            spec["hmr_factor"] = factor
        yield spec


def screen_verdict(rec, ref, temp_slack=5.0, dens_slack=0.02):
    """A stability filter only: drift, temperature blow-up, density blow-up."""
    if rec.get("rejected"):
        return False, "rejected by grompp"
    if rec.get("drift") is None:
        return False, rec.get("reason", "run failed")
    if abs(rec["drift"]) >= fm.DRIFT_LIMIT:
        return False, f"drift {rec['drift']:.3f}"
    if abs(rec["temperature"] - ref["temperature"]) > temp_slack:
        return False, (f"temperature {rec['temperature']:.1f} K vs "
                       f"{ref['temperature']:.1f} K")
    if abs(rec["density"] - ref["density"]) / ref["density"] > dens_slack:
        return False, "density"
    return True, ""


def ref_metrics(rec):
    return {"temperature": rec["temperature"], "density": rec["density"],
            "rdf_g": rec.get("rdf_g")}


def run_reference(gmx, conf, base_top, ps, strategy_dir, sub):
    ref_dir = os.path.join(strategy_dir, "reference_" + sub)
    os.makedirs(ref_dir, exist_ok=True)
    rec = fm.evaluate(gmx, dict(fm.REF_SPEC), conf, base_top, ps, ref_dir, True)
    if rec.get("rejected") or rec.get("drift") is None:
        sys.exit(f"{sub} reference failed: {rec.get('reason', rec)}")
    rec["pass"] = True
    rec["rdf_dev"] = 0.0
    return rec


def run_trials(state, specs, stage, ps, strategy_dir, gmx, conf, base_top,
               ref, certify, state_path):
    for spec in specs:
        label = fm.setting_label(spec)
        rec = {"spec": spec, "label": label, "stage": stage}
        key = stage + ":" + label
        if key in state["done"]:
            continue
        workdir = os.path.join(strategy_dir, label)
        state["current"] = {"label": label, "stage": stage,
                            "workdir": workdir, "started": time.time()}
        save(state, state_path)
        t0 = time.time()
        try:
            rec = fm.evaluate(gmx, spec, conf, base_top, ps, strategy_dir, True)
        except Exception as exc:  # noqa: BLE001 - record, never die mid-search
            rec = {"spec": spec, "label": label, "rejected": False,
                   "drift": None, "reason": f"exception: {exc}"}
        rec["stage"] = stage
        rec["wall"] = time.time() - t0
        if stage == "screen":
            ok, why = screen_verdict(rec, ref)
            rec["screen_pass"] = ok
            rec["screen_reason"] = why
            state["screen_times"].append(rec["wall"])
        else:
            fm.validate(rec, ref_metrics(ref))
            rec["certified"] = rec.get("pass", False)
        rec["key"] = key
        state["trials"].append(rec)
        state["done"][key] = True
        if stage == "screen":
            state["done_screen"] += 1
        else:
            state["done_certify"] += 1
        save(state, state_path)
        status = rec.get("screen_reason") if stage == "screen" else \
            ("PASS" if rec.get("pass") else "; ".join(rec.get("reasons", [])))
        print(f"  [{stage}] {label} -> {status or 'pass'}", flush=True)
    state["current"] = None
    save(state, state_path)


def save(state, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=2, default=str)
    os.replace(tmp, path)


def best_certified(state):
    best = None
    for rec in state["trials"]:
        if rec.get("stage") == "certify" and rec.get("certified"):
            if best is None or (rec.get("ns_per_day") or 0) > \
                    (best.get("ns_per_day") or 0):
                best = rec
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gmx", default=fm.GMX_DEFAULT)
    ap.add_argument("--ntomp", type=int, default=4)
    ap.add_argument("--screen-ps", type=float, default=20.0)
    ap.add_argument("--ref-ps", type=float, default=200.0)
    ap.add_argument("--certify", type=int, default=24)
    ap.add_argument("--dts", default="0.004,0.005,0.006,0.007,0.008,0.009,0.010")
    ap.add_argument("--hmr-factors", default="0,2.5,3,3.5,4,4.5",
                    help="0 means repartition off; above ~4.5 a methyl "
                         "carbon goes negative")
    ap.add_argument("--constraints", default="h-bonds,all-bonds,all-angles")
    ap.add_argument("--mass-scales", default="1.0",
                    help="global solute mass multiplier; not the same physics if >1")
    args = ap.parse_args()
    args.constraints = [c for c in args.constraints.split(",") if c]

    conf = os.path.abspath(args.gro)
    out = os.path.abspath(args.out)
    work = os.path.join(out, "runs")
    os.makedirs(os.path.join(work, "inputs"), exist_ok=True)
    base_conf = os.path.join(work, "inputs", os.path.basename(conf))
    shutil.copy(conf, base_conf)
    base_top = os.path.join(work, "inputs", os.path.basename(args.top))
    shutil.copy(os.path.abspath(args.top), base_top)
    for name, src in fm.local_includes(os.path.abspath(args.top)):
        shutil.copy(src, os.path.join(work, "inputs", name))

    specs = list(candidates(args, fm.is_protein(base_top)))
    state_path = os.path.join(out, "search.json")
    if os.path.exists(state_path):
        state = json.load(open(state_path))
    else:
        state = {"system": conf, "atoms": fm.read_gro_natoms(base_conf),
                 "space": args.__dict__, "trials": [], "done": {},
                 "done_screen": 0, "done_certify": 0, "screen_times": [],
                 "total_screen": len(specs), "total_certify": 0,
                 "current": None, "best": None}
    if not os.path.exists(state_path):
        save(state, state_path)

    gmx = fm.Gmx(args.gmx, work, args.ntomp)
    print(f"search: {state['atoms']} atoms, {len(specs)} screen candidates, "
          f"certify top {args.certify}")

    screen_ref = run_reference(gmx, base_conf, base_top, args.screen_ps, work,
                               "screen")
    state["screen_ref"] = screen_ref
    save(state, state_path)

    run_trials(state, specs, "screen", args.screen_ps, work, gmx, base_conf,
               base_top, screen_ref, False, state_path)

    passing = [r for r in state["trials"]
               if r.get("stage") == "screen" and r.get("screen_pass")]
    passing.sort(key=lambda r: r.get("ns_per_day") or 0.0, reverse=True)
    finalists = passing[:args.certify]
    state["total_certify"] = len(finalists)
    save(state, state_path)
    print(f"screen complete: {len(passing)} passed, certifying {len(finalists)}")

    cert_ref = run_reference(gmx, base_conf, base_top, args.ref_ps, work,
                             "cert")
    state["cert_ref"] = cert_ref
    save(state, state_path)
    run_trials(state, [r["spec"] for r in finalists], "certify", args.ref_ps,
               work, gmx, base_conf, base_top, cert_ref, True, state_path)

    state["best"] = best_certified(state)
    state["finished"] = time.time()
    save(state, state_path)
    if state["best"]:
        b = state["best"]
        speed = b["ns_per_day"] / cert_ref["ns_per_day"]
        print(f"\nBEST {b['label']}  {speed:.2f}x "
              f"({b['ns_per_day']:.1f} vs {cert_ref['ns_per_day']:.1f} ns/day)")
    else:
        print("\nno certified setting passed")


if __name__ == "__main__":
    main()
