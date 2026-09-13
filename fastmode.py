#!/usr/bin/env python3
"""fastmode: find the fastest safe GROMACS setting for a system, and prove it.

A production run usually uses dt = 2 fs and no virtual sites.  A correct setup
reaches 4 to 5 fs, which is a 2x speedup.  The settings are fiddly and a bad
choice fails silently, so this tool measures the claim on the user's own system
instead of asserting it.

The tool runs a 2 fs reference, then sweeps the timestep, hydrogen mass
repartition, multiple time stepping, the neighbour-list interval and the Verlet
buffer tolerance.  Each candidate is validated against the reference on four
checks.  It reports the fastest setting that passes every check, and it reports
every failure with its number.

Hydrogen mass repartition is a mass transform on the topology, so it is done in
Python here.  Hydrogen virtual sites need a pdb2gmx rebuild and live in a
separate stage (see stage2.py).  Water is skipped for repartition because SETTLE
already makes it rigid.

Style follows phase2/localdu.py: plain Python 3, numpy and the standard library.
"""

import argparse
import copy
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
GMX_DEFAULT = os.path.join(PROJECT, "build", "bin", "gmx")
GMXLIB_DEFAULT = os.path.join(PROJECT, "gromacs", "share", "top")

# Mass repartition factor.  Standard HMR value.
HMR_FACTOR = 4.0
# Total-mass conservation tolerance for the repartition transform.
MASS_TOL = 1e-6

# A run is only timed when the load from OTHER work is at or below this.  The
# sweep's own OpenMP threads are subtracted first: on an idle machine the
# 1-minute load average sits near our own thread count, so the literal reading
# would mark every run after the first as untimed.
LOAD_LIMIT = 2.0
# Steps thrown away before timing, because GROMACS tunes PME during them.
WARMUP_STEPS = 2000

# Validation limits, from the handoff.
DRIFT_LIMIT = 0.02          # kJ/mol/ps per atom
TEMP_LIMIT = 1.0            # K from the reference mean
RDF_LIMIT = 0.02            # max absolute deviation
DENSITY_LIMIT = 0.005       # 0.5 percent, relative

# Molecule types treated as water and therefore skipped by repartition.
WATER_NAMES = {"SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5", "SPC", "SPCE"}


def external_load(ntomp):
    """The 1-minute load average from work other than our own run."""
    load = os.getloadavg()[0]
    cores = os.cpu_count() or ntomp
    return max(0.0, load - min(ntomp, cores))


def machine_quiet(ntomp):
    """True when other work on the machine leaves it quiet enough to time."""
    return external_load(ntomp) <= LOAD_LIMIT


# --------------------------------------------------------------- input systems

def read_gro_natoms(path):
    with open(path) as fh:
        fh.readline()
        return int(fh.readline())


def parse_topology(text):
    """Return the inline moleculetypes, their atoms and their bonds.

    Only `[ atoms ]` and `[ bonds ]` blocks written directly in the file are
    returned.  Water and ion moleculetypes usually live in an included `.itp`,
    so they never appear here, which is what we want: repartition skips water.
    """
    lines = text.splitlines()
    mols = {}
    cur = None
    section = None
    for i, raw in enumerate(lines):
        body = raw.split(";", 1)[0]
        s = body.strip()
        if not s:
            continue
        if s.startswith("["):
            section = s.strip("[] ").lower().split()[0]
            continue
        if section == "moleculetype":
            name = s.split()[0]
            cur = name
            mols.setdefault(name, {"atoms": {}, "bonds": [], "order": []})
            continue
        if cur is None:
            continue
        toks = body.split()
        if section == "atoms" and len(toks) >= 6:
            nr = int(toks[0])
            name = toks[4]
            try:
                mass = float(toks[7]) if len(toks) >= 8 else 0.0
            except ValueError:
                mass = 0.0
            # character span of the mass token, so the line can be edited in place
            spans = [(m.start(), m.end()) for m in re.finditer(r"\S+", body)]
            span = spans[7] if len(spans) >= 8 else None
            mols[cur]["atoms"][nr] = {"name": name, "mass": mass,
                                      "line": i, "span": span}
            mols[cur]["order"].append(nr)
        elif section == "bonds" and len(toks) >= 2:
            mols[cur]["bonds"].append((int(toks[0]), int(toks[1])))
    return mols


def is_water(name):
    return name.upper() in WATER_NAMES


def bonded_heavy(atoms, bonds, h_nr):
    """The non-hydrogen partner of hydrogen atom `h_nr`, or None."""
    for a, b in bonds:
        other = b if a == h_nr else (a if b == h_nr else None)
        if other is None:
            continue
        other_name = atoms.get(other, {}).get("name", "")
        if not other_name.upper().startswith("H"):
            return other
    return None


def hmr_transform(top_text, out_path, factor=HMR_FACTOR, scale=1.0):
    """Repartition solute hydrogen mass, then scale every solute atom.

    Each hydrogen mass becomes `factor` times its own value, and the added mass
    is debited from its bonded heavy atom, so repartition alone conserves total
    mass.  `scale` then multiplies every solute atom, which is a pure frequency
    shift and is recorded separately because it also slows the dynamics.

    Writes the result to out_path, then re-parses that file and asserts that
    total solute mass equals `scale` times the original to MASS_TOL.  Returns
    (path, stats).
    """
    mols = parse_topology(top_text)
    lines = top_text.splitlines(keepends=True)

    work = {}
    for name, mol in mols.items():
        if is_water(name):
            continue
        work[name] = copy.deepcopy(mol)

    n_h = 0
    mass_before = 0.0
    for name, mol in work.items():
        atoms = mol["atoms"]
        bonds = mol["bonds"]
        mass_before += sum(a["mass"] for a in atoms.values())
        edits = {}
        for nr in mol["order"]:
            a = atoms[nr]
            if not a["name"].upper().startswith("H"):
                continue
            partner = bonded_heavy(atoms, bonds, nr)
            if partner is None or partner not in atoms:
                continue
            if atoms[partner]["mass"] <= 0.0:
                continue
            added = a["mass"] * (factor - 1.0)
            edits[nr] = a["mass"] * factor
            # a heavy atom can carry several hydrogens, so debit cumulatively
            edits[partner] = edits.get(partner, atoms[partner]["mass"]) - added
            n_h += 1
        if scale != 1.0:
            for nr in mol["order"]:
                edits[nr] = edits.get(nr, atoms[nr]["mass"]) * scale
        for nr, newmass in edits.items():
            a = atoms[nr]
            if a["span"] is None:
                sys.exit(f"cannot repartition {name} atom {nr}: no mass column")
            lines[a["line"]] = (lines[a["line"]][:a["span"][0]]
                                + f"{newmass:.4f}"
                                + lines[a["line"]][a["span"][1]:])

    if n_h == 0:
        return None, {"hydrogens": 0, "rel_mass_change": 0.0,
                      "factor": factor, "scale": scale}

    with open(out_path, "w") as fh:
        fh.writelines(lines)

    # Re-parse the file that was actually written, so the conservation check
    # covers the artifact and not just the arithmetic that produced it.
    after = 0.0
    for name, mol in parse_topology(open(out_path).read()).items():
        if is_water(name):
            continue
        after += sum(a["mass"] for a in mol["atoms"].values())
    expected = mass_before * scale
    rel = abs(after - expected) / expected
    if rel > MASS_TOL:
        sys.exit(f"repartition broke mass conservation: rel change {rel:.3e}")
    return out_path, {"hydrogens": n_h, "rel_mass_change": rel,
                      "factor": factor, "scale": scale,
                      "mass_before": mass_before, "mass_after": after}



def local_includes(top_path):
    """Bare `#include "name.itp"` files that sit beside the topology."""
    base = os.path.dirname(os.path.abspath(top_path))
    found = []
    with open(top_path) as fh:
        for line in fh:
            m = re.match(r'\s*#include\s+"([^"/]+)"', line)
            if not m:
                continue
            p = os.path.join(base, m.group(1))
            if os.path.exists(p):
                found.append((m.group(1), p))
    return found


# --------------------------------------------------------------------- mdp I/O

def output_interval(nprod, levels, want):
    """Round `want` up to a multiple of `levels`, never below `levels`."""
    return levels * max(1, math.ceil(want / levels))


# RDF statistics need dense coordinates.  Saving every 0.5 ps leaves a max
# deviation noise floor near 0.02 between independent runs; every 0.05 ps drops
# it to under 0.01.  Keep the interval physical, not step-count based.
RDF_SAVE_PS = 0.05


def build_mdp(spec, nsteps, ref_ps):
    """An mdp string for one candidate.  MTS dictates output intervals."""
    f = spec["mts"] if spec["mts"] else 1
    nstcalc = f if spec["mts"] else 100
    nstenergy = output_interval(nsteps, nstcalc, nsteps // 100 or 1)
    nstlog = output_interval(nsteps, nstcalc, nsteps // 100 or 1)
    save = max(1, int(round(RDF_SAVE_PS / spec["dt"])))
    if spec["mts"]:
        save = max(f, int(round(save / f)) * f)
    nstxout = save
    lines = [
        "integrator      = md",
        f"dt              = {spec['dt']}",
        f"nsteps          = {nsteps}",
        f"nstlist         = {spec['nstlist']}",
        "cutoff-scheme   = Verlet",
        f"verlet-buffer-tolerance = {spec['tol']}",
        "coulombtype     = PME",
        "rcoulomb        = 1.0",
        "rvdw            = 1.0",
        "fourierspacing  = 0.12",
        "pme-order       = 4",
        "tcoupl          = v-rescale",
        "tc-grps         = System",
        # tau-t 0.1 matches bench/base.mdp, the setup behind the project's
        # measured numbers.  At 0.5, grompp treats the buffer-tolerance 0.05
        # warning as fatal, which would make that sweep point untestable.
        "tau-t           = 0.1",
        "ref-t           = 300",
        "pcoupl          = C-rescale",
        "pcoupltype      = isotropic",
        "tau-p           = 2.0",
        "ref-p           = 1.0",
        "compressibility = 4.5e-5",
        "constraints     = " + spec.get("constraints", "h-bonds"),
        "constraint-algorithm = lincs",
        f"nstcalcenergy   = {nstcalc}",
        f"nstenergy       = {nstenergy}",
        f"nstlog          = {nstlog}",
        "nstxout         = 0",
        "nstvout         = 0",
        f"nstxout-compressed = {nstxout}",
        "compressed-x-grps = System",
    ]
    if spec.get("seed") is None:
        lines.append("gen-vel         = no")
    else:
        # A fresh velocity draw makes repeated runs independent samples, which
        # is what `selfcheck` needs to measure the verifier's noise floor.
        lines.append("gen-vel         = yes")
        lines.append("gen-temp        = 300")
        lines.append(f"gen-seed        = {int(spec['seed'])}")
    if spec["mts"]:
        lines += [
            "mts = yes",
            "mts-levels = 2",
            "mts-level2-forces = longrange-nonbonded",
            f"mts-level2-factor = {spec['mts']}",
        ]
    return "\n".join(lines) + "\n", nstxout


def setting_label(spec):
    hmr = "on" if spec["hmr"] else "off"
    mts = "off" if not spec["mts"] else f"pme{spec['mts']}"
    label = (f"dt{int(round(spec['dt'] * 1000))}fs_hmr{hmr}"
             f"_mts{mts}"
             f"_nstlist{spec['nstlist']}_tol{spec['tol']}")
    if spec.get("hmr") and spec.get("hmr_factor", HMR_FACTOR) != HMR_FACTOR:
        label += f"_f{spec['hmr_factor']:g}"
    if spec.get("seed") is not None:
        label += f"_seed{int(spec['seed'])}"
    cons = spec.get("constraints", "h-bonds")
    if cons != "h-bonds":
        label += "_" + cons.replace("-", "")
    if spec.get("mass_scale", 1.0) != 1.0:
        label += f"_s{spec['mass_scale']:g}"
    return label


# ------------------------------------------------------------------- gromacs

class Gmx:
    def __init__(self, binary, workdir, ntomp):
        self.binary = os.path.abspath(binary)
        self.workdir = workdir
        self.ntomp = ntomp
        self.env = dict(os.environ, GMXLIB=os.environ.get("GMXLIB", GMXLIB_DEFAULT))

    def _path(self, name):
        return os.path.join(self.workdir, name)

    def grompp(self, mdp_path, conf, top, tag):
        out = subprocess.run(
            [self.binary, "grompp", "-f", mdp_path, "-c", conf, "-p", top,
             "-o", self._path(tag + ".tpr")],
            cwd=self.workdir, env=self.env, text=True, capture_output=True)
        return out

    def mdrun(self, tag):
        out = subprocess.run(
            [self.binary, "mdrun", "-s", self._path(tag + ".tpr"),
             "-deffnm", self._path(tag), "-ntmpi", "1", "-ntomp", str(self.ntomp),
             "-nb", "cpu", "-pin", "off", "-resetstep", str(WARMUP_STEPS)],
            cwd=self.workdir, env=self.env, text=True, capture_output=True)
        return out

    def energy(self, tag, warmup_ps):
        out = subprocess.run(
            [self.binary, "energy", "-f", self._path(tag + ".edr"),
             "-b", f"{warmup_ps:.4f}", "-o", self._path(tag + "_ener.xvg")],
            cwd=self.workdir, env=self.env, text=True, capture_output=True,
            input="Temperature\nDensity\n\n")
        return out

    def rdf(self, tag, warmup_ps, rmax=1.0, binwidth=0.01):
        out = subprocess.run(
            [self.binary, "rdf", "-f", self._path(tag + ".xtc"),
             "-s", self._path(tag + ".tpr"),
             "-ref", "name OW", "-sel", "name OW",
             "-rmax", str(rmax), "-bin", str(binwidth),
             "-b", f"{warmup_ps:.4f}", "-o", self._path(tag + "_rdf.xvg")],
            cwd=self.workdir, env=self.env, text=True, capture_output=True)
        return out


def performance_ns_per_day(log_path):
    with open(log_path, errors="ignore") as fh:
        for line in fh:
            if line.strip().startswith("Performance:"):
                return float(line.split()[1])
    return None


def conserved_drift(log_path):
    with open(log_path, errors="ignore") as fh:
        text = fh.read()
    m = re.search(r"Conserved energy drift:\s*([-+0-9.eE]+)", text)
    return float(m.group(1)) if m else None


def read_xvg(path):
    """Return a list of numpy columns from an xvg file."""
    cols, legends = [], []
    rows = {i: [] for i in range(16)}
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = re.match(r'\s*@\s+s\d+\s+legend\s+"(.*)"', line)
            if m:
                legends.append(m.group(1))
            elif line[:1] not in ("#", "@") and line.strip():
                for i, v in enumerate(line.split()):
                    rows.setdefault(i, []).append(float(v))
    indexes = sorted(rows)
    for i in indexes:
        cols.append(np.array(rows[i]))
    return cols, legends


# ------------------------------------------------------------------ evaluation

def derive_top(spec, base_top, workdir):
    """Topology path for a candidate: the repartition copy when it is on."""
    if not spec["hmr"] and spec.get("mass_scale", 1.0) == 1.0:
        return base_top
    out_top = os.path.join(workdir, "top_hmr.top")
    top_text = open(base_top).read()
    path, stats = hmr_transform(top_text, out_top,
                                spec.get("hmr_factor", HMR_FACTOR),
                                spec.get("mass_scale", 1.0))
    if path is None:
        return base_top
    for name, src in local_includes(base_top):
        shutil.copy(src, os.path.join(workdir, name))
    return out_top


def evaluate(gmx, spec, conf, base_top, ref_ps, strategy_dir, timed):
    """Run one candidate and its checks.  Returns a record dict."""
    workdir = os.path.join(strategy_dir, setting_label(spec))
    os.makedirs(workdir, exist_ok=True)
    tag = "run"
    # dt and ref_ps are both in ps, so steps = ps / (ps per step)
    nprod = int(round(ref_ps / spec["dt"]))
    nsteps = WARMUP_STEPS + nprod
    warmup_ps = WARMUP_STEPS * spec["dt"] / 1000.0

    mdp_text, _ = build_mdp(spec, nsteps, ref_ps)
    mdp_path = os.path.join(workdir, "run.mdp")
    previous_mdp = open(mdp_path).read() if os.path.exists(mdp_path) else None
    if previous_mdp != mdp_text:
        with open(mdp_path, "w") as fh:
            fh.write(mdp_text)

    top = derive_top(spec, base_top, workdir)
    g = Gmx(gmx.binary, workdir, gmx.ntomp)

    rec = {"spec": dict(spec), "label": setting_label(spec),
           "workdir": workdir, "timed": timed,
           "rejected": False, "pass": False, "nsteps": nsteps,
           "reasons": []}

    log = os.path.join(workdir, tag + ".log")
    complete = (previous_mdp == mdp_text
                and os.path.exists(log)
                and os.path.exists(os.path.join(workdir, tag + ".edr"))
                and os.path.exists(os.path.join(workdir, tag + ".xtc"))
                and "Finished mdrun" in open(log, errors="ignore").read())
    if not complete:
        gp = g.grompp(mdp_path, os.path.abspath(conf), os.path.abspath(top), tag)
        if gp.returncode != 0:
            rec["rejected"] = True
            rec["reason"] = last_error(gp.stderr + "\n" + gp.stdout)
            return rec
        md = g.mdrun(tag)
        finished = os.path.exists(log) and "Finished mdrun" in \
            open(log, errors="ignore").read()
        if md.returncode != 0 or not finished:
            rec["reason"] = "mdrun failed: " + last_error(md.stderr + "\n"
                                                          + md.stdout)
            return rec
    rec["ns_per_day"] = performance_ns_per_day(log)
    rec["drift"] = conserved_drift(log)
    if rec["ns_per_day"] is None or rec["drift"] is None:
        rec["reason"] = "log has no Performance or Conserved energy drift line"
        return rec

    en = g.energy(tag, warmup_ps)
    if en.returncode != 0:
        rec["reasons"].append("energy extraction failed")
        rec["reason"] = last_error(en.stderr)
        return rec
    cols, legends = read_xvg(os.path.join(workdir, tag + "_ener.xvg"))
    rec["temperature"] = float(np.mean(cols[1]))
    rec["density"] = float(np.mean(cols[2]))

    rd = g.rdf(tag, warmup_ps)
    if rd.returncode != 0:
        rec["reasons"].append("rdf failed")
        rec["reason"] = last_error(rd.stderr)
        return rec
    rcols, _ = read_xvg(os.path.join(workdir, tag + "_rdf.xvg"))
    rec["rdf_r"] = rcols[0].tolist()
    rec["rdf_g"] = rcols[1].tolist()
    return rec


def last_error(stderr):
    """The message after the last ERROR or Fatal error line, not a separator."""
    def real(line):
        s = line.strip()
        return bool(s) and not set(s) <= set("-=_ ")
    lines = stderr.splitlines()
    for i in range(len(lines) - 1, -1, -1):
        if "Fatal error" in lines[i] or re.match(r"\s*ERROR\b", lines[i]):
            tail = [l.strip() for l in lines[i + 1:i + 4] if real(l)]
            msg = " ".join(tail) if tail else lines[i].strip()
            return msg[:300]
    for line in reversed(lines):
        if real(line):
            return line.strip()[:200]
    return "grompp failed"


def validate(rec, ref):
    """Apply the four checks in place.  Returns True if every check passes."""
    reasons = list(rec.get("reasons", []))
    if rec.get("rejected") or "drift" not in rec or rec.get("drift") is None:
        rec["pass"] = False
        return False
    if abs(rec["drift"]) >= DRIFT_LIMIT:
        reasons.append(f"drift {rec['drift']:.3f} exceeds limit {DRIFT_LIMIT}")
    if abs(rec["temperature"] - ref["temperature"]) > TEMP_LIMIT:
        reasons.append(f"temperature {rec['temperature']:.2f} K vs "
                       f"{ref['temperature']:.2f} K exceeds {TEMP_LIMIT} K")
    if abs(rec["density"] - ref["density"]) / ref["density"] > DENSITY_LIMIT:
        pct = 100 * abs(rec["density"] - ref["density"]) / ref["density"]
        reasons.append(f"density {pct:.2f}% from reference exceeds "
                       f"{100 * DENSITY_LIMIT:.1f}%")
    if "rdf_g" in rec and ref.get("rdf_g") is not None:
        a = np.array(rec["rdf_g"])
        b = np.array(ref["rdf_g"])
        n = min(a.size, b.size)
        dev = float(np.max(np.abs(a[:n] - b[:n])))
        rec["rdf_dev"] = dev
        if dev >= RDF_LIMIT:
            reasons.append(f"rdf max deviation {dev:.4f} exceeds {RDF_LIMIT}")
    else:
        rec["rdf_dev"] = None
        reasons.append("rdf check did not run")
    rec["reasons"] = reasons
    rec["pass"] = len(reasons) == 0
    return rec["pass"]


# ----------------------------------------------------------------- the sweeps

REF_SPEC = {"dt": 0.002, "hmr": False, "mts": 0, "nstlist": 10, "tol": 0.005}
DT_VALUES = [0.002, 0.003, 0.004, 0.005]
MTS_VALUES = [0, 2, 3]
NSTLIST_VALUES = [10, 25, 50, 100]
TOL_VALUES = [0.005, 0.05]


def vary(spec, **kw):
    out = dict(spec)
    out.update(kw)
    return out


def is_protein(base_top):
    mols = parse_topology(open(base_top).read())
    return any(not is_water(name) for name in mols)


def grid_specs(base, protein):
    dts = DT_VALUES
    hmrs = [False, True] if protein else [False]
    for dt, hmr, mts, nst, tol in itertools.product(
            dts, hmrs, MTS_VALUES, NSTLIST_VALUES, TOL_VALUES):
        yield vary(base, dt=dt, hmr=hmr, mts=mts, nstlist=nst, tol=tol)


# --------------------------------------------------------------------- report

def fmt_run(rec):
    if rec.get("rejected"):
        return f"  {rec['label']:<52} rejected by grompp: {rec.get('reason','')}"
    if "drift" not in rec or rec.get("drift") is None:
        return f"  {rec['label']:<52} run failed: {rec.get('reason','')}"
    speed = (f"{rec['ns_per_day']:8.1f} ns/day" if rec.get("ns_per_day")
             else "   untimed")
    verdict = "PASS" if rec["pass"] else "FAIL"
    rdf = rec.get("rdf_dev")
    rdf_s = f"{rdf:.4f}" if rdf is not None else "  n/a "
    line = (f"  {rec['label']:<52} drift {rec['drift']:9.4f}  "
            f"T {rec['temperature']:7.2f}  rdfdev {rdf_s}  "
            f"rho {rec['density']:8.4f}  {speed}  {verdict}")
    if rec["reasons"]:
        line += "\n" + "\n".join(f"        -> {r}" for r in rec["reasons"])
    return line


def write_report(path, args, natoms, ref, runs, best):
    lines = []
    lines.append("fastmode audit report")
    lines.append("=====================")
    lines.append(f"system       : {args.gro}")
    lines.append(f"topology     : {args.top}")
    lines.append(f"atoms        : {natoms}")
    lines.append(f"reference    : dt 2 fs, hmr off, mts off, nstlist 10, tol 0.005")
    lines.append(f"production   : {args.ref_ps:.0f} ps per run, {WARMUP_STEPS} warmup steps")
    lines.append(f"threads      : {args.ntomp} OpenMP, 1 MPI")
    lines.append(f"strategy     : {args.strategy}")
    any_untimed = any(not r.get("rejected") and not r.get("timed") for r in runs)
    if any_untimed:
        lines.append("timing       : load average above 2.0 for at least one run, so")
        lines.append("               every performance number below is untimed and the")
        lines.append("               speedup is withheld.  Re-run on a quiet machine.")
    lines.append("")
    lines.append("Reference, measured:")
    lines.append(f"  drift {ref['drift']:.4f} kJ/mol/ps/atom   "
                 f"temperature {ref['temperature']:.2f} K   "
                 f"density {ref['density']:.4f}")
    lines.append("")
    lines.append("Candidates:")
    for rec in runs:
        lines.append(fmt_run(rec))
    lines.append("")
    lines.append("Result")
    lines.append("------")
    if best is None:
        lines.append("  no setting passed every check.  The reference settings are the")
        lines.append("  only safe choice on this system.  See the failures above.")
    else:
        lines.append(f"  fastest passing setting: {best['label']}")
        if best.get("ns_per_day") and ref.get("ns_per_day") and not any_untimed:
            lines.append(f"  speedup over the 2 fs reference : "
                         f"{best['ns_per_day'] / ref['ns_per_day']:.2f}x")
        else:
            lines.append("  speedup over the 2 fs reference : withheld (untimed)")
        lines.append(f"  conserved drift  {best['drift']:.4f} kJ/mol/ps/atom "
                     f"(limit {DRIFT_LIMIT})")
        lines.append(f"  temperature      {best['temperature']:.2f} K vs "
                     f"{ref['temperature']:.2f} K")
        rdf = best.get("rdf_dev")
        if rdf is not None:
            lines.append(f"  rdf max deviation {rdf:.4f} (limit {RDF_LIMIT})")
        pct = 100 * abs(best["density"] - ref["density"]) / ref["density"]
        lines.append(f"  density          {pct:.2f}% from reference "
                     f"(limit {100 * DENSITY_LIMIT:.1f}%)")
    lines.append("")
    lines.append("Assumptions")
    lines.append("-----------")
    lines.append("  - Timestep changes scale the production length, so every run")
    lines.append("    covers the same physical time as the reference.")
    lines.append("  - Hydrogen repartition is skipped for pure water: SETTLE already")
    lines.append("    makes the water rigid, so repartition buys nothing there.")
    lines.append("  - The thermostat coupling is tau-t 0.1, matching bench/base.mdp.")
    lines.append("    At 0.5, grompp rejects verlet-buffer-tolerance 0.05 as fatal.")
    lines.append("  - The reference itself carries finite-run statistical error.")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return "\n".join(lines)


# ----------------------------------------------------------------------- main

def audit(args):
    conf = os.path.abspath(args.gro)
    base_top = os.path.abspath(args.top)
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "runs")
    os.makedirs(work, exist_ok=True)
    os.makedirs(os.path.join(work, "inputs"), exist_ok=True)
    shutil.copy(conf, os.path.join(work, "inputs", os.path.basename(conf)))
    top_copy = os.path.join(work, "inputs", os.path.basename(base_top))
    shutil.copy(base_top, top_copy)
    for name, src in local_includes(base_top):
        shutil.copy(src, os.path.join(work, "inputs", name))
    base_top = top_copy

    natoms = read_gro_natoms(conf)
    protein = is_protein(base_top)
    gmx = Gmx(args.gmx, work, args.ntomp)
    load = external_load(args.ntomp)
    timed = load <= LOAD_LIMIT
    print(f"fastmode: {natoms} atoms, protein={protein}, "
          f"other load {load:.2f}, timed={timed}")

    ref_spec = dict(REF_SPEC)
    ref_dir = os.path.join(work, "reference")
    os.makedirs(ref_dir, exist_ok=True)
    print("running 2 fs reference ...")
    ref = evaluate(gmx, ref_spec, conf, base_top, args.ref_ps, ref_dir, timed)
    if ref.get("rejected") or ref.get("drift") is None:
        sys.exit(f"reference run failed: {ref.get('reason', ref)}")

    ref["pass"] = True
    ref["rdf_dev"] = 0.0
    runs = [ref]
    ref_metrics = {"temperature": ref["temperature"], "density": ref["density"],
                   "rdf_g": ref["rdf_g"]}
    seen = {setting_label(ref_spec)}

    def run_one(spec):
        label = setting_label(spec)
        if label in seen:
            return None
        seen.add(label)
        load = external_load(args.ntomp)
        run_timed = load <= LOAD_LIMIT
        print(f"  {label} ... other load {load:.2f}", flush=True)
        rec = evaluate(gmx, spec, conf, base_top, args.ref_ps,
                       os.path.join(work, "runs"), run_timed)
        if rec.get("rejected"):
            print(f"    rejected: {rec.get('reason','')}")
        else:
            validate(rec, ref_metrics)
            print(f"    {'PASS' if rec['pass'] else 'fail'}"
                  + ("" if rec["pass"]
                     else f" ({'; '.join(rec.get('reasons', []))})"))
        runs.append(rec)
        return rec

    def fastest_passing(records):
        passing = [r for r in records if r.get("pass")]
        if not passing:
            return None
        return max(passing, key=lambda r: r.get("ns_per_day") or 0.0)

    if args.strategy == "grid":
        for spec in grid_specs(ref_spec, protein):
            run_one(spec)
    else:
        # Topology and timestep interact.  Hydrogen repartition does not speed
        # a run up at a fixed timestep; it enables a larger one.  So sweep the
        # full timestep ladder for each topology and keep the faster topology,
        # instead of choosing a topology at 2 fs where both are equal.
        branches = [False, True] if protein else [False]
        branch_best = []
        for hmr in branches:
            dts = DT_VALUES if hmr else DT_VALUES[1:]
            records = []
            for dt in dts:
                rec = run_one(vary(ref_spec, hmr=hmr, dt=dt))
                if rec is not None:
                    records.append(rec)
            winner = fastest_passing(records)
            if winner is not None:
                branch_best.append(winner)

        base = dict(ref_spec)
        chosen = fastest_passing(branch_best)
        if chosen is not None:
            base = dict(chosen["spec"])

        # Coordinate descent on the remaining axes from the chosen base.
        for param, values in [("mts", MTS_VALUES),
                              ("nstlist", NSTLIST_VALUES),
                              ("tol", TOL_VALUES)]:
            records = []
            for value in values:
                if value == base[param]:
                    continue
                rec = run_one(vary(base, **{param: value}))
                if rec is not None:
                    records.append(rec)
            better = fastest_passing(records)
            if better is not None:
                base = dict(better["spec"])

    # The reference must not compete with itself: its checks are exact by
    # construction.  It remains the fallback when nothing faster passes.
    passing = [r for r in runs if r.get("pass") and r is not ref]
    # rank by measured ns/day; this is only a selection key, never a claim
    best = max(passing, key=lambda r: r.get("ns_per_day") or 0.0) if passing else None

    report_path = os.path.join(out, "report.txt")
    text = write_report(report_path, args, natoms, ref, runs, best)
    with open(os.path.join(out, "results.json"), "w") as fh:
        json.dump({"system": conf, "atoms": natoms, "reference": ref,
                   "runs": runs, "best": best}, fh, indent=2, default=str)
    print()
    print(text)


def hmr_check(args):
    text = open(args.top).read()
    out = args.out or os.path.join(HERE, "top_hmr.top")
    path, stats = hmr_transform(text, out)
    print(f"hydrogens repartitioned : {stats['hydrogens']}")
    print(f"relative mass change    : {stats['rel_mass_change']:.3e} "
          f"(limit {MASS_TOL:.0e})")
    if path:
        print(f"wrote {path}")
    else:
        print("no solute hydrogens: repartition not applicable")


def selfcheck(args):
    """Measure the verifier's own noise floor from independent 2 fs samples.

    Two or more runs of the SAME 2 fs settings differ only in their initial
    velocity draw.  The spread between them is the smallest difference the
    four checks can resolve.  This command reports that spread.  It never
    moves a threshold.
    """
    conf = os.path.abspath(args.gro)
    base_top = os.path.abspath(args.top)
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "runs")
    os.makedirs(os.path.join(work, "inputs"), exist_ok=True)
    shutil.copy(conf, os.path.join(work, "inputs", os.path.basename(conf)))
    top_copy = os.path.join(work, "inputs", os.path.basename(base_top))
    shutil.copy(base_top, top_copy)
    for name, src in local_includes(base_top):
        shutil.copy(src, os.path.join(work, "inputs", name))
    base_top = top_copy

    natoms = read_gro_natoms(conf)
    gmx = Gmx(args.gmx, work, args.ntomp)
    load = external_load(args.ntomp)
    timed = load <= LOAD_LIMIT
    print(f"selfcheck: {natoms} atoms, {args.runs} samples of "
          f"{args.ref_ps:.0f} ps, other load {load:.2f}, timed={timed}")

    recs = []
    for i in range(args.runs):
        spec = dict(REF_SPEC)
        spec["seed"] = args.seed + i
        print(f"  sample {i + 1}/{args.runs} seed {spec['seed']} ...",
              flush=True)
        rec = evaluate(gmx, spec, conf, base_top, args.ref_ps,
                       os.path.join(work, "samples"), timed)
        if rec.get("rejected") or rec.get("drift") is None:
            sys.exit(f"sample failed: {rec.get('reason', rec)}")
        recs.append(rec)

    def pair_max(key, fn=lambda a, b: abs(a - b)):
        best = 0.0
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                best = max(best, fn(recs[i][key], recs[j][key]))
        return best

    t_floor = pair_max("temperature")
    rho_floor = 100.0 * pair_max("density", lambda a, b: abs(a - b) / min(a, b))
    rdf_floor = 0.0
    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            gi = np.asarray(recs[i]["rdf_g"])
            gj = np.asarray(recs[j]["rdf_g"])
            n = min(len(gi), len(gj))
            rdf_floor = max(rdf_floor, float(np.max(np.abs(gi[:n] - gj[:n]))))

    def verdict(value, limit):
        return "resolves" if value < limit else "AT OR ABOVE LIMIT"

    lines = []
    lines.append("fastmode selfcheck report")
    lines.append("=========================")
    lines.append(f"system       : {args.gro}")
    lines.append(f"topology     : {args.top}")
    lines.append(f"atoms        : {natoms}")
    lines.append(f"samples      : {args.runs} independent 2 fs runs, "
                 f"{args.ref_ps:.0f} ps each")
    lines.append(f"seeds        : {args.seed} .. {args.seed + args.runs - 1}")
    lines.append(f"threads      : {args.ntomp} OpenMP, 1 MPI")
    if not timed:
        lines.append("timing       : load average above 2.0; the selfcheck does "
                     "not use")
        lines.append("               timings, so this does not affect the result.")
    lines.append("")
    lines.append("Samples:")
    for rec in recs:
        lines.append(f"  seed {rec['spec']['seed']:<3} "
                     f"drift {rec['drift']:.4f}  "
                     f"T {rec['temperature']:7.2f}  "
                     f"rho {rec['density']:8.4f}")
    lines.append("")
    lines.append("Noise floor (largest difference between any two samples):")
    lines.append(f"  temperature      {t_floor:.3f} K   "
                 f"(check limit {TEMP_LIMIT})   {verdict(t_floor, TEMP_LIMIT)}")
    lines.append(f"  density          {rho_floor:.3f} %   "
                 f"(check limit {100 * DENSITY_LIMIT:.1f})   "
                 f"{verdict(rho_floor, 100.0 * DENSITY_LIMIT)}")
    lines.append(f"  rdf max dev      {rdf_floor:.4f}   "
                 f"(check limit {RDF_LIMIT})   {verdict(rdf_floor, RDF_LIMIT)}")
    lines.append("")
    lines.append("Reading")
    lines.append("-------")
    lines.append("  A candidate difference smaller than the floor is not resolved")
    lines.append("  by this verifier at this run length.  The floors fall as the")
    lines.append("  run gets longer and, for the RDF, as the output interval")
    lines.append(f"  gets denser (currently {RDF_SAVE_PS} ps).  The thresholds are")
    lines.append("  fixed; this report only says whether each check can see a")
    lines.append("  difference of one floor.")
    path = os.path.join(out, "selfcheck.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print()
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("audit", help="sweep and validate settings for a system")
    a.add_argument("--gro", required=True)
    a.add_argument("--top", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--gmx", default=GMX_DEFAULT)
    a.add_argument("--ntomp", type=int, default=4)
    # 200 ps keeps the temperature noise floor near 0.5 K, below the 1 K check.
    # At 50 ps the floor is about 1.1 K and correct settings fail.
    a.add_argument("--ref-ps", type=float, default=200.0)
    a.add_argument("--strategy", choices=["staged", "grid"], default="staged")
    a.set_defaults(func=audit)

    h = sub.add_parser("hmr-check", help="prove the mass transform conserves mass")
    h.add_argument("--top", required=True)
    h.add_argument("--out", default=None)
    h.set_defaults(func=hmr_check)

    s = sub.add_parser("selfcheck",
                       help="measure the verifier's own noise floor")
    s.add_argument("--gro", required=True)
    s.add_argument("--top", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--gmx", default=GMX_DEFAULT)
    s.add_argument("--ntomp", type=int, default=4)
    s.add_argument("--ref-ps", type=float, default=200.0)
    s.add_argument("--runs", type=int, default=3)
    s.add_argument("--seed", type=int, default=1)
    s.set_defaults(func=selfcheck)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
