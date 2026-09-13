"""Idea 4 probe: is the MTS failure a correctable bias or a resonance?

MTS updates the long-range (PME) force only every `factor` steps.  If the
failure is a smooth force bias, a learned correction can remove it.  If it is a
resonance, the slow-force sampling interval excites a natural frequency and no
smooth force correction fixes it without changing the splitting.  Map drift
against (dt, factor) and look at the pattern.

Writes mtsres.txt.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import fastmode as fm

CONF = os.path.join(fm.PROJECT, "phase0/systems/npt_3.0.gro")
TOP = os.path.join(fm.PROJECT, "phase0/systems/topol_3.0.top")
OUT = os.path.join(HERE, "mtsres", "runs")
REF_PS = 20.0
DTS = [0.002, 0.003, 0.004, 0.005]
FACTORS = [0, 2, 3, 4, 5]
BASE = {"hmr": False, "nstlist": 10, "tol": 0.005, "constraints": "h-bonds"}


def main():
    os.makedirs(OUT, exist_ok=True)
    gmx = fm.Gmx(fm.GMX_DEFAULT, OUT, 4)
    timed = fm.machine_quiet(4)
    print(f"other load quiet: {timed}")
    rows = []
    for dt in DTS:
        for factor in FACTORS:
            spec = dict(BASE, dt=dt, mts=factor)
            rec = fm.evaluate(gmx, spec, CONF, TOP, REF_PS, OUT, timed)
            drift = rec.get("drift")
            temp = rec.get("temperature")
            note = ""
            if rec.get("rejected"):
                note = "grompp reject"
            elif drift is None:
                note = rec.get("reason", "failed")
            print(f"dt{dt*1000:g} f{factor}: drift={drift} T={temp} {note}")
            rows.append((dt, factor, drift, temp, rec.get("rejected"),
                         rec.get("reason", "")))
    write(rows)


def write(rows):
    path = os.path.join(HERE, "mtsres.txt")
    with open(path, "w") as fh:
        fh.write("MTS resonance map, water small, 20 ps production, 4 threads\n")
        fh.write("drift in kJ/mol/ps/atom; limit 0.02; note 'grompp reject' "
                 "means grompp refused the setting\n\n")
        fh.write(f"{'dt(fs)':>7} {'factor':>7} {'drift':>10} {'T(K)':>8} "
                 f"{'verdict':>10}  note\n")
        for dt, factor, drift, temp, rejected, reason in rows:
            if rejected:
                verdict, d = "rejected", "n/a"
            elif drift is None:
                verdict, d = "failed", "n/a"
            else:
                verdict = "PASS" if abs(drift) < fm.DRIFT_LIMIT else "FAIL"
                d = f"{drift:+.4f}"
            t = f"{temp:.2f}" if temp is not None else "n/a"
            fh.write(f"{dt*1000:7.0f} {factor:7d} {d:>10} {t:>8} "
                     f"{verdict:>10}  {reason[:40]}\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
