"""Local water-patch Metropolis probe (phase0/ml).

Scores a local proposal by the exact GROMACS potential energy of the full box.
A proposal translates and rotates whole water molecules (rigid, so SETTLE stays
satisfied). The move is a symmetric random walk, so the Metropolis acceptance is
min(1, exp(-dU / kT)). The baseline is the phase2/RESULTS.md table: acceptance
falls as exp(-1.6 n_waters) for a random patch move.

Milestone 1: reproduce that baseline on the water-small box and record the number
a learned conditional proposal must beat.

Writes patch.txt. Scratch in patch/.
"""

import argparse
import math
import os
import re
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
GMX = os.path.join(PROJECT, "build", "bin", "gmx")
GMXLIB = os.path.join(PROJECT, "gromacs", "share", "top")
KT = 8.314462618e-3 * 300.0
WORK = os.path.join(HERE, "patch")

MDP = """\
integrator              = md
nsteps                  = 1
dt                      = 0.001
cutoff-scheme           = Verlet
verlet-buffer-tolerance = 0.005
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
fourierspacing          = 0.12
pme-order               = 4
constraints             = h-bonds
constraint-algorithm    = lincs
tcoupl                  = no
pcoupl                  = no
nstlist                 = 10
nstcalcenergy           = 1
nstenergy               = 1
nstxout                 = 0
nstvout                 = 0
nstlog                  = 1
"""


def read_gro(path):
    with open(path) as handle:
        lines = handle.readlines()
    natoms = int(lines[1])
    names, coords = [], []
    for line in lines[2:2 + natoms]:
        names.append(line[10:15].strip())
        coords.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
    box = [float(x) for x in lines[2 + natoms].split()[:3]]
    return names, np.array(coords, float), np.array(box, float)


def write_gro_frame(handle, title, names, coords, box):
    n = len(names)
    handle.write(f"{title}\n")
    handle.write(f"{n:5d}\n")
    for i, (name, (x, y, z)) in enumerate(zip(names, coords), start=1):
        handle.write(f"{i % 100000:5d}{'WAT':>5}{name:>5}{i % 100000:5d}"
                     f"{x:8.3f}{y:8.3f}{z:8.3f}\n")
    handle.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


def water_index(names):
    groups, current = [], []
    for i, name in enumerate(names):
        if name == "OW":
            if current:
                groups.append(current)
            current = [i]
        else:
            current.append(i)
    if current:
        groups.append(current)
    return [g for g in groups if len(g) == 3]


def small_rotation(rng, sigma):
    vec = rng.normal(scale=sigma, size=3)
    theta = float(np.linalg.norm(vec))
    if theta < 1e-12:
        return np.eye(3)
    k = vec / theta
    cross = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(theta) * cross + (1 - math.cos(theta)) * (cross @ cross)


def min_image(delta, box):
    return delta - box * np.round(delta / box)


def make_proposal(coords, waters, oxygens, rng, sigma_t, sigma_r):
    out = coords.copy()
    patch = [w for w in waters if any(o == w[0] for o in oxygens)]
    for w in patch:
        oxygen = coords[w[0]]
        rot = small_rotation(rng, sigma_r)
        shift = rng.normal(scale=sigma_t, size=3)
        out[w] = (coords[w] - oxygen) @ rot.T + oxygen + shift
    return out, len(patch)


def run(cmd, cwd, stdin=None):
    return subprocess.run(cmd, cwd=cwd, env=dict(os.environ, GMXLIB=GMXLIB),
                          input=stdin, capture_output=True, text=True)


def prepare(system_gro, system_top):
    os.makedirs(WORK, exist_ok=True)
    mdp = os.path.join(WORK, "rerun.mdp")
    with open(mdp, "w") as handle:
        handle.write(MDP)
    top = os.path.join(WORK, "topol.top")
    if not os.path.exists(top):
        with open(system_top) as src, open(top, "w") as dst:
            dst.write(src.read())
    tpr = os.path.join(WORK, "rerun.tpr")
    if not os.path.exists(tpr):
        res = run([GMX, "grompp", "-f", "rerun.mdp", "-c", system_gro,
                   "-p", "topol.top", "-o", "rerun.tpr", "-po", "mdout.mdp"], WORK)
        if res.returncode != 0:
            raise SystemExit(res.stderr[-2000:])
    return tpr


def energies(frames_path):
    trr = os.path.join(WORK, "frames.trr")
    res = run([GMX, "trjconv", "-f", frames_path, "-s", "rerun.tpr",
               "-o", "frames.trr"], WORK, stdin="0\n")
    if res.returncode != 0:
        raise SystemExit(res.stderr[-2000:])
    res = run([GMX, "mdrun", "-s", "rerun.tpr", "-rerun", "frames.trr",
               "-deffnm", "rerun", "-ntmpi", "1", "-ntomp", "4",
               "-nb", "cpu", "-pin", "off"], WORK)
    if res.returncode != 0:
        raise SystemExit(res.stderr[-2000:])
    res = run([GMX, "energy", "-f", "rerun.edr", "-o", "energy.xvg"],
              WORK, stdin="Potential\n\n")
    if res.returncode != 0:
        raise SystemExit(res.stderr[-2000:])
    return read_xvg(os.path.join(WORK, "energy.xvg"))


def read_xvg(path):
    cols = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            cols.append([float(p) for p in parts[1:]])
    return np.array(cols)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", default=os.path.join(PROJECT, "phase0", "systems", "npt_3.0.gro"))
    parser.add_argument("--top", default=os.path.join(PROJECT, "phase0", "systems", "topol_3.0.top"))
    parser.add_argument("--radii", default="0.2,0.3,0.4,0.5,0.6,0.7")
    parser.add_argument("--proposals", type=int, default=40)
    parser.add_argument("--sigma-t", type=float, default=0.01)
    parser.add_argument("--sigma-r", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    prepare(args.gro, args.top)
    names, coords, box = read_gro(args.gro)
    waters = water_index(names)
    oxygen = np.array([w[0] for w in waters])
    print(f"atoms {len(names)}  waters {len(waters)}  box {box}")

    rng = np.random.default_rng(args.seed)
    radii = [float(r) for r in args.radii.split(",")]

    frames, labels = [], []
    frames.append(coords.copy())
    labels.append((0.0, 0.0))
    for radius in radii:
        inner = rng.uniform(radius + 0.3, min(box) - radius - 0.3, size=3)
        dist = np.linalg.norm(min_image(coords[oxygen] - inner, box), axis=1)
        patch_oxy = oxygen[dist <= radius]
        for _ in range(args.proposals):
            new, n = make_proposal(coords, waters, patch_oxy, rng,
                                   args.sigma_t, args.sigma_r)
            frames.append(new)
            labels.append((radius, n))

    path = os.path.join(WORK, "frames.gro")
    with open(path, "w") as handle:
        for i, frame in enumerate(frames):
            write_gro_frame(handle, f"frame {i}", names, frame, box)

    pot = energies(path)[:, 0]
    base = pot[0]
    print(f"frames {len(pot)}  base potential {base:.3f} kJ/mol")
    print(f"repeat-frame spread {pot[0] - base:.4f}")

    lines = []
    lines.append(f"system {os.path.basename(args.gro)}  waters {len(waters)}")
    lines.append(f"sigma_t {args.sigma_t} nm  sigma_r {args.sigma_r} rad  kT {KT:.4f} kJ/mol")
    lines.append("")
    lines.append(f"{'radius':>7} {'n_waters':>8} {'dU/kT med':>10} {'A_med':>10} {'A_mean':>10}")
    by_radius = {}
    idx = 1
    for radius in radii:
        subset = labels[idx:idx + args.proposals]
        idx += args.proposals
        du = np.array([pot[idx - args.proposals + j] - base for j in range(args.proposals)])
        n = np.mean([s[1] for s in subset])
        med = float(np.median(du / KT))
        acc = np.mean(np.minimum(1.0, np.exp(-du / KT)))
        by_radius[radius] = (n, med, acc)
        lines.append(f"{radius:7.2f} {n:8.1f} {med:10.3f} {math.exp(-max(0.0, med)):10.3e} {acc:10.3e}")

    ns = np.array([by_radius[r][0] for r in radii])
    dus = np.array([by_radius[r][1] for r in radii])
    slope, intercept = np.polyfit(ns, dus, 1)
    lines.append("")
    lines.append(f"fit dU/kT = {slope:.3f} * n {intercept:+.3f}")

    out = os.path.join(HERE, "patch.txt")
    with open(out, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
