"""Aggregate patchfrontier.txt into median frontiers and plot (phase0/ml)."""

import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    rows = []
    radius = None
    for line in open(path):
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "radius":
            radius = float(parts[1])
            continue
        if parts[0] in ("move", "score", "train", "proposals") or len(parts) < 5:
            continue
        move = parts[0]
        step = float(parts[1])
        rows.append((radius, move, step, float(parts[2]), float(parts[3])))
    return rows


def frontiers(rows):
    groups = defaultdict(list)
    for radius, move, step, rms, du in rows:
        groups[(radius, move, step)].append((rms, du))
    out = defaultdict(list)
    for (radius, move, step), vals in groups.items():
        rms = np.mean([v[0] for v in vals])
        med = float(np.median([v[1] for v in vals]))
        out[(radius, move)].append((rms, med, math.exp(-max(0.0, med))))
    for key in out:
        out[key].sort()
    return out


def main():
    path = os.path.join(HERE, "patchfrontier.txt")
    rows = load(path)
    front = frontiers(rows)
    radii = sorted({r for r, _ in front})
    waters = {}
    for line in open(path):
        if line.startswith("radius") and "waters" in line:
            parts = line.split()
            waters[float(parts[1])] = int(parts[3])

    print(f"{'waters':>7} {'move':>8} {'RMS nm':>9} {'dU/kT':>9} {'accept':>10}")
    for radius in radii:
        for move in ("jiggle", "rotate", "learned"):
            pts = front.get((radius, move), [])
            print(f"--- {waters.get(radius, radius):>4} waters  r={radius:.2f}  {move}")
            for rms, med, acc in pts:
                print(f"{'':7} {move:>8} {rms:9.4f} {med:9.3f} {acc:10.3e}")
    print()
    print("Frontier check: does learned sit above rotate at the same RMS?")
    for radius in radii:
        rot = front.get((radius, "rotate"), [])
        lea = front.get((radius, "learned"), [])
        if not rot or not lea:
            continue
        rms_lea = np.array([p[0] for p in lea])
        rms_rot = np.array([p[0] for p in rot])
        acc_rot = np.array([p[2] for p in rot])
        interp = np.interp(rms_lea, rms_rot, np.log(np.maximum(acc_rot, 1e-300)),
                           left=np.nan, right=np.nan)
        ok = ~np.isnan(interp)
        if ok.any():
            ratio = np.exp(np.array([p[2] for p in lea])[ok] - interp[ok])
            print(f"  {waters.get(radius, radius):>5} waters: "
                  f"learned/rotate acceptance ratio median {np.median(ratio):.3g} "
                  f"(range {ratio.min():.3g} to {ratio.max():.3g})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    moves = [("jiggle", "o", "tab:gray"), ("rotate", "s", "tab:blue"),
             ("learned", "^", "tab:red")]
    fig, axes = plt.subplots(1, len(radii), figsize=(3.6 * len(radii), 4.0), sharey=True)
    if len(radii) == 1:
        axes = [axes]
    for ax, radius in zip(axes, radii):
        for move, marker, colour in moves:
            pts = front.get((radius, move), [])
            if not pts:
                continue
            ax.plot([p[0] for p in pts], [max(p[2], 1e-12) for p in pts],
                    marker=marker, color=colour, label=move, lw=1.5)
        ax.set_yscale("log")
        ax.set_ylim(1e-6, 1.5)
        ax.set_xlabel("RMS displacement (nm)")
        ax.set_title(f"{waters.get(radius, radius)} waters")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("acceptance  exp(-median dU/kT)")
    axes[0].legend()
    fig.suptitle("Local water move efficiency: exact GROMACS scorer, matched patches")
    fig.tight_layout()
    out = os.path.join(HERE, "patchfrontier.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
