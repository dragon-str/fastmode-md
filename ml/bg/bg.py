"""Boltzmann generator prototype for flexible water.

Representation: min-image displacement of every atom from atom 0 (anchor).
The anchor is fixed at the box centre, so a configuration is 3*(N-1) numbers
on a bounded cube. A RealNVP flow is trained on MD frames. The flow is used as
an independence proposal in Metropolis-Hastings, scored by the exact GROMACS
potential. Equilibrium properties only.
"""
import argparse
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
GMX = os.environ.get("GMX", os.path.join(PROJECT, "build", "bin", "gmx"))
GMXLIB = os.environ.get("GMXLIB", os.path.join(PROJECT, "gromacs", "share", "top"))
KT = 2.4943


def run(cmd, cwd, stdin=None):
    env = dict(os.environ)
    env["GMXLIB"] = GMXLIB
    return subprocess.run(cmd, cwd=cwd, input=stdin, text=True,
                          capture_output=True, env=env)


def parse_gro(path):
    lines = open(path).read().splitlines()
    n = int(lines[1].split()[0])
    resname, name, resid, coords = [], [], [], []
    for ln in lines[2:2 + n]:
        resname.append(ln[5:10].strip())
        name.append(ln[10:15].strip())
        resid.append(int(ln[0:5]))
        coords.append([float(ln[20:28]), float(ln[28:36]), float(ln[36:44])])
    box = [float(v) for v in lines[2 + n].split()[:3]]
    return resname, name, resid, np.array(coords), np.array(box)


def read_frames(path):
    lines = open(path).read().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i + 1].split()[0])
        coords = []
        for ln in lines[i + 2:i + 2 + n]:
            coords.append([float(ln[20:28]), float(ln[28:36]), float(ln[36:44])])
        box = [float(v) for v in lines[i + 2 + n].split()[:3]]
        frames.append((np.array(coords), np.array(box)))
        i += 3 + n
    return frames


def to_repr(coords, box):
    d = coords[1:] - coords[0]
    d -= box * np.round(d / box)
    return d.ravel()


def from_repr(z, box):
    d = z.reshape(-1, 3)
    anchor = box / 2.0
    coords = np.empty((len(d) + 1, 3))
    coords[0] = anchor
    coords[1:] = anchor + d
    return coords


def write_gro_frame(fh, resname, name, resid, coords, box, title="frame"):
    n = len(name)
    fh.write(f"{title}\n{n}\n")
    for i in range(n):
        fh.write(f"{resid[i]:5d}{resname[i]:<5s}{name[i]:>5s}{i + 1:5d}"
                 f"{coords[i,0]:8.3f}{coords[i,1]:8.3f}{coords[i,2]:8.3f}\n")
    fh.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


def energies(workdir, frames, resname, name, resid, box):
    """Exact GROMACS potential for each frame, one mdrun -rerun on a batch."""
    gro = os.path.join(workdir, "batch.gro")
    with open(gro, "w") as fh:
        for c in frames:
            write_gro_frame(fh, resname, name, resid, c, box)
    g = run([GMX, "trjconv", "-f", "batch.gro", "-s", "run.tpr",
             "-o", "batch.trr"], workdir, stdin="0\n")
    if g.returncode != 0:
        raise RuntimeError(g.stderr[-400:])
    m = run([GMX, "mdrun", "-s", "run.tpr", "-rerun", "batch.trr",
             "-deffnm", "rerun", "-ntmpi", "1", "-ntomp", "4", "-nb", "cpu",
             "-pin", "off"], workdir)
    if m.returncode != 0:
        raise RuntimeError(m.stderr[-400:])
    e = run([GMX, "energy", "-f", "rerun.edr", "-o", "batch_ener.xvg"],
            workdir, stdin="Potential\n\n")
    if e.returncode != 0:
        raise RuntimeError(e.stderr[-400:])
    vals = []
    for ln in open(os.path.join(workdir, "batch_ener.xvg")):
        if ln[0] in "@#":
            continue
        parts = ln.split()
        if len(parts) >= 2:
            vals.append(float(parts[1]))
    return np.array(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(HERE, "build"))
    ap.add_argument("--work", default=os.path.join(HERE, "bgwork"))
    ap.add_argument("--frames", default=os.path.join(HERE, "traj.gro"))
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--hidden", type=int, default=256)
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--proposals", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--clip", type=float, default=None)
    ap.add_argument("--load", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    import shutil
    shutil.copy(os.path.join(args.build, "prod.tpr"), os.path.join(args.work, "run.tpr"))

    resname, name, resid, coords, box = parse_gro(os.path.join(args.build, "em.gro"))
    fr = read_frames(args.frames)
    print(f"{len(fr)} frames, {len(name)} atoms, box {box[0]:.3f}")
    X = np.array([to_repr(c, b) for c, b in fr])
    X = X[np.isfinite(X).all(1)]
    mu = X.mean(0)
    sd = X.std(0) + 1e-6
    Z = (X - mu) / sd
    D = Z.shape[1]
    print(f"repr dim {D}")

    import torch
    import torch.nn as nn
    torch.manual_seed(args.seed)

    class Coupling(nn.Module):
        def __init__(self, dim, mask, hidden):
            super().__init__()
            self.register_buffer("mask", mask)
            inn = int(mask.sum())
            out = dim - inn
            self.net = nn.Sequential(
                nn.Linear(inn, hidden), nn.ReLU(),
                nn.Linear(hidden, hidden), nn.ReLU(),
                nn.Linear(hidden, 2 * out))
        def forward(self, x):
            xa = x[:, self.mask]
            st = self.net(xa)
            out = x.shape[1] - xa.shape[1]
            s, t = st[:, :out], st[:, out:]
            s = torch.tanh(s) * 0.7
            xb = x[:, ~self.mask]
            yb = xb * torch.exp(s) + t
            y = torch.empty_like(x)
            y[:, self.mask] = xa
            y[:, ~self.mask] = yb
            return y, s.sum(1)

    class Flow(nn.Module):
        def __init__(self, dim, hidden, layers):
            super().__init__()
            self.coups = nn.ModuleList()
            for i in range(layers):
                mask = torch.zeros(dim, dtype=torch.bool)
                mask[i % 2::2] = True
                self.coups.append(Coupling(dim, mask, hidden))
        def forward(self, x):
            logdet = torch.zeros(x.shape[0])
            for c in self.coups:
                x, ld = c(x)
                logdet = logdet + ld
            return x, logdet
        def log_prob(self, x):
            z, logdet = self.forward(x)
            base = -0.5 * (z ** 2).sum(1) - 0.5 * x.shape[1] * np.log(2 * np.pi)
            return base + logdet
        def sample(self, n):
            z = torch.randn(n, self_dim(self))
            x, _ = self.forward_inv(z)
            return x
        def forward_inv(self, z):
            for c in reversed(self.coups):
                z = invert(c, z)
            return z, None

    def self_dim(m):
        return m.coups[0].mask.shape[0]

    def invert(c, y):
        ya = y[:, c.mask]
        st = c.net(ya)
        out = y.shape[1] - ya.shape[1]
        s, t = st[:, :out], st[:, out:]
        s = torch.tanh(s) * 0.7
        yb = y[:, ~c.mask]
        xb = (yb - t) * torch.exp(-s)
        x = torch.empty_like(y)
        x[:, c.mask] = ya
        x[:, ~c.mask] = xb
        return x

    flow = Flow(D, args.hidden, args.layers)
    ckpt = os.path.join(args.work, "flow.pt")
    if args.load and os.path.exists(ckpt):
        ck = torch.load(ckpt, weights_only=False)
        flow.load_state_dict(ck["state"])
        print("loaded flow.pt")
    opt = torch.optim.Adam(flow.parameters(), lr=1e-3, weight_decay=1e-5)
    data = torch.tensor(Z, dtype=torch.float32)
    n = len(data)
    bs = 128
    n_ep = 0 if (args.load and os.path.exists(ckpt)) else args.epochs
    for ep in range(n_ep):
        perm = torch.randperm(n)
        tot = 0.0
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            loss = -flow.log_prob(data[idx]).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss) * len(idx)
        if ep % 25 == 0 or ep == args.epochs - 1:
            print(f"epoch {ep} nll {tot / n:.3f}")
    with torch.no_grad():
        xs = torch.tensor(Z[:8], dtype=torch.float32)
        zs, _ = flow.forward(xs)
        back, _ = flow.forward_inv(zs)
        print("roundtrip max err", float((back - xs).abs().max()))
        for i, c in enumerate(flow.coups):
            st = c.net(xs[:, c.mask])
            out = st.shape[1] // 2
            print(f"layer {i} s[min,max] ({float(st[:,:out].min()):.2f},{float(st[:,:out].max()):.2f}) t_absmax {float(st[:,out:].abs().max()):.3f}")
        ss = flow.sample(500)
        print(f"sample std mean {float(ss.std(0).mean()):.3f} max {float(ss.std(0).max()):.3f} absmax {float(ss.abs().max()):.1f}")
    torch.save({"state": flow.state_dict(), "mu": mu, "sd": sd, "D": D,
                "hidden": args.hidden, "layers": args.layers},
               os.path.join(args.work, "flow.pt"))

    # Metropolis-Hastings, independence sampler
    with torch.no_grad():
        Zp = flow.sample(args.proposals)
        if args.clip:
            Zp = Zp.clamp(-args.clip, args.clip)
        Zp = Zp.numpy()
    Xp = Zp * sd + mu
    cand = [from_repr(z, box) for z in Xp]
    U = energies(args.work, cand, resname, name, resid, box)
    print(f"evaluated {len(U)} proposals")

    rng = np.random.default_rng(args.seed)
    with torch.no_grad():
        logq_cand = flow.log_prob(torch.tensor(Zp, dtype=torch.float32)).numpy()

    # start from an MD frame
    z0 = Z[0]
    x0 = from_repr(z0 * sd + mu, box)
    u0 = energies(args.work, [x0], resname, name, resid, box)[0]
    with torch.no_grad():
        logq0 = float(flow.log_prob(torch.tensor(z0[None], dtype=torch.float32))[0])
    accepted = 0
    chain = [x0]
    for j in range(len(U)):
        loga = -(U[j] - u0) / KT + (logq0 - logq_cand[j])
        if np.log(rng.random()) < loga:
            u0 = U[j]
            logq0 = logq_cand[j]
            z0 = Zp[j]
            accepted += 1
            chain.append(cand[j])
    print(f"acceptance {accepted}/{len(U)} = {accepted / max(1, len(U)):.3f}")

    # energy comparison
    print(f"proposal U mean {U.mean():.1f} std {U.std():.1f}")
    if len(chain) > 1:
        cu = energies(args.work, chain, resname, name, resid, box)
        print(f"chain    U mean {cu.mean():.1f} std {cu.std():.1f} n {len(chain)}")
    np.save(os.path.join(args.work, "proposal_U.npy"), U)
    np.save(os.path.join(args.work, "chain.npy"), np.array(chain[:200]))


if __name__ == "__main__":
    main()
