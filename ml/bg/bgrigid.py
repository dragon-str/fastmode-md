"""Rigid-body Boltzmann generator prototype.

Each water is an oxygen position (3) plus an orientation rotation vector (3);
hydrogens are rebuilt from the fixed TIP3P geometry. This enforces bond
geometry by construction. The orientation change of variables contributes
J(omega) = 2(1-cos|omega|)/|omega|^2 to the Metropolis ratio.
"""
import argparse
import os
import shutil

import numpy as np
import bg

R_OH = 0.09572
ANG = np.deg2rad(104.52)
_C, _S = np.cos(ANG / 2), np.sin(ANG / 2)
U_REF = np.array([[_C, _C, 0.0],
                  [_S, -_S, 0.0],
                  [0.0, 0.0, -np.sin(ANG)]])


def log_map(R):
    tr = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    th = np.arccos(tr)
    if th < 1e-8:
        return np.zeros(3)
    vee = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    if th < np.pi - 1e-6:
        return th / (2.0 * np.sin(th)) * vee
    A = (R + np.eye(3)) / 2.0
    axis = np.sqrt(np.clip(np.diag(A), 0.0, None))
    k = int(np.argmax(axis))
    v = np.zeros(3)
    v[k] = axis[k]
    for j in range(3):
        if j != k and axis[j] > 1e-8:
            v[j] = A[k, j] / v[k]
    v = v / (np.linalg.norm(v) + 1e-12)
    return np.pi * v


def exp_map(w):
    th = np.linalg.norm(w)
    if th < 1e-8:
        return np.eye(3)
    k = w / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def extract(coords):
    x = coords.reshape(-1, 3, 3)
    O = x[:, 0].copy()
    u = x[:, 1] - O
    v = x[:, 2] - O
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    w = np.cross(u, v)
    w /= np.linalg.norm(w, axis=1, keepdims=True)
    Ulab = np.stack([u, v, w], axis=-1)
    R = Ulab @ U_REF.T
    om = np.array([log_map(r) for r in R])
    return O, om


def rebuild(reprv, box, nw):
    d = reprv[:3 * (nw - 1)].reshape(nw - 1, 3)
    om = reprv[3 * (nw - 1):].reshape(nw, 3)
    O = np.empty((nw, 3))
    O[0] = box / 2.0
    O[1:] = O[0] + d
    R = np.array([exp_map(w) for w in om])
    Ulab = R @ U_REF
    H1 = O + R_OH * Ulab[:, :, 0]
    H2 = O + R_OH * Ulab[:, :, 1]
    return np.stack([O, H1, H2], axis=1).reshape(-1, 3)


def logJ(om):
    th = np.linalg.norm(om, axis=1)
    small = th < 1e-6
    J = np.where(small, 1.0, 2.0 * (1 - np.cos(np.where(small, 1, th))) / np.where(small, 1, th) ** 2)
    return np.log(J).sum()


def to_repr(O, om, box):
    d = O[1:] - O[0]
    d -= box * np.round(d / box)
    return np.concatenate([d.ravel(), om.ravel()])


def make_flow(D, hidden, layers, seed):
    import torch
    import torch.nn as nn
    torch.manual_seed(seed)

    class Coupling(nn.Module):
        def __init__(self, dim, mask, hidden):
            super().__init__()
            self.register_buffer("mask", mask)
            inn = int(mask.sum())
            outd = dim - inn
            self.net = nn.Sequential(
                nn.Linear(inn, hidden), nn.ReLU(),
                nn.Linear(hidden, hidden), nn.ReLU(),
                nn.Linear(hidden, 2 * outd))
        def forward(self, x):
            xa = x[:, self.mask]
            st = self.net(xa)
            outd = x.shape[1] - xa.shape[1]
            s, t = st[:, :outd], st[:, outd:]
            s = torch.tanh(s) * 0.7
            y = torch.empty_like(x)
            y[:, self.mask] = xa
            y[:, ~self.mask] = x[:, ~self.mask] * torch.exp(s) + t
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
            ld = torch.zeros(x.shape[0])
            for c in self.coups:
                x, s = c(x)
                ld = ld + s
            return x, ld
        def log_prob(self, x):
            z, ld = self.forward(x)
            base = -0.5 * (z ** 2).sum(1) - 0.5 * x.shape[1] * np.log(2 * np.pi)
            return base + ld
        def sample(self, n):
            z = torch.randn(n, self.coups[0].mask.shape[0])
            x = z
            for c in reversed(self.coups):
                xa = x[:, c.mask]
                st = c.net(xa)
                outd = x.shape[1] - xa.shape[1]
                s, t = st[:, :outd], st[:, outd:]
                s = torch.tanh(s) * 0.7
                y = torch.empty_like(x)
                y[:, c.mask] = xa
                y[:, ~c.mask] = (x[:, ~c.mask] - t) * torch.exp(-s)
                x = y
            return x

    return Flow(D, hidden, layers)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(bg.HERE, "rigid"))
    ap.add_argument("--work", default=os.path.join(bg.HERE, "rigidwork"))
    ap.add_argument("--frames", default=os.path.join(bg.HERE, "rigid", "traj.gro"))
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--hidden", type=int, default=256)
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--proposals", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--clip", type=float, default=4.0)
    ap.add_argument("--load", default=None)
    args = ap.parse_args()
    os.makedirs(args.work, exist_ok=True)
    shutil.copy(os.path.join(args.build, "rigid.tpr"), os.path.join(args.work, "run.tpr"))

    resname, name, resid, coords, box = bg.parse_gro(os.path.join(args.build, "em.gro"))
    nw = len(name) // 3
    fr = bg.read_frames(args.frames)
    print(f"{len(fr)} frames, {len(name)} atoms, {nw} waters, box {box[0]:.3f}")
    X = []
    for c, b in fr:
        O, om = extract(c)
        X.append(to_repr(O, om, b))
    X = np.array(X)
    mu = X.mean(0)
    sd = X.std(0) + 1e-6
    Z = ((X - mu) / sd).astype(np.float32)
    D = Z.shape[1]
    print(f"repr dim {D}")

    import torch
    flow = make_flow(D, args.hidden, args.layers, args.seed)
    if args.load:
        ck = torch.load(args.load, weights_only=False)
        flow.load_state_dict(ck["state"])
    else:
        opt = torch.optim.Adam(flow.parameters(), lr=1e-3, weight_decay=1e-5)
        data = torch.tensor(Z)
        n = len(data)
        bs = 128
        for ep in range(args.epochs):
            perm = torch.randperm(n)
            tot = 0.0
            for i in range(0, n, bs):
                idx = perm[i:i + bs]
                loss = -flow.log_prob(data[idx]).mean()
                opt.zero_grad()
                loss.backward()
                opt.step()
                tot += float(loss) * len(idx)
            if ep % 50 == 0 or ep == args.epochs - 1:
                print(f"epoch {ep} nll {tot / n:.3f}")
        torch.save({"state": flow.state_dict(), "mu": mu, "sd": sd, "D": D,
                    "hidden": args.hidden, "layers": args.layers},
                   os.path.join(args.work, "flow.pt"))

    with torch.no_grad():
        Zp = flow.sample(args.proposals)
        Zp = torch.clamp(Zp, -args.clip, args.clip).numpy()
    Xp = Zp * sd + mu
    cand = [rebuild(x, box, nw) for x in Xp]
    U = bg.energies(args.work, cand, resname, name, resid, box)
    with torch.no_grad():
        logq = flow.log_prob(torch.tensor(Zp)).numpy()
    logj = np.array([logJ(x[3 * (nw - 1):].reshape(nw, 3)) for x in Xp])

    rng = np.random.default_rng(args.seed)
    z0 = ((X[0] - mu) / sd).astype(np.float32)
    x0 = rebuild(z0 * sd + mu, box, nw)
    u0 = bg.energies(args.work, [x0], resname, name, resid, box)[0]
    with torch.no_grad():
        logq0 = float(flow.log_prob(torch.tensor(z0[None]))[0])
    logj0 = logJ(z0[3 * (nw - 1):].reshape(nw, 3))
    acc = 0
    chain = [x0]
    for j in range(len(U)):
        loga = -(U[j] - u0) / bg.KT + (logq0 - logq[j]) + (logj[j] - logj0)
        if np.log(rng.random()) < loga:
            u0 = U[j]
            logq0 = logq[j]
            logj0 = logj[j]
            acc += 1
            chain.append(cand[j])
    print(f"acceptance {acc}/{len(U)} = {acc / max(1, len(U)):.4f}")
    print(f"proposal U mean {U.mean():.1f} std {U.std():.1f}")
    np.save(os.path.join(args.work, "proposal_U.npy"), U)
    np.save(os.path.join(args.work, "chain.npy"), np.array(chain))


if __name__ == "__main__":
    main()
