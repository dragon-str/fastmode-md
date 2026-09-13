#!/usr/bin/env python3
"""search_dashboard: a live view of a running `search.py` sweep.

Pure standard library.  Reads `search.json` from one or more search roots and
renders an auto-refreshing page with:

  - the best certified speedup so far, against the 2 fs reference;
  - elapsed time, ETA and the two progress bars (screen, certify);
  - the live step progress of the trial that is running now;
  - a map of the search space (timestep against mass factor, one panel for each
    constraint set) where each tested point is coloured by its status;
  - a scatter of ns/day against timestep for every tested point.

Usage:  python3 search_dashboard.py --root <search out dir> [--root ...]
        [--port 8780]
"""

import argparse
import glob
import html
import json
import os
import re
import time
from http.server import BaseHTTPRequestHandler, HTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))

GREEN = "#2e9e4f"
LTGREEN = "#8fd19e"
AMBER = "#e0a400"
RED = "#d1495b"
GRAY = "#9aa0a6"
BLUE = "#3b7dd8"

STEP_RE = re.compile(r"^\s*(\d+)\s+([0-9.]+)\s*$")


def load_roots(paths):
    out = []
    for p in paths:
        p = os.path.abspath(p)
        if p.endswith(".json") and os.path.exists(p):
            p = os.path.dirname(p)
        if os.path.isdir(p):
            out.append(p)
    return out


def read_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def live_step(workdir, stage, space):
    """Last step number and target from the running run.log."""
    if not workdir:
        return None
    log = os.path.join(workdir, "run.log")
    if not os.path.exists(log):
        return None
    last = None
    try:
        with open(log, errors="ignore") as fh:
            for line in fh:
                m = STEP_RE.match(line)
                if m:
                    last = int(m.group(1))
    except OSError:
        return None
    if last is None:
        return None
    return last


def target_steps(label, stage, space):
    m = re.search(r"dt(\d+)fs", label or "")
    if not m:
        return None
    dt = int(m.group(1)) / 1000.0
    ps = space.get("screen_ps", 20.0) if stage == "screen" \
        else space.get("ref_ps", 200.0)
    return 2000 + int(round(ps / dt))


def fmt_duration(seconds):
    if seconds is None or seconds < 0:
        return "--"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h {m:02d}m"
    if m:
        return f"{m}m {s:02d}s"
    return f"{s}s"


def status_of(rec):
    if rec.get("stage") == "certify":
        if rec.get("certified"):
            return "certified"
        return "cert-fail"
    if rec.get("rejected"):
        return "rejected"
    if rec.get("screen_pass"):
        return "screen-pass"
    return "screen-fail"


COLORS = {"certified": GREEN, "cert-fail": RED, "screen-pass": LTGREEN,
          "screen-fail": "#e8a0a8", "rejected": GRAY}


def eta(state):
    times = state.get("screen_times") or []
    avg = sum(times) / len(times) if times else 30.0
    space = state.get("space", {})
    ratio = space.get("ref_ps", 200.0) / space.get("screen_ps", 20.0)
    rem_screen = max(0, state.get("total_screen", 0) - state.get("done_screen", 0))
    rem_cert = max(0, state.get("total_certify", 0) - state.get("done_certify", 0))
    return rem_screen * avg + rem_cert * avg * ratio


def best_speedup(state):
    best = state.get("best")
    ref = state.get("cert_ref") or state.get("screen_ref")
    if not best or not ref or not best.get("ns_per_day"):
        return None
    return best["ns_per_day"] / ref["ns_per_day"]


def svg_space_map(trials, space):
    """Search space: dt on x, mass factor on y, one panel per constraint."""
    dts = sorted({round(t["spec"]["dt"] * 1000) for t in trials}) or [4, 5]
    factors = sorted({int(t["spec"].get("hmr_factor", 0))
                      for t in trials}) or [0]
    cons = sorted({t["spec"].get("constraints", "h-bonds") for t in trials}) \
        or ["h-bonds"]
    pw, ph = 260, 190
    pad = 34
    W = pad + len(cons) * (pw + 12)
    H = ph + 56
    x0 = lambda r: pad + r
    fy = lambda f: H - 40 - (factors.index(f) + 0.5) * ((ph - 20) / max(1, len(factors)))
    cx = lambda dt: x0(0) + 14 + (dts.index(dt) + 0.5) * ((pw - 28) / max(1, len(dts)))
    parts = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    parts.append(f'<text x="{pad}" y="16" class="svgt">search space: timestep x mass factor</text>')
    for ci, c in enumerate(cons):
        px = pad + ci * (pw + 12)
        parts.append(f'<rect x="{px}" y="24" width="{pw}" height="{ph}" rx="6" '
                     f'fill="#fbfbfd" stroke="#e3e5ea"/>')
        parts.append(f'<text x="{px + pw/2}" y="{ph + 14}" text-anchor="middle" '
                     f'class="svgl">{html.escape(c)}</text>')
        parts.append(f'<text x="{px}" y="{ph + 30}" text-anchor="start" '
                     f'class="svgax">mass factor vs dt (fs)</text>')
        for dt in dts:
            gx = px + 14 + (dts.index(dt) + 0.5) * ((pw - 28) / max(1, len(dts)))
            parts.append(f'<text x="{gx}" y="{24 + ph - 4}" text-anchor="middle" '
                         f'class="svgax">{dt}</text>')
        for f in factors:
            gy = 24 + (ph - 20) - (factors.index(f) + 0.5) * ((ph - 20) / max(1, len(factors)))
            parts.append(f'<text x="{px - 4}" y="{gy + 3}" text-anchor="end" '
                         f'class="svgax">{f or "off"}</text>')
        for t in trials:
            if t["spec"].get("constraints", "h-bonds") != c:
                continue
            dt = round(t["spec"]["dt"] * 1000)
            f = int(t["spec"].get("hmr_factor", 0))
            if dt not in dts or f not in factors:
                continue
            gx = px + 14 + (dts.index(dt) + 0.5) * ((pw - 28) / max(1, len(dts)))
            gy = 24 + (ph - 20) - (factors.index(f) + 0.5) * ((ph - 20) / max(1, len(factors)))
            st = status_of(t)
            r = 5
            if t.get("stage") == "certify" and t.get("certified"):
                r = 7
            parts.append(f'<circle cx="{gx}" cy="{gy}" r="{r}" '
                         f'fill="{COLORS[st]}" stroke="#333" stroke-width="0.5">'
                         f'<title>{html.escape(t["label"])} {st}</title></circle>')
    parts.append("</svg>")
    return "".join(parts)


def svg_perf(trials, ref):
    """ns/day against dt, all tested points, with the reference line."""
    pts = [t for t in trials if t.get("ns_per_day")]
    W, H, pad = 640, 260, 46
    if not pts:
        return '<p class="muted">no completed runs yet</p>'
    xmax = max(t["spec"]["dt"] * 1000 for t in pts)
    ymax = max([t["ns_per_day"] for t in pts] + ([ref["ns_per_day"]] if ref else [0])) * 1.1
    xmin = min(t["spec"]["dt"] * 1000 for t in pts)
    cx = lambda dt: pad + (dt - xmin) / max(1e-9, xmax - xmin) * (W - pad - 20)
    cy = lambda v: H - pad - v / max(1e-9, ymax) * (H - pad - 20)
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    p.append(f'<text x="{pad}" y="16" class="svgt">performance: ns/day against timestep</text>')
    p.append(f'<line x1="{pad}" y1="{H-pad}" x2="{W-20}" y2="{H-pad}" stroke="#ccc"/>')
    p.append(f'<line x1="{pad}" y1="24" x2="{pad}" y2="{H-pad}" stroke="#ccc"/>')
    for gx in range(int(xmin), int(xmax) + 1):
        p.append(f'<text x="{cx(gx)}" y="{H-pad+16}" text-anchor="middle" '
                 f'class="svgax">{gx}</text>')
    if ref and ref.get("ns_per_day"):
        y = cy(ref["ns_per_day"])
        p.append(f'<line x1="{pad}" y1="{y}" x2="{W-20}" y2="{y}" '
                 f'stroke="{BLUE}" stroke-dasharray="5 4"/>')
        p.append(f'<text x="{W-22}" y="{y-4}" text-anchor="end" class="svgax" '
                 f'fill="{BLUE}">2 fs reference {ref["ns_per_day"]:.0f}</text>')
    for t in pts:
        st = status_of(t)
        gx = cx(t["spec"]["dt"] * 1000)
        gy = cy(t["ns_per_day"])
        r = 6 if st == "certified" else 4
        p.append(f'<circle cx="{gx:.1f}" cy="{gy:.1f}" r="{r}" fill="{COLORS[st]}" '
                 f'fill-opacity="0.85" stroke="#333" stroke-width="0.4">'
                 f'<title>{html.escape(t["label"])}: {t["ns_per_day"]:.0f} ns/day</title>'
                 f'</circle>')
    p.append("</svg>")
    return "".join(p)


def render(root):
    state = read_json(os.path.join(root, "search.json"))
    name = html.escape(os.path.basename(root))
    if state is None:
        return (f'<div class="card"><h2>{name}</h2>'
                f'<p class="muted">waiting for search.json ...</p></div>')
    trials = state.get("trials", [])
    ref = state.get("cert_ref") or state.get("screen_ref") or {}
    sp = best_speedup(state)
    eta_s = eta(state)
    started = min([t.get("wall", 0) for t in trials] or [0])
    elapsed = None
    if state.get("current") and state["current"].get("started"):
        first = state["current"]["started"]
        elapsed = time.time() - first
    done = state.get("done_screen", 0) + state.get("done_certify", 0)
    total = state.get("total_screen", 0) + state.get("total_certify", 0)
    pct = 100.0 * state.get("done_screen", 0) / max(1, state.get("total_screen", 1))

    cur = state.get("current")
    step_html = '<p class="muted">idle</p>'
    if cur:
        lab = cur.get("label", "")
        stage = cur.get("stage", "")
        last = live_step(cur.get("workdir"), stage, state.get("space", {}))
        tgt = target_steps(lab, stage, state.get("space", {}))
        if last and tgt:
            fpct = min(100.0, 100.0 * last / tgt)
            step_html = (f'<div class="lbl">{html.escape(lab)} '
                         f'<span class="muted">({stage})</span></div>'
                         f'<div class="bar"><i style="width:{fpct:.1f}%"></i></div>'
                         f'<div class="muted">{last} / {tgt} steps</div>')
        else:
            step_html = (f'<div class="lbl">{html.escape(lab)}</div>'
                         f'<div class="muted">starting (grompp or first step) ...</div>')

    card = [f'<div class="card"><h2>{name}</h2>']
    if sp is not None:
        b = state["best"]
        card.append(f'<div class="hero"><span class="big">{sp:.2f}x</span>'
                    f'<span class="sub">best certified so far</span></div>')
        card.append(f'<p class="mono">{html.escape(b["label"])}</p>')
        card.append(f'<p class="muted">{b.get("ns_per_day", 0):.1f} ns/day vs '
                    f'reference {ref.get("ns_per_day", 0):.1f}</p>')
    else:
        card.append(f'<div class="hero"><span class="big">--</span>'
                    f'<span class="sub">no certified result yet</span></div>')
    card.append(f'<div class="stats">'
                f'<div><span class="k">elapsed</span><span class="v">'
                f'{fmt_duration(elapsed)}</span></div>'
                f'<div><span class="k">ETA (approx)</span><span class="v">'
                f'{fmt_duration(eta_s)}</span></div>'
                f'<div><span class="k">trials</span><span class="v">'
                f'{done} / {total}</span></div>'
                f'<div><span class="k">screen</span><span class="v">'
                f'{state.get("done_screen",0)}/{state.get("total_screen",0)}</span></div>'
                f'<div><span class="k">certify</span><span class="v">'
                f'{state.get("done_certify",0)}/{state.get("total_certify",0)}</span></div>'
                f'</div>')
    card.append(f'<div class="bar big"><i style="width:{pct:.1f}%"></i></div>')
    card.append(f'<h3>now running</h3>{step_html}')
    card.append(f'<h3>search space</h3>{svg_space_map(trials, state.get("space", {}))}')
    card.append(f'<div class="legend">'
                f'<span><i style="background:{GREEN}"></i>certified pass</span>'
                f'<span><i style="background:{LTGREEN}"></i>screen pass</span>'
                f'<span><i style="background:{RED}"></i>certify fail</span>'
                f'<span><i style="background:#e8a0a8"></i>screen fail</span>'
                f'<span><i style="background:{GRAY}"></i>rejected</span>'
                f'</div>')
    card.append(f'<h3>performance</h3>{svg_perf(trials, ref)}')
    tbl = ['<table><tr><th>setting</th><th>stage</th><th>ns/day</th><th>speedup</th>'
           '<th>drift</th><th>T</th><th>rdf</th><th>status</th></tr>']
    for t in sorted(trials, key=lambda r: (r.get("ns_per_day") or 0), reverse=True)[:40]:
        st = status_of(t)
        nd = f'{t["ns_per_day"]:.1f}' if t.get("ns_per_day") else "--"
        spd = ""
        if t.get("ns_per_day") and ref.get("ns_per_day"):
            spd = f'{t["ns_per_day"]/ref["ns_per_day"]:.2f}x'
        drift = f'{t["drift"]:.4f}' if t.get("drift") is not None else "--"
        temp = f'{t["temperature"]:.1f}' if t.get("temperature") else "--"
        rdf = f'{t["rdf_dev"]:.4f}' if t.get("rdf_dev") is not None else "--"
        tbl.append(f'<tr><td class="mono">{html.escape(t["label"])}</td>'
                   f'<td>{t.get("stage","")}</td><td>{nd}</td><td>{spd}</td>'
                   f'<td>{drift}</td><td>{temp}</td><td>{rdf}</td>'
                   f'<td><span class="dot" style="background:{COLORS[st]}"></span>'
                   f'{st}</td></tr>')
    tbl.append("</table>")
    card.append('<details open><summary><h3 style="display:inline">all trials</h3>'
                '</summary>' + "".join(tbl) + '</details>')
    card.append("</div>")
    return "".join(card)


STYLE = """
body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;background:#f2f3f7;
color:#1b1d22;margin:0;padding:24px}
.card{background:#fff;border:1px solid #e3e5ea;border-radius:12px;padding:18px 20px;
margin:0 auto 20px;max-width:980px;box-shadow:0 1px 3px rgba(0,0,0,.04)}
h2{margin:0 0 10px;font-size:16px}h3{margin:18px 0 8px;font-size:13px;
text-transform:uppercase;letter-spacing:.05em;color:#5b6070}
.hero{display:flex;align-items:baseline;gap:12px}
.big{font-size:44px;font-weight:700;color:#1b7f3b}.sub{color:#6b7180}
.stats{display:flex;flex-wrap:wrap;gap:18px;margin:12px 0}
.stats .k{display:block;font-size:11px;color:#8a90a0;text-transform:uppercase}
.stats .v{font-size:18px;font-weight:600}
.bar{background:#eceef3;border-radius:6px;height:10px;overflow:hidden;margin:6px 0}
.bar.big{height:14px}.bar i{display:block;height:100%;background:#3b7dd8}
.mono{font-family:ui-monospace,Menlo,monospace;font-size:12px;word-break:break-all}
.muted{color:#8a90a0;font-size:12px}
.lbl{font-size:13px;font-weight:600;word-break:break-all}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;margin:10px 0}
.legend i,.dot{display:inline-block;width:10px;height:10px;border-radius:50%;
margin-right:5px;vertical-align:middle}
table{border-collapse:collapse;width:100%;font-size:12px}
th,td{text-align:left;padding:5px 8px;border-bottom:1px solid #eef0f4}
th{color:#8a90a0;font-weight:600}
.svgt{font-size:12px;font-weight:700;fill:#1b1d22}
.svgl{font-size:12px;font-weight:600;fill:#5b6070}
.svgax{font-size:10px;fill:#8a90a0}
"""


def _seed_progress(s):
    nsteps = None
    mdp = os.path.join(s, "run.mdp")
    if os.path.exists(mdp):
        for line in open(mdp, errors="ignore"):
            m = re.match(r"\s*nsteps\s*=\s*(\d+)", line)
            if m:
                nsteps = int(m.group(1))
    step = None
    log = os.path.join(s, "run.log")
    if os.path.exists(log):
        for line in open(log, errors="ignore"):
            m = re.match(r"\s*(\d+)\s+\d+\.\d+\s*$", line)
            if m:
                step = int(m.group(1))
    return step, nsteps


def render_validation(root):
    name = os.path.basename(root.rstrip("/")) or root
    cfg_dirs = sorted(d for d in glob.glob(os.path.join(root, "*"))
                      if os.path.isdir(d) and glob.glob(os.path.join(d, "seed_*")))
    if not cfg_dirs:
        return (f'<div class="card"><h2>{html.escape(name)}</h2>'
                '<p class="muted">waiting for validation runs ...</p></div>')
    configs = {}
    oldest = newest = None
    for d in cfg_dirs:
        seeds = sorted(glob.glob(os.path.join(d, "seed_*")))
        rows = []
        for s in seeds:
            log = os.path.join(s, "run.log")
            fin = os.path.exists(log) and "Finished mdrun" in open(log, errors="ignore").read()
            mt = os.path.getmtime(log) if os.path.exists(log) else os.path.getmtime(s)
            rows.append((os.path.basename(s), fin, s, mt))
            oldest = mt if oldest is None or mt < oldest else oldest
            newest = mt if newest is None or mt > newest else newest
        configs[os.path.basename(d)] = rows
    total = len(configs) * max(len(v) for v in configs.values())
    done = sum(1 for v in configs.values() for r in v if r[1])
    running = None
    for cn, rows in configs.items():
        for rs, fin, s, mt in rows:
            if not fin and mt == newest:
                running = (cn, rs, s)
    elapsed = time.time() - oldest if oldest else 0
    eta = elapsed / done * (total - done) if done else 0
    card = [f'<div class="card"><h2>{html.escape(name)} - validation</h2>']
    card.append(f'<div class="stats">'
                f'<div><span class="k">done</span><span class="v">{done}/{total}</span></div>'
                f'<div><span class="k">elapsed</span><span class="v">{fmt_duration(elapsed)}</span></div>'
                f'<div><span class="k">eta</span><span class="v">{fmt_duration(eta)}</span></div>'
                f'</div>')
    pct = 100.0 * done / total if total else 0
    card.append(f'<div class="bar big"><i style="width:{pct:.1f}%"></i></div>')
    if running:
        cn, rs, s = running
        step, nsteps = _seed_progress(s)
        card.append(f'<h3>now running</h3><p class="lbl">{html.escape(cn)} / {html.escape(rs)}</p>')
        if step is not None and nsteps:
            sp = 100.0 * step / nsteps
            card.append(f'<div class="bar"><i style="width:{sp:.1f}%"></i></div>'
                        f'<p class="muted">{step}/{nsteps} steps ({sp:.0f}%)</p>')
        else:
            card.append('<p class="muted">starting (grompp or first step)</p>')
    card.append('<h3>configs</h3><table><tr><th>config</th><th>seeds done</th>'
                '<th>last seed</th><th>status</th></tr>')
    for cn, rows in configs.items():
        nd = sum(1 for r in rows if r[1])
        last = rows[-1][0] if rows else "--"
        st = "running" if running and running[0] == cn else ("done" if nd == len(rows) and rows else "--")
        card.append(f'<tr><td class="mono">{html.escape(cn)}</td><td>{nd}/{len(rows)}</td>'
                    f'<td>{html.escape(last)}</td><td>{st}</td></tr>')
    card.append('</table>')
    card.append('<p class="muted">each seed: 0.3 ns; observables compared against the '
                'periodic 2 fs reference.</p></div>')
    return "".join(card)


class Handler(BaseHTTPRequestHandler):
    roots = []

    def do_GET(self):
        if self.path.startswith("/search.json"):
            root = self.roots[0] if self.roots else None
            data = read_json(os.path.join(root, "search.json")) if root else None
            body = json.dumps(data).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(body)
            return
        body = ("<!doctype html><html><head><meta charset='utf-8'>"
                "<meta http-equiv='refresh' content='3'>"
                "<title>fastmode search</title><style>" + STYLE + "</style></head>"
                "<body><h1 style='text-align:center;font-size:20px'>"
                "fastmode search</h1>" +
                "".join(render(r) if read_json(os.path.join(r, "search.json"))
                        else render_validation(r) for r in self.roots) +
                "</body></html>").encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *a):
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", required=True)
    ap.add_argument("--port", type=int, default=8780)
    args = ap.parse_args()
    Handler.roots = load_roots(args.root)
    print(f"search dashboard: http://127.0.0.1:{args.port}/  roots={Handler.roots}")
    HTTPServer(("127.0.0.1", args.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
