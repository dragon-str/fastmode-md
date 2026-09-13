#!/usr/bin/env python3
"""Live progress view for fastmode audits.

Serves one auto-refreshing page on localhost.  It reads only files that the
audits already write:

- <root>/results.json          completed candidates and the reference
- <root>.log                   the audit's stdout, for the in-flight candidate
- <root>/runs/*/<label>/run.log  step progress of one candidate

Standard library only.  Start it, then open the printed URL.

    python3 dashboard.py --root probe/villin50 --root probe/audit200
"""

import argparse
import glob
import html
import json
import os
import re
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))

VERDICT_RE = re.compile(r"^\s{4}(PASS|fail|rejected)")
LABEL_RE = re.compile(r"^\s{2}(\S+)\s+\.\.\.(?:\s+load\s+([\d.]+))?")
STEP_RE = re.compile(r"^\s*(\d+)\s+[-\d.eE+]+\s*$")


def read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return ""


def waiting():
    return "waiting for load" in read(os.path.join(HERE, "out",
                                                   "run_timed.status"))


def parse_log(path):
    """Return (events, current).  events: list of [label, verdict, detail]."""
    events = []
    current = None
    for line in read(path).splitlines():
        m = LABEL_RE.match(line)
        if m:
            current = [m.group(1), "running", m.group(2) or ""]
            events.append(current)
            continue
        v = VERDICT_RE.match(line)
        if v and current is not None:
            kind = v.group(1)
            if kind == "rejected":
                current[1] = "rejected"
                current[2] = line.split("rejected:", 1)[-1].strip()
            else:
                current[1] = kind
                rest = line.strip()[len(kind):].strip()
                current[2] = rest.strip("() ")
    return events, current


def candidate_dir(root, label):
    hits = glob.glob(os.path.join(root, "runs", "*", label))
    return hits[0] if hits else None


def progress(root, label):
    d = candidate_dir(root, label)
    if not d:
        return None
    mdp = read(os.path.join(d, "run.mdp"))
    m = re.search(r"^\s*nsteps\s*=\s*(\d+)", mdp, re.M)
    nsteps = int(m.group(1)) if m else None
    log = read(os.path.join(d, "run.log"))
    if "Finished mdrun" in log:
        return {"done": True, "step": nsteps, "nsteps": nsteps,
                "nsday": None}
    step = None
    for line in log.splitlines():
        s = STEP_RE.match(line)
        if s:
            step = int(s.group(1))
    nsday = None
    pm = re.search(r"Performance:\s+([\d.]+)", log)
    if pm:
        nsday = float(pm.group(1))
    if step is None:
        return None
    return {"done": False, "step": step, "nsteps": nsteps, "nsday": nsday}


def metrics_html(rec):
    if rec.get("rejected"):
        return f'<span class="badge rej">rejected</span> {html.escape(rec.get("reasons") or rec.get("reason") or "")}'
    if rec.get("drift") is None:
        return '<span class="badge rej">run failed</span>'
    verdict = rec.get("pass")
    badge = ('<span class="badge pass">PASS</span>' if verdict
             else '<span class="badge fail">FAIL</span>')
    ns = rec.get("ns_per_day")
    ns_s = f"{ns:8.1f}" if ns else "  untimed"
    rdf = rec.get("rdf_dev")
    rdf_s = f"{rdf:.4f}" if rdf is not None else " n/a"
    reasons = "".join(f"<div class='reason'>{html.escape(r)}</div>"
                      for r in (rec.get("reasons") or []))
    return ("<tr><td>" + html.escape(rec["label"]) + "</td>"
            f"<td>{rec.get('drift', float('nan')):9.4f}</td>"
            f"<td>{rec.get('temperature', float('nan')):7.2f}</td>"
            f"<td>{rec.get('density', float('nan')):9.4f}</td>"
            f"<td>{rdf_s}</td><td>{ns_s}</td><td>{badge}{reasons}</td></tr>")


def job_html(name, root):
    root = os.path.abspath(root)
    log_path = os.path.abspath(root) + ".log"
    if not os.path.exists(log_path):
        alt = os.path.join(HERE, os.path.basename(root) + ".log")
        log_path = alt if os.path.exists(alt) else log_path
    results_path = os.path.join(root, "results.json")
    parts = [f"<h2>{html.escape(name)} <span class='path'>{html.escape(os.path.relpath(root, HERE))}</span></h2>"]

    data = None
    if os.path.exists(results_path):
        try:
            data = json.load(open(results_path))
        except (OSError, ValueError):
            data = None

    if data:
        parts.append(
            f"<div class='meta'>system {html.escape(str(data.get('system','')))} &nbsp; "
            f"atoms {data.get('atoms')} &nbsp; {len(data.get('runs',[]))} runs</div>")
        best = data.get("best")
        if best:
            parts.append(f"<div class='best'>fastest passing: "
                         f"{html.escape(best['label'])}</div>")
        parts.append("<table><tr><th>setting</th><th>drift</th><th>T</th>"
                     "<th>rho</th><th>rdfdev</th><th>ns/day</th>"
                     "<th>verdict</th></tr>")
        for rec in data.get("runs", []):
            parts.append(metrics_html(rec))
        parts.append("</table>")

    events, current = parse_log(log_path)
    if current and current[1] == "running":
        pr = progress(root, current[0])
        if pr and pr.get("nsteps"):
            pct = 100.0 * (pr["step"] or 0) / pr["nsteps"]
            bar = (f"<div class='bar'><span style='width:{pct:.1f}%'></span></div>"
                   f"<div class='meta'>step {pr['step']}/{pr['nsteps']} "
                   f"({pct:.1f}%)"
                   + (f" &nbsp; {pr['nsday']:.1f} ns/day" if pr.get("nsday") else "")
                   + "</div>")
        elif pr and pr.get("done"):
            bar = "<div class='meta'>candidate finished, validating</div>"
        else:
            bar = ("<div class='bar'><span style='width:2%'></span></div>"
                   "<div class='meta'>starting (grompp or first step)</div>")
        parts.append(f"<div class='running'>running: "
                     f"{html.escape(current[0])}{bar}</div>")

    if not data and not events:
        done = None
        for candidate in ("report.txt", "selfcheck.txt"):
            if os.path.exists(os.path.join(root, candidate)):
                done = candidate
                break
        if done:
            parts.append(f"<div class='meta'>done; full output in "
                         f"{html.escape(done)} (this tool does not write "
                         f"results.json)</div>")
        elif waiting():
            parts.append("<div class='meta queued'>queued; the timed sweep "
                         "starts when load drops below 2.0</div>")
        else:
            parts.append("<div class='meta'>no results yet</div>")
    return "".join(parts)


PAGE = """<!doctype html><html><head><meta charset="utf-8">
<meta http-equiv="refresh" content="5">
<title>fastmode progress</title>
<style>
 body {{ font: 13px/1.5 -apple-system, system-ui, sans-serif; margin: 20px;
        background:#111; color:#ddd; }}
 h1 {{ font-size: 18px; }} h2 {{ font-size: 15px; margin-top: 22px;
      border-bottom:1px solid #333; padding-bottom:4px; }}
 .path {{ color:#777; font-weight:400; font-size:12px; }}
 .meta {{ color:#999; margin:4px 0; }}
 table {{ border-collapse: collapse; margin:8px 0; font-family: ui-monospace,
          monospace; font-size:12px; }}
 th, td {{ padding:3px 10px; text-align:right; border-bottom:1px solid #222; }}
 th:first-child, td:first-child {{ text-align:left; }}
 th {{ color:#aaa; font-weight:600; }}
 .badge {{ padding:1px 6px; border-radius:3px; font-weight:700;
           font-size:11px; }}
 .pass {{ background:#164; color:#8f8; }} .fail {{ background:#611; color:#f99; }}
 .rej {{ background:#333; color:#aaa; }}
 .reason {{ color:#f99; font-size:11px; }}
 .best {{ color:#8f8; margin:4px 0; }}
 .running {{ margin:8px 0; color:#8cf; }}
 .bar {{ height:8px; background:#222; border-radius:4px; margin:4px 0;
         width:420px; }}
 .bar span {{ display:block; height:100%; background:#4af; border-radius:4px; }}
 .head {{ color:#888; }}
 .banner {{ background:#332; color:#fe9; padding:6px 10px; border-radius:4px;
            margin:10px 0; }}
 .queued {{ color:#fe9; }}
</style></head><body>
<h1>fastmode progress</h1>
<div class="head">{head}</div>
{body}
</body></html>"""


class Handler(BaseHTTPRequestHandler):
    roots = []

    def do_GET(self):
        load = os.getloadavg()
        head = (f"{time.strftime('%H:%M:%S')} &nbsp; load "
                f"{load[0]:.2f} {load[1]:.2f} {load[2]:.2f} &nbsp; "
                f"refresh 5s")
        body = "".join(job_html(name, root) for name, root in self.roots)
        if waiting():
            body = ("<div class='banner'>queued: run_timed.sh waits for load "
                    "below 2.0, then the timed sweep starts by itself</div>"
                    + body)
        page = PAGE.format(head=head, body=body)
        data = page.encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, *a):
        pass


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", action="append", default=[],
                    help="audit output dir or its stdout log; repeatable")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()

    roots = []
    for r in args.root:
        r = os.path.abspath(r)
        if r.endswith(".log"):
            roots.append((os.path.basename(r)[:-4], r[:-4]))
        else:
            roots.append((os.path.basename(r), r))
    if not roots:
        for d in sorted(glob.glob(os.path.join(HERE, "out", "*"))):
            if os.path.isdir(d):
                roots.append((os.path.basename(d), d))
        for d in [os.path.join(HERE, "probe", "villin50"),
                  os.path.join(HERE, "probe", "audit200")]:
            if os.path.isdir(d):
                roots.append((os.path.basename(d), d))
    Handler.roots = roots

    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"dashboard on http://127.0.0.1:{args.port}/  jobs: "
          + ", ".join(n for n, _ in roots), flush=True)
    srv.serve_forever()


if __name__ == "__main__":
    main()
