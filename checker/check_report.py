#!/usr/bin/env python3
"""Mechanical gate for scOpAlign executor reports.

Usage:
    python checker/check_report.py tasks/T01.md reports/T01/report.json [--prev reports/T00/report.json]

Exit 0 = PASS, 1 = FAIL. Prints every failure; never stops at the first.

Checks (all mechanical, no judgment):
  1. report.json validates against templates/report.schema.json
  2. task_id in report matches the task file
  3. gated keys == exactly the names in the task's pre-registration table
  4. every deliverable listed in the task exists in reports/<id>/ and its
     sha256 matches the report
  5. report.json, requirements.lock exist
  6. lockfile sha256 matches report; if --prev given, flag any diff in
     pinned package versions
  7. known_failures == KNOWN_FAILURES in coord/PROJECT.md (both directions)
  8. every gated name has a regeneration command
  9. observations contain no interpretive language (soft check, warns)
"""
import argparse, ast, hashlib, json, re, sys
from pathlib import Path

try:
    import jsonschema
except ImportError:
    sys.exit("pip install jsonschema")

ROOT = Path(__file__).resolve().parent.parent
INTERPRETIVE = re.compile(
    r"\b(suggest|implies|indicates|therefore|because|likely|probably|means that|"
    r"consistent with|supports|confirms|we think|I think|hypothes)\w*", re.I)

def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_task(task_md: Path):
    text = task_md.read_text()
    tid = re.match(r"#\s*(T\d{2}[a-z]?|K\d+[a-z]?)\b", text).group(1)
    # gated table rows: | `name` | ... |
    gated = re.findall(r"^\|\s*`([A-Za-z0-9_]+)`\s*\|", text, re.M)
    # deliverables: lines starting "- `" under ## Deliverables
    m = re.search(r"## Deliverables(.*?)(?:\n## |\Z)", text, re.S)
    deliverables = re.findall(r"-\s*`([^`]+)`", m.group(1)) if m else []
    deliverables = [x for x in deliverables if x != "report.json"]  # cannot checksum itself
    return tid, set(gated), deliverables

def parse_known_failures(project_md: Path):
    # Parse the fenced Python list literally: test IDs contain brackets.
    text = project_md.read_text()
    m = re.search(r"^KNOWN_FAILURES\s*=\s*(\[.*?^\])", text, re.S | re.M)
    if not m:
        m = re.search(r"KNOWN_FAILURES\s*=\s*(\[\])", text)
    if not m:
        raise SystemExit("KNOWN_FAILURES list not found in PROJECT.md")
    return set(ast.literal_eval(m.group(1)))

def lock_versions(p: Path):
    out = {}
    for line in p.read_text().splitlines():
        if "==" in line:
            k, v = line.split("==", 1)
            out[k.strip().lower()] = v.strip()
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("task_md")
    ap.add_argument("report_json")
    ap.add_argument("--prev", help="previous accepted report.json for lockfile diff")
    a = ap.parse_args()

    fails, warns = [], []
    task_md, report_json = Path(a.task_md), Path(a.report_json)
    report_dir = report_json.parent

    # 1 schema
    schema = json.loads((ROOT / "templates" / "report.schema.json").read_text())
    report = json.loads(report_json.read_text())
    try:
        jsonschema.validate(report, schema)
    except jsonschema.ValidationError as e:
        fails.append(f"schema: {e.message}")

    # 2 task id
    tid, task_gated, task_deliv = parse_task(task_md)
    if report.get("task_id") != tid:
        fails.append(f"task_id {report.get('task_id')} != {tid}")

    # 3 gated keys (a pre-registered early stop reports only the numbers reached)
    stopped = report.get("status") == "stopped"
    rep_gated = set(report.get("gated", {}))
    if stopped:
        if not rep_gated or rep_gated - task_gated:
            fails.append(f"stopped report: gated must be a non-empty subset of the task's; extra={rep_gated-task_gated}")
    elif rep_gated != task_gated:
        fails.append(f"gated keys mismatch: missing={task_gated-rep_gated} extra={rep_gated-task_gated}")

    # 4 deliverables exist + checksums
    rep_deliv = {d["path"]: d["sha256"] for d in report.get("deliverables", [])}
    for d in task_deliv:
        p = report_dir / d
        if not p.exists():
            if not (stopped and d not in ("requirements.lock",) and not d.startswith("scripts/")):
                fails.append(f"deliverable missing on disk: {d}")
            continue
        if d not in rep_deliv:
            fails.append(f"deliverable not listed in report: {d}")
        elif sha256(p) != rep_deliv[d]:
            fails.append(f"checksum mismatch: {d}")
    for d in rep_deliv:
        if d not in task_deliv:
            fails.append(f"report lists deliverable not in task: {d}")

    # 5/6 lockfile
    lock = report_dir / "requirements.lock"
    if not lock.exists():
        fails.append("requirements.lock missing")
    else:
        if sha256(lock) != report.get("lockfile_sha256"):
            fails.append("lockfile_sha256 mismatch")
        if a.prev:
            prev_lock = Path(a.prev).parent / "requirements.lock"
            if prev_lock.exists():
                cur, old = lock_versions(lock), lock_versions(prev_lock)
                for k in sorted(set(cur) | set(old)):
                    if cur.get(k) != old.get(k):
                        warns.append(f"lockfile diff: {k} {old.get(k)} -> {cur.get(k)}")

    # 7 known failures
    expected = parse_known_failures(ROOT / "coord" / "PROJECT.md")
    got = set(report.get("known_failures", []))
    if got != expected:
        fails.append(f"known_failures changed: newly_failing={got-expected} newly_passing={expected-got}")

    # 8 regeneration
    regen = report.get("regeneration", {})
    for g in (rep_gated if stopped else task_gated):
        if g not in regen:
            fails.append(f"no regeneration command for gated `{g}`")

    # 9 observations
    for o in report.get("observations", []):
        if INTERPRETIVE.search(o):
            warns.append(f"observation reads as interpretation: {o[:80]}")

    # status
    if report.get("status") == "blocked" and not report.get("blocked_reason"):
        fails.append("status blocked without blocked_reason")
    if stopped and not report.get("stop_rule"):
        fails.append("status stopped without stop_rule (the decision-rule line that fired)")

    for w in warns:
        print("WARN ", w)
    for f in fails:
        print("FAIL ", f)
    print("RESULT", "FAIL" if fails else "PASS", f"({len(fails)} fail, {len(warns)} warn)")
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
