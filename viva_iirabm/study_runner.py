"""Shared harness for the reproduction studies' canonical runners.

Each study's ``sims/run.py`` calls ``execute(slug, build_fn)``: it records a
run event into the workspace ``.pbg/runs.jsonl``, calls ``build_fn`` to produce
the real data (figures + a result dict), writes the interactive viz and a
``result.json`` into the study directory, and records completion.
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
STUDIES = REPO_ROOT / "workspace" / "studies"
INVESTIGATION = "iirabm"


def _append(event: dict):
    try:
        from vivarium_workbench.lib.run_log import append_run_event
        append_run_event(REPO_ROOT, event)
    except Exception as e:  # run-logging is best-effort
        print(f"[run_log] skipped ({e})")


def execute(slug: str, build_fn):
    study_dir = STUDIES / slug
    viz_dir = study_dir / "viz"
    data_dir = study_dir / "data"
    viz_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    run_id = uuid.uuid4().hex
    spec_id = f"viva_iirabm.study.{slug}"
    _append({"run_id": run_id, "event": "started", "spec_id": spec_id,
             "label": slug, "started_at": time.time(), "status": "running",
             "origin": "canonical_run", "study_slug": slug,
             "investigation_slug": INVESTIGATION, "params": {}})
    t0 = time.time()
    try:
        figures, result = build_fn()
        for name, fig in figures.items():
            out = viz_dir / f"{name}.html"
            fig.write_html(str(out), include_plotlyjs="cdn", full_html=True)
            print(f"  viz -> {out.relative_to(REPO_ROOT)}")
        (data_dir / "result.json").write_text(json.dumps(result, indent=2, default=float))
        print(f"  result -> {(data_dir / 'result.json').relative_to(REPO_ROOT)}")
    except Exception:
        _append({"run_id": run_id, "event": "completed", "completed_at": time.time(),
                 "status": "failed", "n_steps": 0})
        raise
    _append({"run_id": run_id, "event": "completed", "completed_at": time.time(),
             "status": "completed", "n_steps": 1})
    print(f"[{slug}] done in {time.time()-t0:.0f}s  ->  {result.get('headline','')}")
    return result
