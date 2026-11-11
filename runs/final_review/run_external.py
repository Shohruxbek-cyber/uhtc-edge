"""External review (9 Nov 2026): run the 2D checks one after the other, detached from the session.

  1. runs/final_review/shoulder_check.py run          (A4 / R1: failure criteria at every surface node)
  2. runs/envelope_2d/run_envelope_2d.py external     (A6 / R2 and A9 / R5: depletion-rate and gate variants)
Progress: runs/final_review/_external.log; runs/final_review/_external.DONE holds the exit codes at the end.
Run: python runs/final_review/run_external.py
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOG = HERE / "_external.log"
env = dict(os.environ, PYTHONUTF8="1", UHTC_WORKERS="4")
steps = [("shoulder", [str(HERE / "shoulder_check.py"), "run"]),
         ("external", [str(ROOT / "runs" / "envelope_2d" / "run_envelope_2d.py"), "external"])]
codes = []
with open(LOG, "a", encoding="utf-8") as log:
    for name, args in steps:
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} start {name}\n")
        log.flush()
        r = subprocess.run([sys.executable] + args, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        codes.append((name, r.returncode))
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} end {name} exit {r.returncode}\n")
        log.flush()
(HERE / "_external.DONE").write_text(f"{codes}\n", encoding="utf-8")
