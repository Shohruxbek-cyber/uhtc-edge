"""Final review (8 Nov 2026): run the 2D boundary searches approved by the author, one stage after the other.

python runs/envelope_2d/run_final_review.py   -> stage "transport" (6 variants), then stage "active" (adds pa_fe).
python runs/envelope_2d/run_final_review.py kinetics   -> second pass (9 Nov): the kinetics stage, whose only new variants
                                                          are kp_all_x0.33 and kp_all_x3 (R4p2-06).
Stages are restartable (finished cases are skipped). Progress: runs/envelope_2d/_final_review.log; a file
runs/envelope_2d/_final_review.DONE (or _final_review_<stages>.DONE) is written at the end with the exit codes.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOG = HERE / "_final_review.log"
env = dict(os.environ, PYTHONUTF8="1", UHTC_WORKERS="4")
stages = sys.argv[1:] or ["transport", "active"]
codes = []
with open(LOG, "a", encoding="utf-8") as log:
    for stage in stages:
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} start {stage}\n")
        log.flush()
        r = subprocess.run([sys.executable, str(HERE / "run_envelope_2d.py"), stage], cwd=HERE.parents[1], env=env,
                           stdout=log, stderr=subprocess.STDOUT)
        codes.append((stage, r.returncode))
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} end {stage} exit {r.returncode}\n")
        log.flush()
done = "_final_review.DONE" if not sys.argv[1:] else f"_final_review_{'_'.join(stages)}.DONE"
(HERE / done).write_text(f"{codes}\n", encoding="utf-8")
