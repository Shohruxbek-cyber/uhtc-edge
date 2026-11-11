# runs/final_review

Runs and calculations of the final review (8–9 Nov 2026) and of the external review response (9 Nov 2026). Each script
writes its outputs next to it; `runs/analysis/trace_numbers.py` (sections `final_review` and `external_review`) checks
every number the paper quotes from them. Run from the repository root with the project Python.

| Script | What it does | Outputs |
|---|---|---|
| `analyses.py` | Final review: rarefaction (Kolodziej Eq. 7), laminar check, Eq. 2 cap, interfacial pressure, time tolerance, converted depth, five-minute boundary, bootstrap intervals of the rank correlations | `analyses.json`, `analyses.md` |
| `checks_1d.py` | Final review: 1D base cases B1 and B2 with the hot-end slope, silica-filled pores, a 50 mm slab, the V2 finest numerics and the O2 budget | `checks_1d.json`, `checks_1d.md` |
| `edge_states.py` | Final review: edge states of the ground tests beyond the 6000 K end of the JANAF tables | `edge_states.json`, `edge_states.md` |
| `ratios_2d.py` | Final review: failure-time ratios of the 2D variants one grid step beyond the baseline boundary | `ratios_2d.md` |
| `v1_case_c_dt.py` | Final review: V1 case C at a 0.25 ms step | `v1_case_c_dt.json` |
| `kd_arrhenius.py` | External review A6 (run R2): Arrhenius fit of the SiC-depletion rate through all of Shugart's depleted-layer depths (1627–1800 °C) | `kd_arrhenius.json` |
| `external_1d.py` | External review A6 and A9 (runs R2, R5): base cases B1 and B2 with the Arrhenius depletion rate (extrapolated, held above 2073 K) and the converted-depth gate of criterion (a) at 0, 1 nm, 0.5 and 2 µm | `external_1d.json`, `external_1d.md` |
| `external_calcs.py` | External review A5 and A11 (runs R4, R6): free-stream and edge Knudsen numbers, Kolodziej rarefaction factor, dynamic pressure at the boundary points and along the X-43A path | `external_calcs.json`, `external_calcs.md` |
| `shoulder_check.py` | External review A4 (run R1): failure criteria at every surface node of the 2D model, for the 15 band columns (highest passing and lowest failing Mach number) and the six 2D X-43A runs, with the stagnation pressure at every node (production) and with the local modified-Newtonian pressure (`run_2d(local_pressure=True)`); criterion (b) margins for the baseline and alternative transition lines | `shoulder_check.jsonl`, `shoulder_check.json`, `shoulder_check.md` |
| `run_external.py` | Detached driver of the external-review 2D runs: `shoulder_check.py run`, then `runs/envelope_2d/run_envelope_2d.py external` (variants `kd_arrh`, `kd_arrh_cap`, `gate_0.5`, `gate_2` on the 15 band columns) | `_external.log`, `_external.DONE`; `runs/envelope_2d/external_*.jsonl/.csv` |
