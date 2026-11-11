# V9: glass-loss clock against the Plasmatron tests of Marschall et al.

Plan, inputs and criterion: `validation/V9_plan.md`, committed before this run. 1D model, no adjusted parameter.

| Bound | Sample | Low plateau (K) [measured] | Peak glass (µm) | Glass at switch (µm) | Glass loss after switch (s) [measured onset] | T at reference time (K) [measured] | Outcome | Timing | Temperature |
|---|---|---|---|---|---|---|---|---|---|
| fc | 1.8 | 2067 [1895] | 93.6 | 68.2 | no loss by 480 s [none by 480 s] | 2387 [2110] | pass | pass | fail |
| fc | 1.5 | 2067 [1825] | 78.8 | 68.0 | no loss by 480 s [68] | 2444 [2210] | fail | fail | fail |
| fc | 1.6 | 2067 [1860] | 71.7 | 68.3 | 71 [21] | 2475 [2250] | pass | fail | fail |
| nc | 1.8 | 1733 [1895] | 3.8 | 3.1 | no loss by 480 s [none by 480 s] | 1872 [2110] | pass | pass | fail |
| nc | 1.5 | 1733 [1825] | 3.8 | 3.1 | no loss by 480 s [68] | 1890 [2210] | fail | fail | fail |
| nc | 1.6 | 1733 [1860] | 26.6 | 3.1 | no loss by 480 s [21] | 1852 [2250] | fail | fail | fail |

V9: **FAIL** (all three samples must pass with the same catalysis bound; fully catalytic False, non-catalytic False).

