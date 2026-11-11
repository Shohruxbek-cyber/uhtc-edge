# NIST-JANAF thermochemical tables (plain-text versions)

**Source:** NIST-JANAF Thermochemical Tables, NIST Standard Reference Database 13 (M. W. Chase Jr., NIST-JANAF Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9, 1998). Web edition: https://janaf.nist.gov/. The site states "Last Update to Data Content: 1998". Dataset DOI 10.18434/T42S31. I checked it at https://api.datacite.org/dois/10.18434/T42S31, which returns the title "NIST-JANAF Thermochemical Tables - SRD 13" (publisher NIST, 2013).
**Status:** NIST reserves copyright on its Standard Reference Data, so the tables are not redistributed with this repository. Run `python literature/janaf/fetch_janaf.py` to download them and check them against the checksums of the tables used for the paper. Dataset DOI: 10.18434/T42S31.
**Retrieved:** 2026-11-02 with `curl` from `https://janaf.nist.gov/tables/<ID>.txt`. Each table ID was looked up in the element index pages `https://janaf.nist.gov/tables/<El>-index.html`.

The files are byte-for-byte as downloaded. Each one keeps the two NIST header lines: line 1 is the name and formula/phase, line 2 holds the column labels. The columns are tab-separated:
`T(K) | Cp | S | -[G-H(Tr)]/T | H-H(Tr) | delta-f H | delta-f G | log Kf`.
Units: J K-1 mol-1 for Cp, S and -[G-H(Tr)]/T; kJ mol-1 for H-H(Tr), delta-f H and delta-f G. Tr = 298.15 K and the standard-state pressure is 0.1 MPa (1 bar).

## Verification

Each file was opened and parsed with `python -I check_janaf.py`, a script kept in the scratch folder outside `literature/`. All 26 files parse. The first data row is T = 0 K and the last is the upper temperature listed below. The "dfH / dfG at 298.15 K" column was read from each file.

| File | NIST table ID | Header line 1 (name, formula/phase) | T range (K) | Data rows | dfH / dfG at 298.15 K (kJ/mol) | Phase transitions marked in table | URL |
|---|---|---|---|---|---|---|---|
| `B2O3_cr.txt` | B-095 | Boron Oxide (B2O3) / B2O3(cr) | 0-1500 | 18 | -1271.936 / -1192.796 | CRYSTAL <--> LIQUID at 723 K | https://janaf.nist.gov/tables/B-095.txt |
| `B2O3_l.txt` | B-096 | Boron Oxide (B2O3) / B2O3(l) | 0-3000 | 33 | -1253.359 / -1181.521 | CRYSTAL <--> LIQUID at 723 K | https://janaf.nist.gov/tables/B-096.txt |
| `B2O3_cr-l.txt` | B-097 | Boron Oxide (B2O3) / B2O3(cr,l) | 0-3000 | 34 | -1271.936 / -1192.796 | CRYSTAL <--> LIQUID at 723 K | https://janaf.nist.gov/tables/B-097.txt |
| `B2O3_g.txt` | B-098 | Boron Oxide (B2O3) / B2O3(g) | 0-6000 | 65 | -835.963 / -825.343 | - | https://janaf.nist.gov/tables/B-098.txt |
| `ZrB2_cr.txt` | B-103 | Zirconium Boride (ZrB2) / B2Zr1(cr) | 0-4000 | 43 | -322.586 / -318.234 | CRYSTAL <--> LIQUID at 3323 K | https://janaf.nist.gov/tables/B-103.txt |
| `ZrB2_l.txt` | B-104 | Zirconium Boride (ZrB2) / B2Zr1(l) | 0-5000 | 55 | -224.832 / -229.067 | GLASS <--> LIQUID at 2200 K; CRYSTAL <--> LIQUID at 3323 K | https://janaf.nist.gov/tables/B-104.txt |
| `ZrB2_cr-l.txt` | B-105 | Zirconium Boride (ZrB2) / B2Zr1(cr,l) | 0-5000 | 54 | -322.586 / -318.234 | CRYSTAL <--> LIQUID at 3323 K | https://janaf.nist.gov/tables/B-105.txt |
| `SiO2_cr_cristobalite-high.txt` | O-035 | Silicon Oxide, Cristobalite, High (SiO2) / O2Si1(cr) | 0-3000 (0 K row is empty) | 32 | -905.489 / -853.636 | - (no melting marked) | https://janaf.nist.gov/tables/O-035.txt |
| `SiO2_cr_cristobalite-low.txt` | O-036 | Silicon Oxide, Cristobalite, Low (SiO2) / O2Si1(cr) | 0-1000 | 12 | -908.346 / -854.509 | - | https://janaf.nist.gov/tables/O-036.txt |
| `SiO2_cr_quartz.txt` | O-037 | Silicon Oxide, Quartz (SiO2) / O2Si1(cr) | 0-2000 | 25 | -910.857 / -856.443 | I <--> II at 847 K; II <--> LIQUID at 1696 K | https://janaf.nist.gov/tables/O-037.txt |
| `SiO2_l.txt` | O-038 | Silicon Oxide (SiO2) / O2Si1(l) | 0-4500 | 48 | -902.661 / -850.175 | II <--> LIQUID at 1696 K | https://janaf.nist.gov/tables/O-038.txt |
| `SiO2_cr-l.txt` | O-039 | Silicon Oxide (SiO2) / O2Si1(cr,l) | 0-4500 | 51 | -910.857 / -856.443 | I <--> II at 847 K; II <--> LIQUID at 1696 K | https://janaf.nist.gov/tables/O-039.txt |
| `SiO2_g.txt` | O-040 | Silicon Oxide (SiO2) / O2Si1(g) | 0-6000 | 65 | -305.432 / -306.926 | - | https://janaf.nist.gov/tables/O-040.txt |
| `SiO_g.txt` | O-012 | Silicon Oxide (SiO) / O1Si1(g) | 0-6000 | 65 | -100.416 / -127.305 | - | https://janaf.nist.gov/tables/O-012.txt |
| `SiC_cr_alpha.txt` | C-100 | Silicon Carbide, Alpha (SiC) / C1Si1(cr) | 0-4000 | 42 | -71.546 / -69.139 | - | https://janaf.nist.gov/tables/C-100.txt |
| `SiC_cr_beta.txt` | C-101 | Silicon Carbide, Beta (SiC) / C1Si1(cr) | 0-4000 | 42 | -73.220 / -70.850 | - | https://janaf.nist.gov/tables/C-101.txt |
| `ZrO2_cr.txt` | O-052 | Zirconium Oxide (ZrO2) / O2Zr1(cr) | 0-3500 | 40 | -1097.463 / -1039.724 | I <--> II at 1478 K; II <--> LIQUID at 2950 K | https://janaf.nist.gov/tables/O-052.txt |
| `ZrO2_l.txt` | O-053 | Zirconium Oxide (ZrO2) / O2Zr1(l) | 0-5000 | 57 | -1023.155 / -972.661 | I <--> II at 1478 K; GLASS <--> LIQUID at 2000 K; II <--> LIQUID at 2950 K | https://janaf.nist.gov/tables/O-053.txt |
| `ZrO2_cr-l.txt` | O-054 | Zirconium Oxide (ZrO2) / O2Zr1(cr,l) | 0-5000 | 56 | -1097.463 / -1039.724 | I <--> II at 1478 K; II <--> LIQUID at 2950 K | https://janaf.nist.gov/tables/O-054.txt |
| `CO_g.txt` | C-093 | Carbon Monoxide (CO) / C1O1(g) | 0-6000 | 62 | -110.527 / -137.163 | - | https://janaf.nist.gov/tables/C-093.txt |
| `CO2_g.txt` | C-095 | Carbon Dioxide (CO2) / C1O2(g) | 0-6000 | 62 | -393.522 / -394.389 | - | https://janaf.nist.gov/tables/C-095.txt |
| `O2_ref.txt` | O-029 | Oxygen (O2) / O2(ref) | 0-6000 | 65 | 0 / 0 | - | https://janaf.nist.gov/tables/O-029.txt |
| `O_g.txt` | O-001 | Oxygen (O) / O1(g) | 0-6000 | 65 | 249.173 / 231.736 | - | https://janaf.nist.gov/tables/O-001.txt |
| `N2_ref.txt` | N-023 | Nitrogen (N2) / N2(ref) | 0-6000 | 65 | 0 / 0 | - | https://janaf.nist.gov/tables/N-023.txt |
| `N_g.txt` | N-002 | Nitrogen (N) / N1(g) | 0-6000 | 65 | 472.683 / 455.540 | - | https://janaf.nist.gov/tables/N-002.txt |
| `NO_g.txt` | N-005 | Nitrogen Oxide (NO) / N1O1(g) | 0-6000 | 65 | 90.291 / 86.600 | - | https://janaf.nist.gov/tables/N-005.txt |

Notes on the choice of tables:
- **O2 and N2:** JANAF lists the "(ref)" tables (O-029, N-023) as the reference-state gas tables. They were saved as `O2_ref.txt` and `N2_ref.txt`; for O2 and N2 these are the gas-phase tables.
- **SiC:** both alpha (C-100) and beta (C-101) tables were saved. No "cr,l" table and no decomposition point is marked in the txt files. SiC decomposition is not a congruent melting point; see INDEX_properties.md.
- **SiO2:** the JANAF cr,l path (O-039) goes through quartz and has quartz II melting (metastably) to liquid at 1696 K. The cristobalite-high table (O-035) runs to 3000 K with no melting marked. Choose the condensed SiO2 phase (cristobalite vs liquid vs glass) deliberately in the model. The open sources give cristobalite melting points of 1713 °C and 1722 °C (PubChem; see `notes/PubChem2026_phase_change_props.md`).
- **ZrO2:** the JANAF cr table (O-052) treats monoclinic as I and tetragonal as II (I <--> II at 1478 K) and melts at 2950 K. No separate cubic phase is tabulated.
- **ZrB2:** the liquid tables (B-104, B-105) are an extrapolated "GLASS <--> LIQUID" construction below the 3323 K melting point.

## Known formatting quirk in the NIST txt files (do not ignore when parsing)

The following rows are reproduced exactly as served by NIST. In them the last three fields are separated by spaces, not tabs, and **the minus signs of delta-f H and delta-f G are missing**. The log Kf value in the same row is positive, which confirms that delta-f G is negative. Example: B2O3(cr) at 800 K gives log Kf = 69.229, so delta-f G = -69.229 x R x 800 x ln10 = -1060.3 kJ/mol. The neighbouring rows are also negative.

| File | Row (as served) |
|---|---|
| B2O3_cr.txt | `800  111.169 140.489 83.017 45.978 1268.779  1060.279 69.229` |
| B2O3_l.txt | `800  129.704 175.883 109.547 53.069 1243.111  1062.926 69.402` |
| SiO2_cr_quartz.txt | `1700 75.981 154.307 97.787 96.083 949.326  609.041 18.714` |
| SiO2_l.txt | `1700 85.772 158.857 102.619 95.604 941.610  609.059 18.714` |
| ZrB2_cr.txt | `3400 96.655 212.387 139.368 248.262 479.212  181.433 2.787` |
| ZrB2_l.txt | `3400 96.232 243.863 168.833 255.105 374.614  183.856 2.825` |
| ZrO2_cr.txt | `3000 74.475 219.844 151.325 205.558 1106.286  543.038 9.455` |
| ZrO2_l.txt | `3000 87.864 249.570 176.587 218.947 1018.589  544.518 9.481` |

A parser should split on any whitespace and should either negate these delta-f values or rebuild delta-f G from log Kf. The files were deliberately **not** edited.

## Quick cross-checks computed from these tables (my calculations, not literature values)

These numbers were computed in the scratch folder with `janaf_checks.py`, using linear interpolation in T, and are given only so the project's own JANAF implementation can be checked against them.
- **B2O3(l) = B2O3(g):** p_B2O3 = exp(-[dfG(g) - dfG(l)]/RT) gives 2.4e-6 bar at 1400 K, 1.47e-4 bar at 1600 K, 3.4e-3 bar at 1800 K, 4.1e-2 bar at 2000 K and 0.31 bar at 2200 K. p = 1 bar is reached at about 2339 K. Published fits are compared in `notes/Hlavac1984_B2O3_vaporization.md`. Near 1600 K JANAF and the fits agree within about 10 %. At 2000 K JANAF is about 2x higher than the fits.
- **SiO2(l) = SiO(g) + 1/2 O2(g):** K = 1.2e-10 bar^1.5 at 1800 K, 2.1e-8 at 2000 K, 1.4e-6 at 2200 K and 4.5e-5 at 2400 K. Over SiO2(l), p(SiO2,g) is 5.5e-9 bar at 1800 K and 2.6e-7 bar at 2000 K.
- **cp(ZrB2, cr):** 48.242 J/(mol K) at 298.15 K (427.5 J/(kg K) with M = 112.846 g/mol) and 82.362 J/(mol K) at 2000 K.
- **cp(SiC, beta):** 26.840 J/(mol K) at 298.15 K and 53.798 J/(mol K) at 2000 K.
