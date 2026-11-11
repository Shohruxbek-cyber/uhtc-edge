"""ZrB2-SiC composite properties for the conduction solvers (sources: uhtc.data and INPUTS.md section 6).

  rho(T)  rule-of-mixtures room-temperature density (ZrB2 6090, SiC 3210 kg/m3) scaled with the measured thermal
          expansion of ZrB2-20SiC (TPSX: 5560 -> 5310 kg/m3, 293 -> 2273 K);
  k(T)    TPSX ManLabs table for ZrB2-20SiC (373-2273 K, held constant outside), times a composition factor
          (1.17 for 10 vol%, 0.75 for 30 vol%; UNSURE) and an optional sensitivity multiplier;
  cp(T)   mass-weighted rule of mixtures of the JANAF heat capacities of ZrB2(cr) and beta-SiC (agrees with the
          TPSX drop-calorimetry data for 20 vol% within about 3 %, INDEX_properties 2.2).
"""
from __future__ import annotations

import numpy as np

from . import data as dt
from . import thermo as th
from .conduction1d import Material


def mass_fractions(f_SiC):
    mZ = (1 - f_SiC) * dt.RHO_ZrB2
    mS = f_SiC * dt.RHO_SiC
    return mZ / (mZ + mS), mS / (mZ + mS)


def rho_rt(f_SiC):
    return (1 - f_SiC) * dt.RHO_ZrB2 + f_SiC * dt.RHO_SiC


def composite(f_SiC=0.2, k_mult=1.0, cp_mult=1.0, T_max=4000.0, name=None) -> Material:
    if f_SiC not in dt.K_SCALE_SiC:
        raise ValueError("f_SiC must be 0.1, 0.15, 0.2 or 0.3 (composition scaling of k is defined for these only)")
    kf = dt.K_SCALE_SiC[f_SiC] * k_mult
    r0 = rho_rt(f_SiC)
    wZ, wS = mass_fractions(f_SiC)
    tZ, tS = th.table("ZrB2"), th.table("SiC")
    MZ, MS = th.MOLAR_MASS["ZrB2"], th.MOLAR_MASS["SiC"]

    def rho(T):
        return r0 * np.interp(T, dt.TPSX_RHO_T, dt.TPSX_RHO) / dt.TPSX_RHO[0]

    def k(T):
        return kf * np.interp(T, dt.TPSX_K_T, dt.TPSX_K)

    def cp(T):
        Tc = np.clip(np.asarray(T, float), 100.0, min(tS.T[-1], tZ.T[-1]))
        return cp_mult * (wZ * np.interp(Tc, tZ.T, tZ.Cp) / MZ + wS * np.interp(Tc, tS.T, tS.Cp) / MS)

    return Material(name or f"ZrB2-{int(round(f_SiC * 100))}SiC", rho=rho, k=k, cp=cp, T_min=150.0, T_max=T_max,
                    sources={"k": "TPSX id 28 (ManLabs) x composition factor", "cp": "JANAF ZrB2 + beta-SiC",
                             "rho": "rule of mixtures x TPSX expansion"})
