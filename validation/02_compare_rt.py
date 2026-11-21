"""Module 2 validation (Rt), Python side: compare with the R estimates.

Run on the Mac (terminal, virtual environment active), from the root of
densnv-py, after validation/02_rt_r.R.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

import densnv

PATH = Path("/Users/fdzul/Library/CloudStorage/Dropbox/dataset/dge/denv/2026/DENGUE2_.txt")
OUT = Path("validation/out")
areas = {
    "mx": dict(level="country"),
    "yuc": dict(level="state", cve_edo="31"),
    "merida": dict(level="municipality", cve_edo="31", cve_mpo="050"),
    "merida_loc": dict(level="locality", cve_edo="31", cve_mpo="050", cve_loc="0001"),
    "acapulco": dict(level="locality", cve_edo="12", cve_mpo="001", cve_loc="0001"),
    "yuc_res": dict(level="state", cve_edo="31", place="RES"),
}
dataset = densnv.read_sinave(PATH)
win_r = pd.read_csv(OUT / "rt_window_r.csv").set_index("area")
curve_r = pd.read_csv(OUT / "rt_curve_r.csv")
ok = True
for k, a in areas.items():
    r = densnv.rt_dengue(dataset, plot=False, window_days=28, **a)
    w = r["rt_ventana"].iloc[0]
    wr = win_r.loc[k]
    same_win = (w["casos"] == wr["casos"] and w["fase"] == wr["fase"]
                and np.allclose([w["rt"], w["rt_inf"], w["rt_sup"]],
                                [wr["rt"], wr["rt_inf"], wr["rt_sup"]], rtol=1e-6, equal_nan=True))
    cr = curve_r[curve_r["area"] == k].set_index("t_end")
    cp = r["rt"].set_index("t_end") if r["rt"] is not None else pd.DataFrame(columns=["mean", "q025", "q975"])
    same_curve = (len(cr) == len(cp) and cr.index.equals(cp.index)
                  and np.allclose(cr[["mean", "q025", "q975"]], cp[["mean", "q025", "q975"]], rtol=1e-6))
    print(f"{k:11s} window: {'OK' if same_win else 'DIFFERENT'}   curve ({len(cp)} windows): "
          f"{'OK' if same_curve else 'DIFFERENT'}")
    if not (same_win and same_curve):
        ok = False
        print("   R :", wr.to_dict())
        print("   Py:", w.to_dict(), "| R curve windows:", len(cr))
print("\nOK: same Rt as R" if ok else "\nThere are differences (see above)")
sys.exit(0 if ok else 1)
