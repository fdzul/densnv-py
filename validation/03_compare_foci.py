"""Module 3 validation (mp_active_foci), Python side: compare with R.

Run on the Mac (terminal, virtual environment active), from the root of
densnv-py, after validation/03_foci_r.R. The deterministic results (foci,
blocks, areas, Knox pairs, coverage, Rt, scenario) must be equal; the Knox
p-values come from random permutations, so they are only compared roughly.
"""
from pathlib import Path
import contextlib
import io
import sys

import numpy as np
import pandas as pd

import densnv

OUT = Path("validation/out")
z = pd.read_csv(OUT / "z_geocoded.csv", dtype=str, encoding="utf-8")
for c in ("long", "lat"):
    z[c] = pd.to_numeric(z[c])
z["ESTATUS_CASO"] = pd.to_numeric(z["ESTATUS_CASO"])
foci_r = pd.read_csv(OUT / "foci_r.csv")
diag_r = pd.read_csv(OUT / "diag_r.csv").set_index("run")

ok = True
for run in diag_r.index:
    edo, loc, method = run.split("|")
    with contextlib.redirect_stdout(io.StringIO()):
        r = densnv.mp_active_foci(z, cve_edo=edo, locality=loc, temporal_scale=28,
                                  buffer_m=0, method=method, n_sim=999, seed=123)
    t = r["diagnostico"]["tabla"].iloc[0]
    k = r["knox"]
    dr = diag_r.loc[run]
    checks = {
        "casos_periodo": t["casos_periodo"] == dr["casos_periodo"],
        "casos_ventana": t["casos_ventana"] == dr["casos_ventana"],
        "knox pares": k["knox"] == dr["knox"] and k["pares_espacio"] == dr["pares_espacio"],
        "focos": t["focos"] == dr["focos"],
        "manzanas": t["manzanas"] == dr["manzanas"],
        "cobertura": abs(t["cobertura_casos"] - dr["cobertura_casos"]) <= 0.1,
        "% en focos": t["pct_casos_en_focos"] == dr["pct_casos_en_focos"],
        "% área focos": abs(t["pct_area_focos"] - dr["pct_area_focos"]) <= 0.1,
        "rt": (pd.isna(t["rt"]) and pd.isna(dr["rt"])) or t["rt"] == dr["rt"],
        "fase": t["fase"] == dr["fase"],
        "escenario": t["escenario"] == dr["escenario"],
    }
    fr = foci_r[foci_r["run"] == run]
    fp = r["foci"]
    key = ["casos", "casos_recientes", "confirmados", "knox_pares", "manzanas"]
    if fp is None:
        checks["foci"] = fr.empty
    else:
        a = fp[key].astype(int).sort_values(key).to_numpy()
        b = fr[key].astype(int).sort_values(key).to_numpy()
        checks["foci"] = a.shape == b.shape and (a == b).all()
    bad = [c for c, v in checks.items() if not v]
    print(f"{run:40s} {'OK' if not bad else 'DIFFERENT: ' + ', '.join(bad)}"
          f"   (Knox p R {dr['knox_p']:.3f} / Py {k['p_value']:.3f})")
    if bad:
        ok = False
        print("   R :", dr[["casos_periodo", "casos_ventana", "knox", "focos", "manzanas",
                          "cobertura_casos", "pct_casos_en_focos", "pct_area_focos", "rt",
                          "fase", "escenario"]].to_dict())
        print("   Py:", {**t[["casos_periodo", "casos_ventana", "focos", "manzanas",
                                "cobertura_casos", "pct_casos_en_focos", "pct_area_focos",
                                "rt", "fase", "escenario"]].to_dict(), "knox": k["knox"]})
print("\nOK: same foci and diagnosis as R" if ok else "\nThere are differences (see above)")
sys.exit(0 if ok else 1)
