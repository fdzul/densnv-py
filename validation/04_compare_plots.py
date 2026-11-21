"""Module 4 validation (graphs and maps), Python side: compare the tables
behind the figures with R. Run after validation/04_plots_r.R."""
from pathlib import Path
import contextlib
import io
import sys

import numpy as np
import pandas as pd

import densnv

PATH = Path("/Users/fdzul/Library/CloudStorage/Dropbox/dataset/dge/denv/2026/DENGUE2_.txt")
OUT = Path("validation/out")
YEAR = 2026
ds = densnv.read_sinave(PATH)
ok = True


def same(name, py, r, keys, val, tol=0.0):
    global ok
    py = py.copy(); r = r.copy()
    for c in keys:
        py[c] = py[c].astype(str); r[c] = r[c].astype(str)
    m = py[keys + [val]].merge(r[keys + [val]], on=keys, how="outer",
                               suffixes=("_py", "_r"), indicator=True)
    only = m[m["_merge"] != "both"]
    both = m[m["_merge"] == "both"]
    diff = both[~np.isclose(both[val + "_py"].astype(float), both[val + "_r"].astype(float),
                            atol=tol)]
    good = only.empty and diff.empty
    ok &= good
    print(f"{name:22s} {len(m):5d} rows  {'OK' if good else 'DIFFERENT'}")
    if not good:
        print(only.head(8).to_string(), diff.head(8).to_string(), sep="\n")


rd = lambda f: pd.read_csv(OUT / f, encoding="utf-8", keep_default_na=False)  # noqa: E731

h = densnv.plt_heatmap(ds, state="YUCATAN", year=YEAR, status=[1, 2], plot=False)["data"]
same("plt_heatmap state", h, rd("heatmap_state_r.csv"), ["DES_MPO_REP", "week"], "value")
h = densnv.plt_heatmap(ds, state=None, year=YEAR, status=2, hospitalizados=True,
                       plot=False)["data"]
same("plt_heatmap country", h, rd("heatmap_country_r.csv"), ["DES_MPO_REP", "week"], "value")

p = densnv.plt_pyramid(ds, plot=False)["data"]
same("plt_pyramid", p, rd("pyramid_r.csv"), ["age_class", "sex"], "count")

t = densnv.plt_treemap(ds, country=True, year=YEAR, plot=False)["data"]
r = rd("treemap_country_r.csv").rename(columns={"DES_EDO_REP": "area"})
same("plt_treemap country", t, r, ["area", "DES_DIAG_FINAL"], "value")
t = densnv.plt_treemap(ds, country=False, cve_edo="31", year=YEAR, plot=False)["data"]
r = rd("treemap_state_r.csv").rename(columns={"DES_MPO_REP": "area"})
same("plt_treemap state", t, r, ["area", "DES_DIAG_FINAL"], "value")

b = densnv.mp_bumpmap(ds, year=YEAR, plot=False)["data"]
same("mp_bumpmap country", b, rd("bump_country_r.csv"), ["NOMGEO"], "value")
same("mp_bumpmap ranking", b, rd("bump_country_r.csv"), ["NOMGEO"], "fine_cap", tol=1e-6)
b = densnv.mp_bumpmap(ds, year=YEAR, state=False, cve_edo="31", plot=False)["data"]
same("mp_bumpmap state", b, rd("bump_state_r.csv"), ["NOMGEO"], "value")

s = densnv.plt_serotype(ds, year=YEAR, plot=False)
r = rd("serotype_diag_r.csv")
r["area"] = densnv._plots._title(r["DES_EDO_REP"])
r["DES_DIAG_FINAL"] = r["DES_DIAG_FINAL"].map({"DENGUE CON SIGNOS DE ALARMA": "DCSA",
                                               "DENGUE GRAVE": "DG", "DENGUE NO GRAVE": "DNG"})
same("plt_serotype diag", s["diagnosis"], r, ["area", "DES_DIAG_FINAL"], "n")
r = rd("serotype_ser_r.csv")
r["area"] = densnv._plots._title(r["DES_EDO_REP"])
same("plt_serotype serotype", s["serotypes"], r, ["area", "DENGUE_SER_TRIPLEX"], "n")

with contextlib.redirect_stdout(io.StringIO()):
    ms = densnv.mp_serotype(PATH, cve_edo="31", plot=False)["serotypes"]
py = ms.melt(id_vars=["CVE_MUN"], value_vars=["D1", "D2", "D3", "D4"],
             var_name="DENGUE_SER_TRIPLEX", value_name="n")
py = py[py["n"] > 0]
py["DENGUE_SER_TRIPLEX"] = py["DENGUE_SER_TRIPLEX"].str[1]
r = rd("mp_serotype_r.csv").rename(columns={"CVE_MPO_REP": "CVE_MUN"})
r["CVE_MUN"] = r["CVE_MUN"].astype(str).str.zfill(3)
same("mp_serotype", py, r, ["CVE_MUN", "DENGUE_SER_TRIPLEX"], "n")

print("\nOK: same tables as R" if ok else "\nThere are differences (see above)")
sys.exit(0 if ok else 1)
