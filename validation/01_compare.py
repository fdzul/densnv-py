"""Module 1 validation, Python side: compare with the R counts.

Run on the Mac (terminal, with the virtual environment active), from the root
of densnv-py, after validation/01_densnv_r.R.
"""
from pathlib import Path
import sys

import pandas as pd

import densnv

PATH = Path("/Users/fdzul/Library/CloudStorage/Dropbox/dataset/dge/denv/2026/DENGUE2_.txt")
OUT = Path("validation/out")
ok = True


def check(name, py, r, keys):
    global ok
    m = py.merge(r, on=keys, how="outer", suffixes=("_py", "_r"), indicator=True)
    only = m[m["_merge"] != "both"]
    cols = [c[:-3] for c in m.columns if c.endswith("_py")]
    diff = m[(m["_merge"] == "both") & (m[[c + "_py" for c in cols]].values
                                         != m[[c + "_r" for c in cols]].values).any(axis=1)]
    print(f"{name}: {len(m)} rows, {len(only)} only in one side, {len(diff)} with different counts")
    if len(only) or len(diff):
        ok = False
        print(only.head(10).to_string(), diff.head(10).to_string(), sep="\n")


loc_py = densnv.cases_by_loc(PATH, "DENV").drop(columns="geometry")
cols = ["cvegeo", "probable", "confirmado", "DNG", "DCSA", "DG", "DENV1", "DENV2", "DENV3", "DENV4"]
loc_r = pd.read_csv(OUT / "loc_r.csv", dtype={"cvegeo": str})
check("cases_by_loc", loc_py[cols], loc_r[cols], ["cvegeo"])

agem_py = densnv.cases_by_agem(PATH.parent, "DENGUE", status_case=[1, 2])
agem_r = pd.read_csv(OUT / "agem_r.csv", dtype={"CVE_ENT": str, "CVE_MUN": str})
check("cases_by_agem", agem_py[["CVE_ENT", "CVE_MUN", "n"]], agem_r, ["CVE_ENT", "CVE_MUN"])

st = densnv.read_dengue_dataset(PATH, "country", status_caso=[1, 2])
country_py = st["CVE_EDO_RES"].value_counts().rename("Freq").reset_index()
country_r = pd.read_csv(OUT / "country_r.csv", dtype={"CVE_EDO_RES": str})
check("read_dengue_dataset", country_py, country_r, ["CVE_EDO_RES"])
extra = set(country_py["CVE_EDO_RES"]) - set(country_r["CVE_EDO_RES"])
if extra:
    print("DES_EDO_RES of the extra states:",
          st.loc[st["CVE_EDO_RES"].isin(extra), "DES_EDO_RES"].map(repr).value_counts().to_dict())

print("\nOK: same results as R" if ok else "\nThere are differences (see above)")
sys.exit(0 if ok else 1)
