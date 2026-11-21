from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

import densnv

FILE = Path(__file__).parent / "data" / "DENGUE2_.txt"


@pytest.fixture
def ds():
    return densnv.read_sinave(FILE)


@pytest.fixture
def fake_rgeomex(monkeypatch):
    import rgeomex
    agee = gpd.GeoDataFrame({"CVEGEO": ["31", "01"], "CVE_ENT": ["31", "01"],
                             "NOMGEO": ["Yucatán", "Aguascalientes"]},
                            geometry=[box(-90.5, 19.5, -87.5, 21.6), box(-102.9, 21.6, -101.8, 22.5)],
                            crs=4326)
    agem = gpd.GeoDataFrame({"CVEGEO": ["31050", "31041"], "CVE_ENT": ["31", "31"],
                             "CVE_MUN": ["050", "041"], "NOMGEO": ["Mérida", "Kanasín"]},
                            geometry=[box(-89.8, 20.8, -89.5, 21.2), box(-89.5, 20.8, -89.4, 21.0)],
                            crs=4326)

    def load(name, cve_edo=None):
        g = {"AGEE_inegi19_mx": agee, "AGEM_inegi19_mx": agem}[name].copy()
        if cve_edo is not None:
            keys = cve_edo if isinstance(cve_edo, list) else [cve_edo]
            g = g[g["CVE_ENT"].isin([str(k).zfill(2) for k in keys])]
        return g
    monkeypatch.setattr(rgeomex, "load", load)


def test_plt_heatmap(ds):
    r = densnv.plt_heatmap(ds, state="YUCATÁN", year=2026, status=[1, 2])
    d = r["data"].set_index(["DES_MPO_REP", "week"])["value"]
    assert d[("Mérida", 10)] == 5 and d[("Mérida", 11)] == 4 and d[("Kanasín", 11)] == 0
    assert r["figure"] is not None
    h = densnv.plt_heatmap(ds, state=None, year=2026, hospitalizados=True)
    assert set(h["data"]["DES_MPO_REP"]) == {"Mérida_31", "Kanasín_31"}


def test_plt_pyramid(ds):
    r = densnv.plt_pyramid(ds)
    t = r["data"]
    assert len(t) == 30 and t["count"].sum() == 4          # confirmed cases, all years
    assert t.set_index(["age_class", "sex"]).loc[("30-34", "Female"), "count"] == 4


def test_plt_serotype_and_treemap(ds):
    r = densnv.plt_serotype(ds, year=2026)
    dg = r["diagnosis"].set_index(["area", "DES_DIAG_FINAL"])["n"]
    assert dg[("Yucatán", "DCSA")] == 1 and dg[("Yucatán", "DG")] == 1
    sr = r["serotypes"].set_index(["area", "DENGUE_SER_TRIPLEX"])["n"]
    assert sr[("Yucatán", 2)] == 2 and sr[("Yucatán", 1)] == 1
    m = densnv.plt_serotype(ds, year=2026, cve_edo="31")["serotypes"]
    assert set(m["area"]) == {"Mérida", "Kanasín"}
    t = densnv.plt_treemap(ds, country=False, cve_edo=31, year=2026)["data"]
    assert t.set_index(["area", "DES_DIAG_FINAL"]).loc[("Mérida", "DNG"), "value"] == 8


def test_maps(ds, fake_rgeomex, monkeypatch):
    b = densnv.mp_bumpmap(ds, year=2026)
    assert b["data"].set_index("NOMGEO").loc["Yucatán", "value"] == 11   # all statuses, as in R
    b2 = densnv.mp_bumpmap(ds, year=2026, state=False, cve_edo=31)
    assert set(b2["data"]["NOMGEO"]) == {"Mérida", "Kanasín"}
    s = densnv.mp_serotype(FILE, cve_edo="31")
    assert set(s["serotypes"]["serotype_combination"]) == {"D2D3", "D1"}
    assert s["cases"].set_index("NOMGEO").loc["Mérida", "n"] == 3
    from densnv import _cases
    it = pd.DataFrame({"entidad": ["31"], "nom_ent": ["Yucatán"], "mun": ["050"],
                       "nom_mun": ["Mérida"], "nom_loc": ["Mérida"], "loc": ["0001"],
                       "longitud": [-89.62], "latitud": [20.97], "altitud": [9.0],
                       "pobtot": [921771]})
    monkeypatch.setattr(_cases, "_iter_loc", lambda: it)
    loc = densnv.cases_by_loc(FILE, "DENV")
    m = densnv.mp_cases_loc(loc, var="probable", cve_edo="31")
    assert m["top"]["n"].iloc[0] == 8 and m["map"] is not None
