from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

import densnv
from densnv import _cases

DATA = Path(__file__).parent / "data"
FILE = DATA / "DENGUE2_.txt"


def test_read_sinave_types_and_encoding():
    x = densnv.read_sinave(FILE)
    assert len(x) == 14
    assert (x["DES_EDO_RES"] == "OTROS PAISES").sum() == 1
    assert x["CVE_EDO_REP"].iloc[0] == "31"
    assert x["CVE_MPO_REP"].iloc[0] == "050"
    assert x["CVE_LOC_REP"].iloc[0] == "0001"
    assert x["DES_EDO_RES"].iloc[0] == "YUCATÁN"
    assert str(x["ESTATUS_CASO"].dtype) == "Int64"
    assert x["FEC_INI_SIGNOS_SINT"].iloc[0] == pd.Timestamp("2026-03-02")
    assert x.loc[x["ANO"] == 2025, "FEC_INI_SIGNOS_SINT"].iloc[0] == pd.Timestamp("2025-07-15")


def test_read_folder_and_columns():
    full = densnv.read(DATA, vbd="DENGUE", complete=True)
    arbo = densnv.read(DATA, vbd="DENGUE", complete=False)
    assert "OTRA_COLUMNA" in full and "OTRA_COLUMNA" not in arbo
    assert set(arbo.columns) <= set(densnv.ARBO_COLUMNS)
    with pytest.raises(FileNotFoundError):
        densnv.read(DATA, vbd="ZIKA")


def test_read_dengue_dataset_filters():
    c = densnv.read_dengue_dataset(FILE, "country", status_caso=[1, 2])
    assert len(densnv.read_dengue_dataset(FILE, "country", status_caso=3)) == 2
    assert len(c) == 10         # 11 probable/confirmed minus the foreign residence
    s = densnv.read_dengue_dataset(FILE, "state", status_caso=2, des_edo_res="YUCATÁN")
    assert len(s) == 4
    m = densnv.read_dengue_dataset(FILE, "municipality", status_caso=[1, 2],
                                   des_edo_res="YUCATÁN", des_mpo_res="KANASÍN")
    assert len(m) == 1


@pytest.fixture
def fake_iter(monkeypatch):
    it = pd.DataFrame({"entidad": ["31", "31"], "nom_ent": ["Yucatán"] * 2,
                       "mun": ["050", "041"], "nom_mun": ["Mérida", "Kanasín"],
                       "nom_loc": ["Mérida", "Kanasín"], "loc": ["0001", "0001"],
                       "longitud": [-89.62, -89.56], "latitud": [20.97, 20.93],
                       "altitud": [9.0, 10.0], "pobtot": [921771, 77240]})
    monkeypatch.setattr(_cases, "_iter_loc", lambda: it)


def test_cases_by_loc(fake_iter):
    x = densnv.cases_by_loc(FILE, "DENV")
    assert isinstance(x, gpd.GeoDataFrame) and x.crs.to_epsg() == 4326
    m = x.set_index("cvegeo").loc["310500001"]
    # 2026, probable+confirmed, reported in Merida: 4 + 3 + foreign residence = 8
    assert m["probable"] == 8 and m["confirmado"] == 3
    assert (m["DNG"], m["DCSA"], m["DG"]) == (6, 1, 1)
    assert (m["DENV1"], m["DENV2"], m["DENV3"], m["DENV4"]) == (0, 2, 1, 0)
    k = x.set_index("cvegeo").loc["310410001"]
    assert k["probable"] == 1 and k["DENV1"] == 1
    w = densnv.cases_by_loc(FILE, "DENV", week=11)
    assert w.set_index("cvegeo").loc["310500001", "probable"] == 4
    y = densnv.cases_by_loc(FILE, "ZIKV", year=2025)
    assert list(y.columns[:6]) == ["cvegeo", "nom_ent", "nom_mun", "nom_loc",
                                   "probable", "confirmado"]
    assert y["probable"].sum() == 1


def test_cases_by_agem(monkeypatch):
    agem = gpd.GeoDataFrame({"CVEGEO": ["31050", "31041", "01001"],
                             "CVE_ENT": ["31", "31", "01"],
                             "CVE_MUN": ["050", "041", "001"],
                             "NOMGEO": ["Mérida", "Kanasín", "Aguascalientes"]},
                            geometry=[box(0, 0, 1, 1)] * 3, crs=4326)
    monkeypatch.setattr(_cases, "_agem", lambda: agem)
    x = densnv.cases_by_agem(DATA, "DENGUE", status_case=[1, 2])
    assert dict(zip(x["CVE_MUN"], x["n"])) == {"050": 10, "041": 1, "001": 0}
    y = densnv.cases_by_agem(DATA, "DENGUE", status_case=2, country=False, cve_edo="31")
    assert len(y) == 2 and y["n"].sum() == 4


def test_cases_week():
    r = densnv.cases_week(DATA, mpo="MÉRIDA", year=2026, plot=True)
    man = r["manejo"].set_index(["SEM", "MANEJO"])["n"]
    assert man[(11, "Hospitalizado")] == 2
    est = r["estatus"].set_index(["SEM", "status"])["n"]
    # week 10: probable 4, confirmed 1, discarded 1
    # positivity round(1/2*100) = 50 -> estimated 0.5*4 + 1 = 3
    assert est[(10, "probable")] == 4 and est[(10, "estimados")] == 3
    assert r["figure"] is not None


@pytest.mark.online
def test_cases_by_loc_with_rgeomex():
    pytest.importorskip("rgeomex")
    x = densnv.cases_by_loc(FILE, "DENV")
    assert set(x["cvegeo"]) == {"310500001", "310410001"}
    assert x.set_index("cvegeo").loc["310500001", "nom_loc"] == "Mérida"
