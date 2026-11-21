import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

import densnv
from densnv import dbscan_nn, frnn

# a synthetic locality around Merida (no real cases)
LON0, LAT0 = -89.62, 20.97
LOC = gpd.GeoDataFrame({"CVEGEO": ["310500001"], "CVE_ENT": ["31"], "CVE_MUN": ["050"],
                        "CVE_LOC": ["0001"], "NOMGEO": ["Mérida"], "AMBITO": ["Urbana"]},
                       geometry=[box(LON0 - 0.05, LAT0 - 0.05, LON0 + 0.05, LAT0 + 0.05)],
                       crs=4326)


def _blocks():
    step = 0.001   # ~100 m
    cells = [box(LON0 - 0.05 + i * step, LAT0 - 0.05 + j * step,
                 LON0 - 0.05 + (i + 0.9) * step, LAT0 - 0.05 + (j + 0.9) * step)
             for i in range(100) for j in range(100)]
    return gpd.GeoDataFrame({"entidad": 31.0, "manzana": np.arange(len(cells))},
                            geometry=cells, crs=4326)


@pytest.fixture
def fake_rgeomex(monkeypatch):
    import rgeomex
    blocks = _blocks()
    monkeypatch.setattr(rgeomex, "extract_locality", lambda cve_edo, locality: LOC.copy())
    monkeypatch.setattr(rgeomex, "load", lambda name, cve_edo=None: blocks.copy())


def _cases(seed=1):
    rng = np.random.default_rng(seed)
    ref = pd.Timestamp("2026-10-31")
    rows = []
    # background: 300 cases over a year, scattered in the locality
    for _ in range(300):
        rows.append((LON0 + rng.uniform(-0.045, 0.045), LAT0 + rng.uniform(-0.045, 0.045),
                     ref - pd.Timedelta(days=int(rng.integers(0, 365)))))
    # focus 1: 8 recent cases within ~150 m; focus 2: 4 recent cases
    for _ in range(8):
        rows.append((LON0 + rng.normal(0, 0.0007), LAT0 + rng.normal(0, 0.0007),
                     ref - pd.Timedelta(days=int(rng.integers(0, 14)))))
    for _ in range(4):
        rows.append((LON0 + 0.03 + rng.normal(0, 0.0005), LAT0 + 0.03 + rng.normal(0, 0.0005),
                     ref - pd.Timedelta(days=int(rng.integers(0, 14)))))
    d = pd.DataFrame(rows, columns=["long", "lat", "FEC_INI_SIGNOS_SINT"])
    d["ESTATUS_CASO"] = 1
    d.loc[::5, "ESTATUS_CASO"] = 2
    d["CVE_EDO_REP"], d["CVE_MPO_REP"], d["CVE_LOC_REP"] = "31", "050", "0001"
    d["DES_LOC_REP"] = "MÉRIDA"
    return d


def test_dbscan_matches_definition():
    # two groups and one isolated point on a line (metres)
    xy = np.array([[0, 0], [100, 0], [200, 0], [1000, 0], [1100, 0], [1150, 0], [5000, 0]])
    nn = frnn(xy, 150)
    lab = dbscan_nn(nn, 3)
    assert lab.tolist() == [1, 1, 1, 2, 2, 2, 0]
    # border point reached by two clusters joins the first one
    xy = np.array([[0, 0], [10, 0], [20, 0], [30, 0], [115, 0],
                   [200, 0], [210, 0], [220, 0], [230, 0]])
    lab = dbscan_nn(frnn(xy, 90), 4)
    assert lab.tolist() == [1, 1, 1, 1, 1, 2, 2, 2, 2]


def test_dbscan_same_partition_as_sklearn():
    sk = pytest.importorskip("sklearn.cluster")
    rng = np.random.default_rng(3)
    xy = rng.uniform(0, 5000, size=(400, 2))
    lab = dbscan_nn(frnn(xy, 300), 4)
    ref = sk.DBSCAN(eps=300, min_samples=4).fit(xy).labels_
    core = np.array([len(v) + 1 >= 4 for v in frnn(xy, 300)])
    # core points: same partition (labels up to renaming)
    pairs = set(zip(lab[core], ref[core]))
    assert len(pairs) == len(set(lab[core])) == len(set(ref[core]))


def test_mp_active_foci_window(fake_rgeomex):
    r = densnv.mp_active_foci(_cases(), cve_edo="31", locality="Merida", n_sim=199, seed=1,
                              temporal_scale=21, buffer_m=200, print_map=True)
    foci = r["foci"]
    assert len(foci) >= 2
    big = foci.sort_values("casos", ascending=False).iloc[0]
    assert big["casos"] >= 8 and big["manzanas"] > 0
    assert r["blocks"]["foco"].notna().all()
    assert set(r["cases"]["grupo"]) <= {"Caso en foco", "Caso aislado"}
    k = r["knox"]
    assert k["n_casos"] == 312 and len(k["sim"]) == 199 and 0 < k["p_value"] <= 1
    t = r["diagnostico"]["tabla"].iloc[0]
    assert t["casos_periodo"] == 312
    assert t["escenario"].startswith("Transmisión")
    assert r["map"] is not None


def test_mp_active_foci_st_dbscan_and_no_knox(fake_rgeomex):
    r = densnv.mp_active_foci(_cases(), cve_edo=31, locality="Merida", method="st-dbscan",
                              n_sim=0, ref_date="2026-10-31")
    assert r["knox"] is None
    assert set(r["foci"]["tipo"]) <= {"Nuevo", "Persistente"}
    assert r["diagnostico"]["tabla"]["interaccion"].iloc[0] == "Sin prueba"


def test_holm_and_signif():
    from densnv._foci import _holm, _signif
    assert np.allclose(_holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
    assert _signif(0.012345) == 0.012 and _signif(123.4) == 120
