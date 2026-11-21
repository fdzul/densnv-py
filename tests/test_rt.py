import numpy as np
import pandas as pd
import pytest

import densnv


def test_discr_si_is_a_distribution():
    w = densnv.discr_si(np.arange(200), 21, 6.1)
    assert w[0] == 0
    assert w.sum() == pytest.approx(1, abs=1e-6)
    assert (w * np.arange(200)).sum() == pytest.approx(21, abs=1e-3)


def _dataset(daily, start="2025-01-01", edo="31", mpo="050", loc="0001"):
    dates = pd.date_range(start, periods=len(daily), freq="D")
    onset = np.repeat(dates, daily)
    return pd.DataFrame({"FEC_INI_SIGNOS_SINT": onset, "ESTATUS_CASO": 1,
                         "CVE_EDO_REP": edo, "CVE_MPO_REP": mpo, "CVE_LOC_REP": loc,
                         "DES_EDO_REP": "YUCATÁN", "DES_MPO_REP": "MÉRIDA",
                         "DES_LOC_REP": "MÉRIDA"})


def test_constant_incidence_gives_rt_one():
    r = densnv.rt_dengue(_dataset([20] * 400), level="locality", cve_edo="31",
                         cve_mpo=50, cve_loc=1, plot=False)
    w = r["rt_ventana"].iloc[0]
    assert w["rt"] == pytest.approx(1, abs=0.05)
    assert w["fase"] == "Estable"
    assert r["area"] == "MÉRIDA"
    # the curve starts once there is a full history of cases inside the period
    assert r["rt"]["mean"].between(0.95, 1.05).all()


def test_growth_and_decline():
    t = np.arange(400)
    up = densnv.rt_dengue(_dataset(np.round(5 * np.exp(0.01 * t)).astype(int)),
                          level="state", cve_edo=31, plot=False)
    assert up["rt_ventana"]["fase"].iloc[0] == "En ascenso"
    down = densnv.rt_dengue(_dataset(np.round(300 * np.exp(-0.01 * t)).astype(int)),
                            level="country", plot=True)
    assert down["rt_ventana"]["fase"].iloc[0] == "En descenso"
    assert down["figure"] is not None


def test_no_history_and_errors():
    # the window starts less than mean_si + 2 std_si days after the first onset
    r = densnv.rt_dengue(_dataset([20] * 30), plot=False)
    assert r["rt_ventana"]["fase"].iloc[0] == "Sin estimar"
    assert r["rt"] is None
    with pytest.raises(ValueError):
        densnv.rt_dengue(_dataset([1]), level="state")
    with pytest.raises(ValueError):
        densnv.rt_dengue(_dataset([1]), level="locality", cve_edo=31, cve_mpo=[50, 41], cve_loc=1)
