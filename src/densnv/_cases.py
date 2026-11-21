"""Aggregate the SINAVE cases by municipality, locality and week."""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd

from ._read import _as_list, read, read_sinave

_SEROTYPES = {1: "DENV1", 2: "DENV2", 3: "DENV3", 4: "DENV4"}
_DIAGNOSES = {"DNG": "DENGUE NO GRAVE", "DCSA": "DENGUE CON SIGNOS DE ALARMA",
              "DG": "DENGUE GRAVE"}


def _agem() -> gpd.GeoDataFrame:
    import rgeomex
    return rgeomex.load("AGEM_inegi19_mx")


def _iter_loc() -> pd.DataFrame:
    import rgeomex
    return rgeomex.load("iter_2020_loc")


def cases_by_agem(path_vbd: str | Path, vbd: str, status_case=(1, 2),
                  country: bool = True, cve_edo=None) -> gpd.GeoDataFrame:
    """Cases by municipality (AGEM) of report.

    Parameters
    ----------
    path_vbd : str or Path
        Folder with the SINAVE txt files.
    vbd : str
        Text that identifies the files of the disease, e.g. ``"DENGUE"``.
    status_case : int or list of int
        1 probable, 2 confirmed, 3 discarded.
    country : bool
        ``True`` returns all the municipalities of Mexico; ``False`` only those
        of ``cve_edo``.
    cve_edo : str or list of str, optional
        INEGI state key(s), used when ``country=False``.

    Returns
    -------
    geopandas.GeoDataFrame
        The municipalities (INEGI 2019) with the number of cases ``n``
        (0 where there are no cases).
    """
    x = read(path_vbd, vbd=vbd, complete=True)
    x = x[x["ESTATUS_CASO"].isin(_as_list(status_case))]
    n = (x.groupby(["CVE_EDO_REP", "CVE_MPO_REP"]).size()
          .rename("n").reset_index())
    xy = _agem().merge(n, how="left", left_on=["CVE_ENT", "CVE_MUN"],
                       right_on=["CVE_EDO_REP", "CVE_MPO_REP"])
    xy = xy.drop(columns=["CVE_EDO_REP", "CVE_MPO_REP"])
    xy["n"] = xy["n"].fillna(0).astype(int)
    if not country:
        keys = [str(k).zfill(2) for k in _as_list(cve_edo)]
        xy = xy[xy["CVE_ENT"].isin(keys)]
    return xy.reset_index(drop=True)


def cases_by_loc(path: str | Path, snv_dataset: str = "DENV", year=None,
                 week=None) -> gpd.GeoDataFrame:
    """Arbovirus cases by locality of report.

    Counts the probable and confirmed cases by locality of report and joins
    them to the INEGI 2020 locality catalogue (ITER 2020): one point per
    locality with cases.

    Parameters
    ----------
    path : str or Path
        The SINAVE txt file (``"DENGUE2_.txt"``, ``"ZIKA2_.txt"`` or
        ``"CHIKUNGUNYA2_.txt"``).
    snv_dataset : {"DENV", "ZIKV", "CHIKV"}
    year : int or list of int, optional
        Year(s) to keep (``ANO``). The most recent year by default.
    week : int or list of int, optional
        Epidemiological weeks to keep (``SEM``). All by default.

    Returns
    -------
    geopandas.GeoDataFrame
        Points (EPSG:4326) with ``cvegeo``, ``nom_ent``, ``nom_mun``,
        ``nom_loc``, ``probable`` (probable + confirmed), ``confirmado``,
        ``altitud`` and ``pobtot``; for DENV also ``DNG``, ``DCSA``, ``DG``
        and ``DENV1`` to ``DENV4``.
    """
    if snv_dataset not in ("DENV", "ZIKV", "CHIKV"):
        raise ValueError("snv_dataset must be 'DENV', 'ZIKV' or 'CHIKV'")
    x = read_sinave(path)
    x = x[x["ESTATUS_CASO"].isin([1, 2])]
    x = x[(x["CVE_EDO_REP"] != "00") & (x["CVE_MPO_REP"] != "000")
          & (x["CVE_LOC_REP"] != "0000")]
    x = x.assign(cvegeo=x["CVE_EDO_REP"] + x["CVE_MPO_REP"] + x["CVE_LOC_REP"])
    if year is None:
        year = int(x["ANO"].max())
    x = x[x["ANO"].isin(_as_list(year))]
    if week is not None:
        x = x[x["SEM"].isin(_as_list(week))]
    if x.empty:
        raise ValueError(f"there are no probable or confirmed cases for year {year}"
                         + (f" and weeks {week}" if week is not None else ""))

    g = x.groupby("cvegeo")
    data = pd.DataFrame({"probable": g.size(),
                         "confirmado": g["ESTATUS_CASO"].apply(lambda s: int((s == 2).sum()))})
    if snv_dataset == "DENV":
        for k, v in _DIAGNOSES.items():
            data[k] = g["DES_DIAG_FINAL"].apply(lambda s, v=v: int((s == v).sum()))
        for k, v in _SEROTYPES.items():
            data[v] = g["DENGUE_SER_TRIPLEX"].apply(lambda s, k=k: int((s == k).sum()))
    data = data.reset_index()

    loc = _iter_loc()
    loc = loc.assign(cvegeo=loc["entidad"] + loc["mun"] + loc["loc"])
    loc = loc[["cvegeo", "nom_ent", "nom_mun", "nom_loc", "longitud", "latitud",
               "altitud", "pobtot"]]
    out = data.merge(loc, how="left", on="cvegeo")
    out = out.dropna(subset=["longitud", "latitud"]).reset_index(drop=True)
    geom = gpd.points_from_xy(out["longitud"], out["latitud"], crs=4326)
    cols = ["cvegeo", "nom_ent", "nom_mun", "nom_loc"] + \
        [c for c in data.columns if c != "cvegeo"] + ["altitud", "pobtot"]
    return gpd.GeoDataFrame(out[cols].reset_index(drop=True), geometry=geom,
                            crs=4326)


def cases_week(path: str | Path, mpo, year: int, plot: bool = True) -> dict:
    """Dengue cases by epidemiological week in a municipality.

    Parameters
    ----------
    path : str or Path
        Folder with the SINAVE txt files (the dengue file has ``"DEN"`` in
        its name).
    mpo : str or list of str
        Municipality(ies) of report in uppercase, as in SINAVE
        (``DES_MPO_REP``).
    year : int
        Year (``ANO``).
    plot : bool
        Also draw the figure (matplotlib).

    Returns
    -------
    dict
        ``"manejo"``: probable + confirmed cases by week and management
        (hospitalised / outpatient); ``"estatus"``: probable, confirmed and
        estimated cases by week (estimated = positivity × probable +
        confirmed); ``"figure"``: the matplotlib figure (``None`` when
        ``plot=False``).
    """
    x = read(path, vbd="DEN", complete=True)
    x = x[(x["ANO"] == year) & x["DES_MPO_REP"].isin(_as_list(mpo))]

    a = x[x["ESTATUS_CASO"].isin([1, 2])]
    manejo = (a.assign(MANEJO=a["MANEJO"].map(lambda m: "Hospitalizado" if m == 1 else "Ambulatorio"))
               .groupby(["SEM", "MANEJO"]).size().rename("n").reset_index())

    b = (x.groupby(["DES_MPO_REP", "SEM", "ESTATUS_CASO"]).size()
          .unstack("ESTATUS_CASO", fill_value=0))
    for k in (1, 2, 3):
        if k not in b:
            b[k] = 0
    b = b.rename(columns={1: "probable", 2: "confirmado", 3: "descartado"})
    pos = (b["confirmado"] / (b["confirmado"] + b["descartado"]) * 100).round(0)
    b["estimados"] = pos / 100 * b["probable"] + b["confirmado"]
    estatus = (b[["probable", "confirmado", "estimados"]].reset_index()
                .melt(id_vars=["DES_MPO_REP", "SEM"], var_name="status", value_name="n"))
    estatus.columns.name = None

    fig = _plot_cases_week(manejo, estatus) if plot else None
    return {"manejo": manejo, "estatus": estatus, "figure": fig}


def _plot_cases_week(manejo: pd.DataFrame, estatus: pd.DataFrame):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for (name, d), col in zip(manejo.groupby("MANEJO"), ("#D7DF71", "#0E5F76")):
        ax[0].plot(d["SEM"], d["n"], color=col, lw=2, marker="o", ms=8,
                   mec="white", label=name)
    cols = {"confirmado": "#E01E5A", "estimados": "#2EB67D", "probable": "#ECB32D"}
    for name, d in estatus.groupby("status"):
        ax[1].plot(d["SEM"], d["n"], color=cols[name], lw=2, marker="o", ms=8,
                   mec="white", label=name)
    for a in ax:
        a.set_ylabel("Número de Casos")
        a.legend(frameon=False, loc="upper left")
        a.spines[["top", "right"]].set_visible(False)
    ax[1].set_xlabel("Semanas Epidemiológicas")
    ax[1].set_xticks(range(0, 53, 2))
    fig.tight_layout()
    return fig
