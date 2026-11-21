"""Graphs and maps of the SINAVE dengue cases (module 4).

Each function returns a dict with the data behind the figure (the same table
that densnv builds in R, so it can be checked) and the figure itself.
"""
from __future__ import annotations

import re
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd

from ._read import _FOREIGN, _as_list

_SUFFIX = re.compile(r" Quint Roo| Yuc| Ver| Mich| Nvo Leon| Tab")
_DIAG = ["DENGUE CON SIGNOS DE ALARMA", "DENGUE NO GRAVE", "DENGUE GRAVE"]
_SEROTYPE_COLORS = {1: "#FAB255", 2: "#43B284", 3: "#0F7BA2", 4: "#DD5129"}


def _title(s: pd.Series) -> pd.Series:
    """Like stringr::str_to_title() (ICU word breaks): the first letter of each
    word in uppercase, where letters joined by ".", "'" or "_" are one word
    ("D.F." -> "D.f.", "O'HIGGINS" -> "O'higgins")."""
    return s.astype(str).str.replace(
        r"[\w.'’]+", lambda m: m.group(0)[:1].upper() + m.group(0)[1:].lower(),
        regex=True)


def _int(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").astype("Int64")


# --------------------------------------------------------------------------
# plt_heatmap
# --------------------------------------------------------------------------
def plt_heatmap(dataset: pd.DataFrame, state=None, year: int | None = None,
                status=(1, 2), hospitalizados: bool = False, plot: bool = True) -> dict:
    """Heatmap of the cases by municipality (rows) and epidemiological week.

    Parameters
    ----------
    dataset : pandas.DataFrame
        SINAVE dengue dataset (``read_sinave()``).
    state : str or list, optional
        State(s) of report, as in ``DES_EDO_REP`` (e.g. ``"GUERRERO"``).
        ``None`` uses all the municipalities of the country.
    year : int
        Year (``ANO``).
    status : int or list of int
        Case status (ignored with ``hospitalizados=True``).
    hospitalizados : bool
        Use the hospitalised cases (with date of admission) instead of the
        status.
    plot : bool
        Also build the plotly figure.

    Returns
    -------
    dict
        ``data`` (municipality, week, value) and ``figure`` (plotly).
    """
    x = dataset[_int(dataset["ANO"]) == year].copy()
    # as in R, cases without the name of the municipality are kept (empty name)
    x["DES_MPO_REP"] = x["DES_MPO_REP"].fillna("")
    if state is not None:
        x = x[x["DES_EDO_REP"].isin(_as_list(state))]
    else:
        x["DES_MPO_REP"] = x["DES_MPO_REP"].astype(str) + "_" + \
            _int(x["CVE_EDO_REP"]).astype(str)
    if hospitalizados:
        x = x[pd.to_datetime(x["FEC_INGRESO"], errors="coerce").notna()]
    else:
        x = x[_int(x["ESTATUS_CASO"]).isin(_as_list(status))]

    order = x.groupby("DES_MPO_REP").size().sort_values(kind="stable")
    levels = _title(pd.Series(order.index)).str.replace(_SUFFIX, "", regex=True)
    tab = (x.groupby(["DES_MPO_REP", "SEM"]).size().unstack("SEM", fill_value=0)
             .stack().rename("value").reset_index().rename(columns={"SEM": "week"}))
    tab["DES_MPO_REP"] = _title(tab["DES_MPO_REP"]).str.replace(_SUFFIX, "", regex=True)
    tab["week"] = tab["week"].astype(int)
    tab["DES_MPO_REP"] = pd.Categorical(tab["DES_MPO_REP"],
                                        categories=list(dict.fromkeys(levels)))
    tab = tab.sort_values(["value", "DES_MPO_REP"], ascending=False).reset_index(drop=True)

    fig = None
    if plot:
        try:
            import plotly.graph_objects as go
        except ImportError:
            warnings.warn("install plotly to draw the heatmap", stacklevel=2)
        else:
            w = tab.pivot_table(index="DES_MPO_REP", columns="week", values="value",
                                observed=False, fill_value=0)
            fig = go.Figure(go.Heatmap(z=w.to_numpy(), x=w.columns, y=w.index.astype(str),
                                       colorscale="Viridis"))
            title = "Casos Hospitalizados" if hospitalizados else "Casos Confirmados de Dengue"
            fig.update_layout(
                title={"text": f"<b>{title}</b>", "font": {"color": "#E01E5A", "size": 14}},
                xaxis={"title": "<b>Semanas Epidemiológicas</b>", "color": "#E01E5A",
                       "tickfont": {"color": "#006AFF"},
                       "rangeslider": {"visible": state is not None}},
                yaxis={"title": "<b>Municipios</b>", "color": "#E01E5A",
                       "tickfont": {"color": "#006AFF"}})
    return {"data": tab, "figure": fig}


# --------------------------------------------------------------------------
# plt_pyramid
# --------------------------------------------------------------------------
_AGE_LABELS = ["<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
               "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65+"]
_AGE_BREAKS = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65]


def plt_pyramid(data: pd.DataFrame, col_male: str = "#0072B2",
                col_female: str = "#CC79A7", gap: float = 0.18,
                show_zero: bool = False, title: str | None = None,
                subtitle: str | None = None,
                caption: str = "Fuente: SINAVE/DGE, Secretaría de Salud.",
                plot: bool = True) -> dict:
    """Age-sex pyramid of the confirmed cases (``ESTATUS_CASO == 2``).

    Males on the left, females on the right, by five-year age group.

    Returns
    -------
    dict
        ``data`` (age_class, sex, count; groups without cases with 0) and
        ``figure`` (matplotlib).
    """
    d = data.copy()
    d.columns = [c.lower() for c in d.columns]
    d = d[_int(d["estatus_caso"]) == 2]
    sex = _int(d["ide_sex"]).map({1: "Female", 2: "Male"})
    age = pd.to_numeric(d["ide_eda_ano"], errors="coerce")
    cls = pd.cut(age, bins=_AGE_BREAKS + [np.inf], labels=_AGE_LABELS, right=False)
    tab = pd.DataFrame({"age_class": cls, "sex": sex}).dropna()
    tab = (tab.groupby(["age_class", "sex"], observed=False).size()
              .rename("count").reset_index())
    full = pd.MultiIndex.from_product([_AGE_LABELS, ["Female", "Male"]],
                                      names=["age_class", "sex"])
    tab = (tab.set_index(["age_class", "sex"])["count"].reindex(full, fill_value=0)
              .reset_index())
    tab["age_class"] = pd.Categorical(tab["age_class"], categories=_AGE_LABELS, ordered=True)

    fig = None
    if plot:
        import matplotlib.pyplot as plt
        max_n = max(int(tab["count"].max()), 1)
        g = max_n * gap
        fig, ax = plt.subplots(figsize=(8, 7))
        for i, lab in enumerate(_AGE_LABELS, start=1):
            for sx, col in (("Male", col_male), ("Female", col_female)):
                n = int(tab.loc[(tab["age_class"] == lab) & (tab["sex"] == sx), "count"].iloc[0])
                left = -g - n if sx == "Male" else g
                ax.barh(i, n, left=left, height=0.84, color=col, edgecolor="white")
                if n > 0 or show_zero:
                    inside = n >= max_n * 0.10
                    if inside:
                        xl, ha, tc = (left + n / 2), "center", "white"
                    elif sx == "Male":
                        xl, ha, tc = -g - n - max_n * 0.03, "right", col
                    else:
                        xl, ha, tc = g + n + max_n * 0.03, "left", col
                    ax.text(xl, i, str(n), ha=ha, va="center", color=tc,
                            fontsize=8, fontweight="bold")
            ax.text(0, i, lab, ha="center", va="center", color="#404040",
                    fontsize=9, fontweight="bold")
        tot_m = int(tab.loc[tab["sex"] == "Male", "count"].sum())
        tot_f = int(tab.loc[tab["sex"] == "Female", "count"].sum())
        top = len(_AGE_LABELS) + 0.9
        ax.text(-g - max_n / 2, top, f"Hombres (n = {tot_m})", ha="center",
                color=col_male, fontweight="bold")
        ax.text(g + max_n / 2, top, f"Mujeres (n = {tot_f})", ha="center",
                color=col_female, fontweight="bold")
        ax.text(0, top, "Edad", ha="center", color="#404040", fontweight="bold")
        brks = [b for b in np.linspace(0, max_n, 5).round().astype(int)]
        ax.set_xticks([-g - b for b in brks] + [g + b for b in brks])
        ax.set_xticklabels(brks + brks)
        ax.set_xlabel("Número de casos")
        ax.set_yticks([])
        ax.set_ylim(0.4, len(_AGE_LABELS) + 1.5)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        if title:
            fig.suptitle(title, fontweight="bold")
        if subtitle:
            ax.set_title(subtitle, color="gray")
        fig.text(0.99, 0.01, caption, ha="right", fontsize=7, color="gray")
        fig.tight_layout()
    return {"data": tab, "figure": fig}


# --------------------------------------------------------------------------
# plt_serotype
# --------------------------------------------------------------------------
def plt_serotype(data: pd.DataFrame, year: int, cve_edo=None,
                 inset: tuple = (0.58, 0.16, 0.38, 0.38), plot: bool = True) -> dict:
    """Confirmed cases by final diagnosis and by serotype, by state (or by
    municipality of a state).

    Parameters
    ----------
    data : pandas.DataFrame
        SINAVE dengue dataset.
    year : int
        Year (``ANO``).
    cve_edo : str or int, optional
        State; with it the bars are the municipalities of the state.
    inset : tuple
        Position of the serotype plot (left, bottom, width, height, as a
        fraction of the figure); replaces ``x_serotype``, ``y_serotype`` and
        ``scale_serotype`` of R.

    Returns
    -------
    dict
        ``diagnosis`` (area, DES_DIAG_FINAL as DNG/DCSA/DG, n),
        ``serotypes`` (area, DENGUE_SER_TRIPLEX, n) and ``figure``.
    """
    x = data[(_int(data["ANO"]) == year) & (_int(data["ESTATUS_CASO"]) == 2)]
    if cve_edo is not None:
        x = x[_int(x["CVE_EDO_REP"]).isin([int(v) for v in _as_list(cve_edo)])]
        area = "DES_MPO_REP"
    else:
        area = "DES_EDO_REP"
    dg = (x[x["DES_DIAG_FINAL"].isin(_DIAG)].groupby([area, "DES_DIAG_FINAL"]).size()
          .rename("n").reset_index())
    dg = dg[dg[area] != ""]
    dg["area"] = _title(dg[area])
    dg["DES_DIAG_FINAL"] = dg["DES_DIAG_FINAL"].map(
        {"DENGUE CON SIGNOS DE ALARMA": "DCSA", "DENGUE GRAVE": "DG", "DENGUE NO GRAVE": "DNG"})
    sr = x[_int(x["DENGUE_SER_TRIPLEX"]).isin([1, 2, 3, 4])]
    sr = (sr.assign(DENGUE_SER_TRIPLEX=_int(sr["DENGUE_SER_TRIPLEX"]))
            .groupby([area, "DENGUE_SER_TRIPLEX"]).size().rename("n").reset_index())
    sr = sr[sr[area] != ""]
    sr["area"] = _title(sr[area])
    dg = dg[["area", "DES_DIAG_FINAL", "n"]].reset_index(drop=True)
    sr = sr[["area", "DENGUE_SER_TRIPLEX", "n"]].reset_index(drop=True)

    fig = None
    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_axes([0.2, 0.08, 0.75, 0.86])
        _stacked(ax, dg, "DES_DIAG_FINAL", ["DNG", "DCSA", "DG"],
                 {"DNG": "#FFEDA0", "DCSA": "#FEB24C", "DG": "#F03B20"})
        ax.set_xlabel("Casos de Dengue")
        ax.set_ylabel("Municipios" if cve_edo is not None else "Estados")
        if not sr.empty:
            ax2 = fig.add_axes(list(inset))
            _stacked(ax2, sr, "DENGUE_SER_TRIPLEX", [1, 2, 3, 4], _SEROTYPE_COLORS,
                     labels={k: f"DENV{k}" for k in range(1, 5)}, small=True)
            ax2.set_title("Serotipos", fontsize=10, fontstyle="italic", fontweight="bold")
    return {"diagnosis": dg, "serotypes": sr, "figure": fig}


def _stacked(ax, d, col, cats, colors, labels=None, small=False):
    w = d.pivot_table(index="area", columns=col, values="n", aggfunc="sum", fill_value=0)
    w = w.reindex(columns=[c for c in cats if c in w.columns])
    w = w.loc[w.sum(axis=1).sort_values(kind="stable").index]
    left = np.zeros(len(w))
    for c in w.columns:
        v = w[c].to_numpy()
        ax.barh(w.index, v, left=left, color=colors[c], edgecolor="white",
                label=(labels or {}).get(c, c))
        for i, (l, n) in enumerate(zip(left, v)):
            if n > 0:
                ax.text(l + n / 2, i, str(int(n)), ha="center", va="center",
                        fontsize=5 if small else 6)
        left += v
    ax.tick_params(axis="y", labelsize=5 if small else 7)
    if small:
        ax.legend(frameon=False, fontsize=7, loc="lower right")
    else:
        ax.legend(frameon=False, fontsize=10, loc="lower right", ncol=len(w.columns),
                  bbox_to_anchor=(1, 1.0))
    ax.spines[["top", "right"]].set_visible(False)


# --------------------------------------------------------------------------
# plt_treemap
# --------------------------------------------------------------------------
def plt_treemap(snv_dataset: pd.DataFrame, country: bool = True, cve_edo=None,
                year: int | None = None, plot: bool = True) -> dict:
    """Treemap of the cases by final diagnosis, by state (``country=True``) or
    by municipality of a state.

    Returns
    -------
    dict
        ``data`` (area, DES_DIAG_FINAL as DSA/DG/DNG, value) and ``figure``.
    """
    x = snv_dataset[(_int(snv_dataset["ANO"]) == year)
                    & ~snv_dataset["DES_EDO_RES"].isin(_FOREIGN)
                    & snv_dataset["DES_DIAG_FINAL"].isin(_DIAG)]
    if country:
        area = "DES_EDO_REP"
    else:
        x = x[_int(x["CVE_EDO_REP"]).isin([int(v) for v in _as_list(cve_edo)])]
        area = "DES_MPO_REP"
    # as in R, cases without the name of the area are kept (empty name)
    x = x.assign(**{area: x[area].fillna("")})
    tab = x.groupby([area, "DES_DIAG_FINAL"]).size().rename("value").reset_index()
    tab["area"] = _title(tab[area])
    tab["DES_DIAG_FINAL"] = tab["DES_DIAG_FINAL"].map(
        {"DENGUE CON SIGNOS DE ALARMA": "DSA", "DENGUE GRAVE": "DG", "DENGUE NO GRAVE": "DNG"})
    tab = tab[["area", "DES_DIAG_FINAL", "value"]].reset_index(drop=True)

    fig = None
    if plot and not tab.empty:
        try:
            import squarify
        except ImportError:
            warnings.warn("install squarify to draw the treemap", stacklevel=2)
        else:
            import matplotlib.pyplot as plt
            from matplotlib import colormaps
            from matplotlib.patches import Rectangle
            W, H = 100.0, 60.0
            tot = tab.groupby("area")["value"].sum().sort_values(ascending=False)
            areas = sorted(tot.index)
            cmap = colormaps["viridis"].resampled(max(len(areas), 2))
            fig, ax = plt.subplots(figsize=(12, 7))
            outer = squarify.squarify(squarify.normalize_sizes(tot.to_numpy(), W, H),
                                      0, 0, W, H)
            for a, r in zip(tot.index, outer):
                sub = tab[tab["area"] == a].sort_values("value", ascending=False)
                inner = squarify.squarify(
                    squarify.normalize_sizes(sub["value"].to_numpy(), r["dx"], r["dy"]),
                    r["x"], r["y"], r["dx"], r["dy"])
                for (_, row), q in zip(sub.iterrows(), inner):
                    ax.add_patch(Rectangle((q["x"], q["y"]), q["dx"], q["dy"],
                                           facecolor=cmap(areas.index(a)),
                                           edgecolor="white", linewidth=0.4))
                    lab = (f"{a}\n" if country else "") + f"{row['DES_DIAG_FINAL']}\n{row['value']}"
                    ax.text(q["x"] + q["dx"] / 2, q["y"] + 0.6, lab, ha="center",
                            va="bottom", fontsize=7, fontstyle="italic",
                            color="white" if country else "black", alpha=0.8)
                ax.add_patch(Rectangle((r["x"], r["y"]), r["dx"], r["dy"], fill=False,
                                       edgecolor="white", linewidth=2.5))
                if not country:
                    ax.text(r["x"] + r["dx"] / 2, r["y"] + r["dy"] / 2, a, ha="center",
                            va="center", color="white", fontsize=12, fontweight="bold")
            ax.set_xlim(0, W)
            ax.set_ylim(0, H)
            ax.axis("off")
            fig.tight_layout()
    return {"data": tab, "figure": fig}


# --------------------------------------------------------------------------
# mp_bumpmap
# --------------------------------------------------------------------------
def _rescale(v, to):
    v = np.asarray(v, dtype=float)
    if v.max() == v.min():
        return np.full(len(v), np.mean(to))
    return to[0] + (v - v.min()) / (v.max() - v.min()) * (to[1] - to[0])


def mp_bumpmap(dataset: pd.DataFrame, year: int, state: bool = True, cve_edo=None,
               pal_vir: str = "viridis", top_n: int | None = None,
               label_spacing: float | None = None, axis_start: float | None = None,
               smooth: float = 10, plot: bool = True) -> dict:
    """Bump map: map of the cases with sigmoid curves to a ranking.

    ``state=True`` ranks the states of Mexico; ``state=False`` the
    municipalities of ``cve_edo`` (optionally the ``top_n`` with most cases).

    Returns
    -------
    dict
        ``data`` (NOMGEO, value, fine_cap = position in the ranking) and
        ``figure`` (matplotlib).
    """
    import rgeomex

    if not state and cve_edo is None:
        raise ValueError('with state=False give cve_edo (e.g. cve_edo="14")')
    keys = None if cve_edo is None else [str(int(v)).zfill(2) for v in _as_list(cve_edo)]
    if axis_start is None:
        axis_start = 5 if state else 0.5
    x = dataset[(_int(dataset["ANO"]) == year) & ~dataset["DES_EDO_RES"].isin(_FOREIGN)
                & dataset["DES_DIAG_FINAL"].isin(_DIAG)]
    if state:
        df = x.groupby("CVE_EDO_REP").size().rename("value").reset_index()
        base = rgeomex.load("AGEE_inegi19_mx")
        sf = base.merge(df, how="left", left_on="CVE_ENT", right_on="CVE_EDO_REP")
        outline = None
    else:
        x = x[x["CVE_EDO_REP"].isin(keys)]
        df = x.groupby(["CVE_EDO_REP", "CVE_MPO_REP"]).size().rename("value").reset_index()
        base = rgeomex.load("AGEM_inegi19_mx", keys)
        sf = base.merge(df, how="left", left_on=["CVE_ENT", "CVE_MUN"],
                        right_on=["CVE_EDO_REP", "CVE_MPO_REP"])
        outline = rgeomex.load("AGEE_inegi19_mx", keys)
    sf = sf[sf["value"].notna()].copy()
    sf["NOMGEO"] = sf["NOMGEO"].str.strip()
    sf["geometry"] = sf.geometry.make_valid()
    if not state and top_n is not None:
        sf = sf.sort_values("value", ascending=False, kind="stable").head(top_n)
    if sf.empty:
        raise ValueError("no dengue cases for the filters (year / cve_edo)")
    sf = sf.reset_index(drop=True)

    pts = sf.geometry.to_crs(6372).representative_point().to_crs(sf.crs)
    n = len(sf)
    if state:
        y_range, xend = (10, 35), -85
    else:
        bb = outline.total_bounds
        if label_spacing is None:
            y_range = (bb[1], bb[3])
        else:
            mid, half = (bb[1] + bb[3]) / 2, (n - 1) * label_spacing / 2
            y_range = (mid - half, mid + half)
        xend = bb[2]
    rank_first = sf["value"].rank(method="first")
    sf["fine_cap"] = _rescale(rank_first, y_range) if n > 1 else np.mean(y_range)
    tab = sf[["NOMGEO", "value", "fine_cap"]].copy()
    tab["value"] = tab["value"].astype(int)

    fig = None
    if plot:
        import matplotlib.pyplot as plt
        from matplotlib import colormaps
        cmap = colormaps[pal_vir]
        norm = plt.Normalize(*y_range)
        fig, ax = plt.subplots(figsize=(12, 7))
        base.plot(ax=ax, color="#D9D9D9", edgecolor="white", linewidth=0.2)
        sf.plot(ax=ax, color=[cmap(norm(v)) for v in sf["fine_cap"]], edgecolor="white",
                linewidth=0.2)
        if outline is not None:
            outline.boundary.plot(ax=ax, color="black", linewidth=0.3)
        t = np.linspace(-smooth, smooth, 100)
        s = np.exp(t) / (np.exp(t) + 1)
        x1 = xend + axis_start
        for (px, py), y1, v, nm in zip(zip(pts.x, pts.y), sf["fine_cap"], sf["value"],
                                       sf["NOMGEO"]):
            col = cmap(norm(y1))
            ax.plot((t + smooth) / (2 * smooth) * (x1 - px) + px, s * (y1 - py) + py,
                    color=col, lw=1)
            ax.plot(px, py, "o", color=col, ms=2)
            ax.text(x1 + 0.2, y1, f"{nm}  {int(v)}", va="center", fontsize=6, color=col)
        ax.axis("off")
        fig.tight_layout()
    return {"data": tab, "figure": fig}


# --------------------------------------------------------------------------
# mp_cases_loc
# --------------------------------------------------------------------------
def mp_cases_loc(data: gpd.GeoDataFrame, var: str = "probable", top_n: int = 10,
                 title: str = "Casos Probables", r_min: float = 4, r_max: float = 20,
                 cve_edo=None, mun_boundaries: bool = True, edo_boundaries: bool = True):
    """Interactive map (folium) of the cases by locality, with a panel of the
    ``top_n`` localities.

    Parameters
    ----------
    data : geopandas.GeoDataFrame
        Output of :func:`cases_by_loc`.
    var : str
        Column to map (``"probable"``, ``"confirmado"``, ``"DENV2"``...).

    Returns
    -------
    dict
        ``top`` (rank, locality, cases) and ``map`` (folium).
    """
    import folium
    import rgeomex
    from matplotlib import colormaps, colors as mcolors

    d = data.to_crs(4326)
    edo = rgeomex.load("AGEE_inegi19_mx")
    mun = rgeomex.load("AGEM_inegi19_mx") if mun_boundaries else None
    if cve_edo is not None:
        keys = [str(int(v)).zfill(2) for v in _as_list(cve_edo)]
        edo = edo[edo["CVE_ENT"].isin(keys)]
        if mun is not None:
            mun = mun[mun["CVE_ENT"].isin(keys)]
        d = d[d.intersects(edo.union_all())]
    d = d.assign(n=d[var])
    d = d[d["n"] >= 1].sort_values("n", ascending=False, kind="stable").reset_index(drop=True)
    if d.empty:
        raise ValueError(f"there are no localities with {var} >= 1")
    nom = d["nom_loc"].fillna(d["cvegeo"]) if "nom_loc" in d else d["cvegeo"]

    cmap = colormaps["viridis"]
    lo, hi = d["n"].min(), d["n"].max()
    col = lambda v: mcolors.to_hex(cmap(0.5 if hi == lo else (v - lo) / (hi - lo)))  # noqa: E731
    rad = (np.full(len(d), (r_min + r_max) / 2) if hi == lo
           else r_min + (d["n"] - lo) / (hi - lo) * (r_max - r_min))

    bb = edo.total_bounds if cve_edo is not None else (-118.4, 14.5, -86.7, 32.7)
    m = folium.Map(tiles=None, control_scale=True)
    for i, (name, url) in enumerate([
            ("Positron", "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/"
                         "World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}"),
            ("Satelite", "https://server.arcgisonline.com/ArcGIS/rest/services/"
                         "World_Imagery/MapServer/tile/{z}/{y}/{x}"),
            ("OpenStreetMap", "https://tile.openstreetmap.org/{z}/{x}/{y}.png")]):
        folium.TileLayer(url, name=name, attr="Esri / OpenStreetMap", overlay=False,
                         show=(i == 0)).add_to(m)
    if mun is not None:
        folium.GeoJson(mun[["NOMGEO", "geometry"]], name="Límite municipal",
                       style_function=lambda f: {"color": "gray", "weight": 0.5,
                                                 "fillOpacity": 0},
                       tooltip=folium.GeoJsonTooltip(["NOMGEO"], labels=False)).add_to(m)
    if edo_boundaries:
        folium.GeoJson(edo[["geometry"]], name="Límite estatal",
                       style_function=lambda f: {"color": "black", "weight": 1.5,
                                                 "fill": False}).add_to(m)
    fg = folium.FeatureGroup(name=title)
    for (_, r), nm, rr in zip(d.iterrows(), nom, rad):
        pop = (f"<b>Municipio:</b> {r.get('nom_mun', '')}<br><b>Localidad:</b> {nm}<br>"
               f"<b>Probables:</b> {r.get('probable', '')}<br>"
               f"<b>Confirmados:</b> {r.get('confirmado', '')}")
        if "DENV1" in d:
            pop += "<br>" + " | ".join(f"<b>DENV{k}:</b> {r[f'DENV{k}']}" for k in range(1, 5))
        folium.CircleMarker([r.geometry.y, r.geometry.x], radius=float(rr), color="white",
                            weight=0.5, fill=True, fill_color=col(r["n"]), fill_opacity=0.8,
                            popup=pop, tooltip=f"{nm}: {r['n']} {title.lower()}").add_to(fg)
    fg.add_to(m)
    m.fit_bounds([[bb[1], bb[0]], [bb[3], bb[2]]])
    top = pd.DataFrame({"rank": np.arange(1, len(d) + 1), "loc": nom, "n": d["n"]}).head(top_n)
    rows = "".join(f"<tr><td>{r.rank}</td><td><b>{r.loc}</b></td>"
                   f"<td style='background:{col(r.n)};color:white;padding:0 5px'>{r.n}</td></tr>"
                   for r in top.itertuples())
    html = (f"<div style='position:fixed;top:10px;right:60px;z-index:9999;background:white;"
            f"padding:6px;border-radius:4px;font:11px Helvetica'><b>Top {top_n} · {title}</b>"
            f"<table>{rows}</table></div>")
    m.get_root().html.add_child(folium.Element(html))
    folium.LayerControl(collapsed=True).add_to(m)
    return {"top": top.reset_index(drop=True), "map": m}


# --------------------------------------------------------------------------
# mp_serotype
# --------------------------------------------------------------------------
def mp_serotype(path_sinave, cve_edo, palette: str = "viridis", plot: bool = True) -> dict:
    """Serotypes of the confirmed cases by municipality of a state.

    Parameters
    ----------
    path_sinave : str or Path
        The dengue txt file (``"DENGUE2_.txt"``).
    cve_edo : str or int
        State key.

    Returns
    -------
    dict
        ``serotypes`` (municipalities with D1..D4, number of serotypes and
        combination), ``cases`` (municipalities with confirmed cases ``n``)
        and ``map`` (folium).
    """
    import rgeomex

    from ._read import read_sinave

    keys = [str(int(v)).zfill(2) for v in _as_list(cve_edo)]
    x = read_sinave(path_sinave)
    x = x[_int(x["ESTATUS_CASO"]) == 2]
    s = x[_int(x["DENGUE_SER_TRIPLEX"]).isin([1, 2, 3, 4])]
    s = (s.assign(D=_int(s["DENGUE_SER_TRIPLEX"]))
           .groupby(["CVE_EDO_REP", "CVE_MPO_REP", "D"]).size().unstack("D", fill_value=0))
    for k in (1, 2, 3, 4):
        if k not in s:
            s[k] = 0
    s = s[[1, 2, 3, 4]].rename(columns=lambda k: f"D{k}").reset_index()
    s.columns.name = None
    s["n_serotipo"] = (s[["D1", "D2", "D3", "D4"]] > 0).sum(axis=1)
    s["serotype_combination"] = s.apply(
        lambda r: "".join(f"D{k}" for k in range(1, 5) if r[f"D{k}"] > 0), axis=1)
    s = s[s["CVE_EDO_REP"].isin(keys)]
    c = (x[x["CVE_EDO_REP"].isin(keys)].groupby(["CVE_EDO_REP", "CVE_MPO_REP"]).size()
           .rename("n").reset_index())
    y = rgeomex.load("AGEM_inegi19_mx", keys)
    ser = y.merge(s, left_on=["CVE_ENT", "CVE_MUN"], right_on=["CVE_EDO_REP", "CVE_MPO_REP"])
    cas = y.merge(c, left_on=["CVE_ENT", "CVE_MUN"], right_on=["CVE_EDO_REP", "CVE_MPO_REP"])
    ser = ser.drop(columns=["CVE_EDO_REP", "CVE_MPO_REP"])
    cas = cas.drop(columns=["CVE_EDO_REP", "CVE_MPO_REP"])

    m = None
    if plot:
        import folium
        m = y.explore(color="#E6E6E6", name="Municipios", tiles="OpenStreetMap",
                      style_kwds={"fillOpacity": 0.5})
        if not cas.empty:
            cas.explore(m=m, column="n", cmap=palette, name="Número de Casos")
        for k in sorted(ser["n_serotipo"].unique()):
            part = ser[ser["n_serotipo"] == k]
            part.explore(m=m, column="serotype_combination",
                         cmap=palette if k < 4 else None,
                         color="#E01E5A" if k == 4 else None,
                         name=f"{k} Serotipo{'s' if k > 1 else ''}")
        folium.LayerControl().add_to(m)
    return {"serotypes": ser, "cases": cas, "map": m}
