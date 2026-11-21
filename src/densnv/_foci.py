"""Active dengue transmission foci and diagnosis of the epidemiological scenario."""
from __future__ import annotations

import contextlib
import io
import math
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

from ._rt import rt_dengue

CRS_M = 6372          # metres (Mexico ITRF2008 / LCC), as in R
QUAD_SEGS = 30        # segments per quarter circle of the buffers (as sf/GEOS)
MIN_CASOS_DIAG = 50   # minimum cases of the period for a diagnosis
UMBRAL_FOCOS = 0.5    # recent cases in foci: focalised vs dispersed

TXT_BARRIDO = ("Barrido: nebulización con máquina pesada al 100% de la localidad "
               "(4 aplicaciones, una por semana, cada una en menos de 7 días), control "
               "larvario al 100% en máximo 6 semanas (cobertura > 85%, ICP < 3%) y "
               "termonebulización en manzanas calientes (casos de las últimas 2 semanas).")
TXT_MEGA = ("Operativo integral (megaoperativo): nebulización con máquina pesada al 100% "
            "de la localidad (4 aplicaciones semanales), control larvario en el 70% de la "
            "localidad (áreas de mayor transmisión) en máximo 10 días, termonebulización "
            "en manzanas calientes y comando estatal con evaluación semanal.")
TXT_FOCAL = ("Estrategia focalizada en las zonas de control: nebulización con máquina "
             "pesada (3 a 4 aplicaciones, una por semana), control larvario (cobertura > "
             "85%, ICP < 3%) y termonebulización (casos de las últimas 2 semanas); "
             "priorizar los focos significativos (rojos).")


# --------------------------------------------------------------------------
# DBSCAN on precomputed neighbourhoods (as dbscan::dbscan() on a frNN object)
# --------------------------------------------------------------------------
def frnn(xy: np.ndarray, eps: float) -> list[np.ndarray]:
    """Fixed-radius nearest neighbours: for each point, the other points within
    ``eps`` (distance <= eps), as ``dbscan::frNN()``."""
    tree = cKDTree(xy)
    out = []
    for i, ids in enumerate(tree.query_ball_point(xy, r=eps)):
        ids = np.asarray(sorted(ids), dtype=int)
        out.append(ids[ids != i])
    return out


def dbscan_nn(nn: list[np.ndarray], min_pts: int) -> np.ndarray:
    """DBSCAN with border points on precomputed neighbour lists.

    A point is core if it has at least ``min_pts`` points in its neighbourhood
    (itself included). Clusters are numbered 1, 2, ... in the order they are
    found (points in index order); border points join the first cluster that
    reaches them; 0 = noise. Same result as ``dbscan::dbscan(frNN, minPts,
    borderPoints = TRUE)``.
    """
    n = len(nn)
    core = np.array([len(v) + 1 >= min_pts for v in nn])
    label = np.zeros(n, dtype=int)
    k = 0
    for i in range(n):
        if label[i] or not core[i]:
            continue
        k += 1
        label[i] = k
        stack = list(nn[i])
        while stack:
            j = stack.pop()
            if label[j]:
                continue
            label[j] = k
            if core[j]:
                stack.extend(v for v in nn[j] if not label[v])
    return label


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _signif(x, digits=2):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return np.nan
    if x == 0:
        return 0.0
    return round(x, digits - 1 - int(math.floor(math.log10(abs(x)))))


def _r(x, digits=0):
    """round() that keeps NaN."""
    return np.nan if x is None or pd.isna(x) else round(float(x), digits)


def _holm(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    o = np.argsort(p, kind="stable")
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o])
    out = np.empty(m)
    out[o] = np.minimum(adj, 1)
    return out


def _buffer(gs: gpd.GeoSeries, d: float) -> gpd.GeoSeries:
    return gs.buffer(d, quad_segs=QUAD_SEGS)


# --------------------------------------------------------------------------
# main function
# --------------------------------------------------------------------------
def mp_active_foci(geocoded_dataset: pd.DataFrame, cve_edo, locality: str,
                   status_caso=(1, 2), spatial_scale: float = 400,
                   temporal_scale: int = 21, min_cases: int = 3,
                   buffer_m: float = 200, ref_date=None, method: str = "window",
                   n_sim: int = 999, perm_days: int = 365, alpha: float = 0.05,
                   seed: int | None = None, place: str = "REP", mean_si: float = 21,
                   std_si: float = 6.1, capacity_cases: float | None = None,
                   coverage_min: float | None = None, foci_min: float = 0.9,
                   print_map: bool = False, save: bool = False, dir: str = "./") -> dict:
    """Active dengue transmission foci and diagnosis of the epidemiological scenario.

    Python version of ``mp_active_foci()`` of densnv (R). It identifies the
    active transmission foci of a locality with DBSCAN (window) or ST-DBSCAN,
    delimits the transmission and control zones, selects the blocks to treat,
    tests the space-time interaction with the Knox test (global and by focus),
    estimates Rt with :func:`rt_dengue` and diagnoses the epidemiological
    scenario with the operational strategy.

    Parameters
    ----------
    geocoded_dataset : pandas.DataFrame
        Geocoded cases with ``long``, ``lat`` (EPSG:4326), ``ESTATUS_CASO``
        and ``FEC_INI_SIGNOS_SINT``; for Rt also the ``CVE_*_REP`` (or
        ``RES``) keys.
    cve_edo : str or int
        INEGI state key of the locality.
    locality : str
        Name of the locality (as in ``rgeomex.extract_locality()``).
    status_caso : int or list of int
        1 probable, 2 confirmed.
    spatial_scale : float
        Spatial window in metres (flight range of *Aedes aegypti*, 400).
    temporal_scale : int
        Temporal window in days (intrinsic + extrinsic incubation, 21).
    min_cases : int
        Minimum cases in the space-time window to define a focus (3).
    buffer_m : float
        Buffer of the control zone around the transmission zone (200 m).
    ref_date : str or date, optional
        Reference date; by default the latest onset in the locality.
    method : {"window", "st-dbscan"}
        Fixed window or space-time DBSCAN (new and persistent foci).
    n_sim : int
        Permutations of the Knox test (0 skips it).
    perm_days : int
        Days, up to ``ref_date``, of the cases used by the Knox test and by
        ST-DBSCAN (>= ``temporal_scale``).
    alpha : float
        Significance level (Holm-adjusted p-value by focus; Poisson test).
    seed : int, optional
        Seed of the permutations.
    place : {"REP", "RES"}
        Place of the cases used for Rt.
    mean_si, std_si : float
        Serial interval for Rt (21, 6.1 days).
    capacity_cases : float, optional
        Operational capacity: cases per week the brigades can work.
    coverage_min : float, optional
        Coverage of the locality by the recent cases to consider them
        dispersed; by default 0.40 with 21 days, 0.50 with 28 and linear in
        between.
    foci_min : float
        Proportion of the recent cases in foci from which the cases are also
        considered dispersed (0.9).
    print_map : bool
        Return the interactive map (folium) in ``map``.
    save : bool
        Save the selected blocks as a GeoPackage in ``dir``.
    dir : str
        Folder where the blocks are saved.

    Returns
    -------
    dict
        ``foci``, ``transmission``, ``blocks``, ``cases`` (GeoDataFrames in
        EPSG:4326, ``None`` without foci), ``knox`` (global Knox test),
        ``diagnostico`` (``tabla``, ``escenario``, ``accion``, ``detalle``,
        ``rt``, ``rt_figure``, ``casos_ventana``, ``casos_mes``,
        ``pct_pares_tiempo``, ``cobertura_casos``, ``pct_casos_en_focos``,
        ``pct_area_transmision``) and ``map`` (folium, or ``None``).
    """
    import rgeomex

    if method not in ("window", "st-dbscan"):
        raise ValueError("method must be 'window' or 'st-dbscan'")
    if place not in ("REP", "RES"):
        raise ValueError("place must be 'REP' or 'RES'")
    missing = [c for c in ("long", "lat", "ESTATUS_CASO", "FEC_INI_SIGNOS_SINT")
               if c not in geocoded_dataset.columns]
    if missing:
        raise ValueError(f"geocoded_dataset has no column(s): {', '.join(missing)}")
    if perm_days < temporal_scale:
        raise ValueError("perm_days must be greater than or equal to temporal_scale")
    cve_edo = str(int(cve_edo)).zfill(2)
    status = list(status_caso) if isinstance(status_caso, (list, tuple, set, range)) \
        else [status_caso]

    # Step 1. locality
    loc = rgeomex.extract_locality(cve_edo, locality).to_crs(4326)
    loc["geometry"] = loc.geometry.make_valid()
    loc_union = loc.union_all()

    # Step 2. cases of the locality
    d = geocoded_dataset.copy()
    d["onset"] = pd.to_datetime(d["FEC_INI_SIGNOS_SINT"], errors="coerce").dt.normalize()
    d = d[d["long"].notna() & d["lat"].notna()
          & pd.to_numeric(d["ESTATUS_CASO"], errors="coerce").isin(status)
          & d["onset"].notna()]
    z = gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d["long"], d["lat"]), crs=4326)
    z = z[z.intersects(loc_union)].reset_index(drop=True)
    if z.empty:
        raise ValueError(f"there are no geocoded cases with date of onset in {locality}")

    # Step 3. spatial and temporal window
    ref_date = z["onset"].max() if ref_date is None else pd.Timestamp(ref_date).normalize()
    win_ini = ref_date - pd.Timedelta(days=temporal_scale - 1)
    per_ini = ref_date - pd.Timedelta(days=perm_days - 1)
    z_per = z[(z["onset"] >= per_ini) & (z["onset"] <= ref_date)].to_crs(CRS_M) \
        .reset_index(drop=True)

    # Step 4. foci with DBSCAN
    if method == "window":
        z = z[(z["onset"] >= win_ini) & (z["onset"] <= ref_date)].to_crs(CRS_M) \
            .reset_index(drop=True)
        if z.empty:
            raise ValueError(f"there are no cases between {win_ini.date()} and {ref_date.date()}")
        xy = np.column_stack([z.geometry.x, z.geometry.y])
        nn = frnn(xy, spatial_scale)
    else:
        z = z_per.copy()
        if z.empty:
            raise ValueError(f"there are no cases between {per_ini.date()} and {ref_date.date()}")
        xy = np.column_stack([z.geometry.x, z.geometry.y])
        nn = frnn(xy, spatial_scale)
        t = z["onset"].to_numpy().astype("datetime64[D]").astype(int)
        nn = [ids[np.abs(t[ids] - t[i]) <= temporal_scale] for i, ids in enumerate(nn)]

    lab = dbscan_nn(nn, min_cases)
    z["core"] = [len(v) >= min_cases - 1 for v in nn]
    z["foco"] = pd.array(np.where(lab == 0, pd.NA, lab), dtype="Int64")

    if method == "window":
        z["reciente"] = True
    else:
        z["reciente"] = z["onset"] >= win_ini
        activos = sorted(z.loc[z["foco"].notna() & z["reciente"], "foco"].unique())
        renum = {f: i + 1 for i, f in enumerate(activos)}
        z["foco"] = pd.array([renum.get(f, pd.NA) if pd.notna(f) else pd.NA
                              for f in z["foco"]], dtype="Int64")
        z = z[z["reciente"] | z["foco"].notna()].reset_index(drop=True)
        if not z["reciente"].any():
            raise ValueError(f"there are no cases between {win_ini.date()} and {ref_date.date()}")

    # Step 4.1. Knox test
    def prueba_knox(zonas: gpd.GeoDataFrame | None = None):
        rng = np.random.default_rng(seed)
        tt = z_per["onset"].to_numpy().astype("datetime64[D]").astype(int)
        n_p = len(tt)
        nn_p = frnn(np.column_stack([z_per.geometry.x, z_per.geometry.y]), spatial_scale)
        i = np.repeat(np.arange(n_p), [len(v) for v in nn_p])
        j = np.concatenate(nn_p) if n_p else np.array([], dtype=int)
        keep = i < j
        i, j = i[keep], j[keep]
        en_zona = []
        if zonas is not None:
            for g in zonas.geometry:
                inc = z_per.intersects(g).to_numpy()
                en_zona.append(inc[i] & inc[j])

        def cuenta(tv):
            cerca = np.abs(tv[i] - tv[j]) <= temporal_scale
            return np.array([cerca.sum()] + [(cerca & e).sum() for e in en_zona], dtype=float)

        obs = cuenta(tt)
        print(f"Knox test: {n_sim} permutations of {n_p} cases...")
        sim = np.column_stack([cuenta(tt[rng.permutation(n_p)]) for _ in range(n_sim)]) \
            if n_sim > 0 else np.zeros((len(obs), 0))
        p = (1 + (sim >= obs[:, None]).sum(axis=1)) / (n_sim + 1)
        esp = sim.mean(axis=1)
        rr = np.where(esp > 0, obs / np.where(esp > 0, esp, 1), np.nan)
        glob = {"n_sim": n_sim, "periodo": (per_ini, ref_date), "n_casos": n_p,
                "pares_espacio": int(len(i)), "knox": int(obs[0]), "esperados": esp[0],
                "rr": rr[0], "p_value": p[0], "sim": sim[0]}
        focos = None
        if zonas is not None:
            focos = pd.DataFrame({"foco": zonas["foco"].to_numpy(),
                                  "knox_pares": obs[1:].astype(int),
                                  "knox_esperados": np.round(esp[1:], 2),
                                  "knox_rr": np.round(rr[1:], 2),
                                  "p_value": p[1:], "p_holm": _holm(p[1:])})
            focos["significativo"] = focos["p_holm"] < alpha
        return glob, focos

    # Step 5. transmission and control zones
    zf = z[z["foco"].notna() & z["reciente"]]
    transmission = foci = blocks_sel = None
    if zf.empty:
        warnings.warn(f"no active foci between {win_ini.date()} and {ref_date.date()} "
                      f"({min_cases} or more cases within {spatial_scale} m and "
                      f"{temporal_scale} days)", stacklevel=2)
        knox = prueba_knox()[0] if n_sim > 0 else None
    else:
        bufs = gpd.GeoDataFrame({"foco": zf["foco"].astype(int).to_numpy()},
                                geometry=_buffer(zf.geometry, spatial_scale).to_numpy(),
                                crs=CRS_M)
        transmission = bufs.dissolve(by="foco").reset_index()[["foco", "geometry"]]
        control = transmission.copy()
        control["geometry"] = _buffer(transmission.geometry, buffer_m)

        # blocks of the state near the control zones
        bl = rgeomex.load("blocks_ine20_mx", cve_edo)
        bl = bl[pd.to_numeric(bl["entidad"], errors="coerce") == int(cve_edo)]
        xmin, ymin, xmax, ymax = control.to_crs(bl.crs).total_bounds
        bl = bl.cx[xmin:xmax, ymin:ymax].copy()
        bl["geometry"] = bl.geometry.make_valid()
        bl = bl.to_crs(CRS_M).reset_index(drop=True)
        bl["id_manzana"] = np.arange(len(bl))
        blocks_sel = gpd.sjoin(bl, control[["foco", "geometry"]], how="inner",
                               predicate="intersects")
        blocks_sel = (blocks_sel.sort_values(["id_manzana", "foco"])
                      .drop_duplicates("id_manzana")
                      .drop(columns=["id_manzana", "index_right"]).reset_index(drop=True))

        g = z[z["foco"].notna()].groupby("foco")
        resumen = pd.DataFrame({
            "casos": g.size(),
            "casos_recientes": g["reciente"].sum(),
            "confirmados": g["ESTATUS_CASO"].apply(lambda s: int((pd.to_numeric(s) == 2).sum())),
            "primer_inicio": g["onset"].min(),
            "ultimo_inicio": g["onset"].max()}).reset_index()
        resumen["foco"] = resumen["foco"].astype(int)
        resumen["duracion_dias"] = (resumen["ultimo_inicio"] - resumen["primer_inicio"]).dt.days + 1
        resumen["tipo"] = np.where(resumen["primer_inicio"] >= win_ini, "Nuevo", "Persistente") \
            if method == "st-dbscan" else None
        n_manz = blocks_sel.groupby("foco").size().rename("manzanas").reset_index()

        if n_sim > 0:
            knox, kf = prueba_knox(transmission)
            resumen = resumen.merge(kf, on="foco", how="left")
        else:
            knox = None
            for c in ("knox_pares", "knox_esperados", "knox_rr", "p_value", "p_holm"):
                resumen[c] = np.nan
            resumen["significativo"] = pd.NA

        foci = control.copy()
        foci["area_control_km2"] = np.round(control.area.to_numpy() / 1e6, 2)
        foci["area_transmision_km2"] = np.round(transmission.area.to_numpy() / 1e6, 2)
        foci = foci.merge(resumen, on="foco", how="left").merge(n_manz, on="foco", how="left")
        foci["manzanas"] = foci["manzanas"].fillna(0).astype(int)
        cols = ["foco", "tipo", "casos", "casos_recientes", "confirmados", "knox_pares",
                "knox_esperados", "knox_rr", "p_value", "p_holm", "significativo",
                "primer_inicio", "ultimo_inicio", "duracion_dias", "area_transmision_km2",
                "area_control_km2", "manzanas"]
        foci = foci[cols + ["geometry"]]

        transmission = transmission.to_crs(4326)
        foci = foci.to_crs(4326)
        blocks_sel = blocks_sel.to_crs(4326)
        if save:
            f = Path(dir) / f"{locality.lower().replace(' ', '_')}_active_foci_blocks_{ref_date.date()}.gpkg"
            blocks_sel.to_file(f, driver="GPKG")
            print(f"Blocks saved in {f}")

    # Step 5.2. diagnosis of the epidemiological scenario
    umbral_dispersion = coverage_min if coverage_min is not None else \
        min(max(0.40 + (temporal_scale - 21) * 0.10 / 7, 0.40), 0.50)
    loc_m = gpd.GeoSeries([loc_union], crs=4326).to_crs(CRS_M).make_valid().union_all()
    area_loc = loc_m.area

    ts = np.sort(z_per["onset"].to_numpy().astype("datetime64[D]").astype(int))
    n_t = len(ts)
    pct_pares_tiempo = (np.sum(np.searchsorted(ts, ts + temporal_scale, side="right")
                               - np.arange(1, n_t + 1)) / math.comb(n_t, 2)) if n_t > 1 else np.nan
    casos_mes = (z_per["onset"].dt.strftime("%Y-%m").value_counts().sort_index()
                 .rename_axis("mes").rename("casos").reset_index())
    casos_mes["pct"] = np.round(100 * casos_mes["casos"] / casos_mes["casos"].sum(), 1)

    # criterion 1. sustained and significant increase
    n_vent = perm_days // temporal_scale
    k = ((ref_date - z_per["onset"]).dt.days // temporal_scale).to_numpy()
    cv = np.bincount(k[k < n_vent], minlength=n_vent)[:n_vent] if n_vent else np.array([])
    base = float(np.mean(cv[2:])) if n_vent >= 4 else np.nan
    if np.isnan(base):
        p_inc = np.nan
    elif base == 0:
        p_inc = float(cv[0] == 0)
    else:
        p_inc = float(stats.poisson.sf(cv[0] - 1, base))
    razon = np.nan if (np.isnan(base) or base == 0) else cv[0] / base

    try:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rt_loc = rt_dengue(geocoded_dataset, level="locality", place=place,
                               cve_edo=cve_edo, cve_mpo=loc["CVE_MUN"].tolist(),
                               cve_loc=loc["CVE_LOC"].tolist(), status_caso=status,
                               ref_date=ref_date, period_days=perm_days,
                               window_days=temporal_scale, mean_si=mean_si,
                               std_si=std_si, plot=True)
    except Exception as e:  # noqa: BLE001
        print(f"Rt not estimated: {e}")
        rt_loc = None
    rv = None if rt_loc is None else rt_loc["rt_ventana"].iloc[0]
    rt_med = np.nan if rv is None else rv["rt"]
    rt_inf = np.nan if rv is None else rv["rt_inf"]
    rt_sup = np.nan if rv is None else rv["rt_sup"]
    fase = "Sin estimar" if rv is None else rv["fase"]

    incremento = bool(not np.isnan(p_inc) and p_inc < alpha and cv[1] > base
                      and (np.isnan(rt_med) or rt_med > 1))

    # criterion 2. operational capacity
    n_manzanas = 0 if blocks_sel is None else len(blocks_sel)
    casos_semana = cv[0] / (temporal_scale / 7)
    rebasa = None if capacity_cases is None else bool(casos_semana > capacity_cases)

    # criterion 3. dispersion of the recent cases
    zr = z[z["reciente"]]
    cobertura = 0.0 if zr.empty else \
        _buffer(zr.geometry, spatial_scale).union_all().intersection(loc_m).area / area_loc
    n_rec = int(z["reciente"].sum())
    pct_casos_focos = 0.0 if n_rec == 0 else len(zf) / n_rec
    disperso = bool(cobertura >= umbral_dispersion or pct_casos_focos >= foci_min)

    n_focos = 0 if foci is None else len(foci)
    pct_area = 0.0 if transmission is None else \
        transmission.to_crs(CRS_M).union_all().intersection(loc_m).area / area_loc

    if knox is None:
        interaccion = "Sin prueba"
    elif np.isnan(knox["rr"]):
        interaccion = "Sin datos"
    elif knox["p_value"] < alpha:
        interaccion = "Sí"
    elif knox["rr"] > 1.05:
        interaccion = "No concluyente"
    else:
        interaccion = "No"

    generalizada = incremento and (disperso or rebasa is True)
    if len(z_per) < MIN_CASOS_DIAG:
        escenario = "Casos insuficientes"
    elif generalizada:
        escenario = "Transmisión generalizada"
    elif n_focos == 0:
        escenario = "Sin focos activos"
    elif pct_casos_focos >= UMBRAL_FOCOS:
        escenario = "Transmisión focalizada"
    else:
        escenario = "Transmisión focalizada con casos dispersos"
    accion = {
        "Casos insuficientes": f"Vigilancia: menos de {MIN_CASOS_DIAG} casos en el periodo "
                               "para el diagnóstico.",
        "Transmisión generalizada": (TXT_BARRIDO + " Si los casos por semana rebasan la "
                                     "capacidad operativa: operativo integral (indicar "
                                     "capacity_cases).") if rebasa is None else
                                    (TXT_MEGA if rebasa else TXT_BARRIDO),
        "Sin focos activos": "Vigilancia entomológica (ovitrampas, entomovirológica) y "
                             "control alrededor de los casos.",
        "Transmisión focalizada": TXT_FOCAL,
        "Transmisión focalizada con casos dispersos":
            TXT_FOCAL + " Además, control alrededor de los casos aislados.",
    }[escenario]
    alerta = bool(incremento and fase == "En ascenso" and not generalizada
                  and escenario != "Casos insuficientes")
    if alerta:
        accion = "Alerta temprana (incremento significativo, Rt en ascenso). " + accion
    if not generalizada and rebasa is True:
        accion += (f" Los casos por semana ({round(casos_semana, 1)}) rebasan la capacidad "
                   f"operativa ({capacity_cases}).")

    fmt = lambda v: f"{v:.2f}".rstrip("0").rstrip(".")  # noqa: E731
    detalle = [
        (f"Incremento: {cv[0]} casos en la ventana actual y {cv[1]} en la previa; media por "
         f"ventana de {temporal_scale} días en el periodo: {_r(base, 1)}"
         + (f" (razón {fmt(razon)}, p = {_signif(p_inc)})" if not np.isnan(razon) else "")
         + "; Rt de la ventana actual: "
         + ("sin estimar" if np.isnan(rt_med) else
            f"{fmt(rt_med)} (IC95% {fmt(rt_inf)} a {fmt(rt_sup)}), {fase.lower()}")
         + ". " + ("Sostenido y significativo." if incremento else "No.")),
        (f"Capacidad operativa: {round(casos_semana, 1)} casos por semana a trabajar"
         + (" (capacidad no indicada)." if rebasa is None else
            f" contra {capacity_cases} de capacidad. "
            + ("La rebasa." if rebasa else "No la rebasa."))),
        (f"Dispersión: los casos recientes (a {spatial_scale} m) cubren {round(100 * cobertura)}% "
         f"de la localidad (umbral {round(100 * umbral_dispersion)}%); casos recientes en focos: "
         f"{round(100 * pct_casos_focos)}% (umbral {round(100 * foci_min)}%). "
         + ("Dispersos en la localidad." if disperso else "No dispersos en toda la localidad.")),
        (f"Interacción espacio-temporal (Knox): {interaccion}"
         + (f" (RR = {fmt(knox['rr'])}, p = {_signif(knox['p_value'])})"
            if knox is not None and not np.isnan(knox["rr"]) else "") + "."),
    ]
    tabla = pd.DataFrame([{
        "localidad": locality, "ref_date": ref_date, "casos_periodo": len(z_per),
        "casos_ventana": int(cv[0]), "razon_incremento": _r(razon, 2),
        "p_incremento": _signif(p_inc), "rt": _r(rt_med, 2), "rt_inf": _r(rt_inf, 2),
        "rt_sup": _r(rt_sup, 2), "fase": fase, "incremento": incremento,
        "alerta_temprana": alerta, "cobertura_casos": _r(100 * cobertura, 1),
        "umbral_cobertura": round(100 * umbral_dispersion), "dispersos": disperso,
        "pct_casos_en_focos": round(100 * pct_casos_focos),
        "pct_area_focos": _r(100 * pct_area, 1), "focos": n_focos,
        "focos_signif": 0 if n_focos == 0 else int(foci["significativo"].fillna(False).sum()),
        "manzanas": n_manzanas, "casos_semana": _r(casos_semana, 1),
        "capacidad": np.nan if capacity_cases is None else capacity_cases,
        "rebasa_capacidad": rebasa, "interaccion": interaccion,
        "knox_rr": np.nan if knox is None else _r(knox["rr"], 2),
        "knox_p": np.nan if knox is None else _signif(knox["p_value"]),
        "escenario": escenario, "accion": accion}])
    diagnostico = {
        "tabla": tabla, "escenario": escenario, "accion": accion, "detalle": detalle,
        "rt": None if rt_loc is None else rt_loc["rt"],
        "rt_figure": None if rt_loc is None else rt_loc["figure"],
        "casos_ventana": cv, "casos_mes": casos_mes,
        "pct_pares_tiempo": pct_pares_tiempo, "cobertura_casos": cobertura,
        "pct_casos_en_focos": pct_casos_focos, "pct_area_transmision": pct_area}

    cases = z.to_crs(4326)
    cases["grupo"] = np.select(
        [cases["foco"].notna() & cases["reciente"], cases["foco"].notna() & ~cases["reciente"]],
        ["Caso en foco", "Caso previo del foco"], default="Caso aislado")

    n_cases_foci = len(zf)
    msg = (f"{n_focos} active foci, {n_cases_foci} of {len(z)} cases in foci, {n_manzanas} "
           f"blocks to control ({win_ini.date()} to {ref_date.date()}, method = {method}).")
    if knox is not None:
        msg += (f" Knox test ({knox['n_sim']} permutations, {knox['n_casos']} cases from "
                f"{per_ini.date()}): global {knox['knox']} pairs, expected "
                f"{knox['esperados']:.1f}, RR = {knox['rr']:.2f}, p = {_signif(knox['p_value'])}.")
        if foci is not None:
            msg += (f" {int(foci['significativo'].sum())} of {len(foci)} foci with Holm "
                    f"p < {alpha}.")
    print(msg)
    print(f"Diagnóstico: {escenario}\n"
          + "\n".join(f"  {i}. {t}" for i, t in enumerate(detalle, 1))
          + f"\nAcción: {accion}")

    m = _map(loc, foci, transmission, blocks_sel, cases, alpha, ref_date, spatial_scale,
             buffer_m, method) if print_map else None
    return {"map": m, "foci": foci, "transmission": transmission, "blocks": blocks_sel,
            "cases": cases, "knox": knox, "diagnostico": diagnostico}


def _map(loc, foci, transmission, blocks, cases, alpha, ref_date, spatial_scale,
         buffer_m, method):
    """Interactive map (folium). Returns None when folium is not installed."""
    try:
        import folium
    except ImportError:
        warnings.warn("install folium to draw the map: pip install folium", stacklevel=3)
        return None

    b = loc.total_bounds
    m = folium.Map(tiles=None, control_scale=True)
    m.fit_bounds([[b[1], b[0]], [b[3], b[2]]])
    for i, (name, url, attr) in enumerate([
        ("Esri Light Gray", "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/"
         "World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}", "Tiles &copy; Esri"),
        ("Esri Dark Gray", "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/"
         "World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}", "Tiles &copy; Esri"),
        ("OpenStreetMap", "https://tile.openstreetmap.org/{z}/{x}/{y}.png",
         "&copy; OpenStreetMap contributors"),
        ("World Imagery", "https://server.arcgisonline.com/ArcGIS/rest/services/"
         "World_Imagery/MapServer/tile/{z}/{y}/{x}", "Tiles &copy; Esri")]):
        folium.TileLayer(url, name=name, attr=attr, overlay=False,
                         show=(i == 0)).add_to(m)

    folium.GeoJson(loc[["geometry"]], name="Límite de la localidad",
                   style_function=lambda f: {"color": "#444444", "weight": 1.2,
                                             "fill": False}).add_to(m)
    if foci is not None:
        col = {"Significativo": "#D7191C", "RR > 1": "#FFC300", "RR <= 1": "#969696",
               "Sin prueba": "#D7191C"}
        tr = transmission.merge(foci[["foco", "significativo", "knox_rr", "knox_pares"]],
                                on="foco")
        tr["clase"] = np.select(
            [tr["significativo"].isna(), tr["significativo"].fillna(False).astype(bool),
             (tr["knox_rr"] > 1) | (tr["knox_rr"].isna() & (tr["knox_pares"] > 0))],
            ["Sin prueba", "Significativo", "RR > 1"], default="RR <= 1")
        ctrl = foci.copy()
        for c in ("primer_inicio", "ultimo_inicio"):
            ctrl[c] = ctrl[c].dt.strftime("%Y-%m-%d")
        folium.GeoJson(ctrl, name="Zona de control",
                       style_function=lambda f: {"color": "#F46D43", "weight": 1.5,
                                                 "dashArray": "4 4", "fillColor": "#FDAE61",
                                                 "fillOpacity": 0.25},
                       popup=folium.GeoJsonPopup(
                           ["foco", "tipo", "casos", "casos_recientes", "confirmados",
                            "knox_pares", "knox_rr", "p_holm", "manzanas",
                            "area_control_km2"])).add_to(m)
        folium.GeoJson(tr[["foco", "clase", "geometry"]], name="Zona de transmisión",
                       style_function=lambda f: {"color": col[f["properties"]["clase"]],
                                                 "weight": 0,
                                                 "fillColor": col[f["properties"]["clase"]],
                                                 "fillOpacity": 0.4}).add_to(m)
        folium.GeoJson(blocks[["foco", "geometry"]], name="Manzanas a controlar",
                       style_function=lambda f: {"color": "#7F2704", "weight": 0.6,
                                                 "fill": False}).add_to(m)
    colores = {"Caso en foco": "#B2182B", "Caso previo del foco": "#F4A582",
               "Caso aislado": "#4D4D4D"}
    fg = folium.FeatureGroup(name="Casos")
    for _, r in cases.iterrows():
        folium.CircleMarker([r.geometry.y, r.geometry.x], radius=4, weight=0.8,
                            color="white", fill=True, fill_color=colores[r["grupo"]],
                            fill_opacity=1,
                            popup=f"<b>{r['grupo']}</b>"
                                  + ("" if pd.isna(r["foco"]) else f" (foco {r['foco']})")
                                  + f"<br>Inicio de síntomas: {r['onset']:%Y-%m-%d}").add_to(fg)
    fg.add_to(m)
    folium.LayerControl(collapsed=True).add_to(m)
    return m
