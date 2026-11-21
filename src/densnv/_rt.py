"""Effective reproductive number (Rt) of dengue (Cori et al. 2013)."""
from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
from scipy import stats

_QUANTILES = {"q025": 0.025, "q05": 0.05, "q25": 0.25, "median": 0.5,
              "q75": 0.75, "q95": 0.95, "q975": 0.975}


def discr_si(k, mu: float, sigma: float) -> np.ndarray:
    """Discretised serial interval: gamma distribution shifted by one day.

    Same as ``EpiEstim::discr_si()``.
    """
    if sigma < 0:
        raise ValueError("sigma must be >= 0")
    if mu <= 1:
        raise ValueError("mu must be > 1")
    k = np.asarray(k, dtype=float)
    a = ((mu - 1) / sigma) ** 2
    b = sigma ** 2 / (mu - 1)

    def cdf(x, shape):
        return stats.gamma.cdf(x, shape, scale=b)

    res = k * cdf(k, a) + (k - 2) * cdf(k - 2, a) - 2 * (k - 1) * cdf(k - 1, a)
    res = res + a * b * (2 * cdf(k - 1, a + 1) - cdf(k - 2, a + 1) - cdf(k, a + 1))
    return np.maximum(res, 0)


def estimate_r(incidence, t_start, t_end, mean_si: float = 21, std_si: float = 6.1,
               mean_prior: float = 5, std_prior: float = 5) -> pd.DataFrame:
    """Rt with a parametric serial interval (``EpiEstim::estimate_R()``,
    method ``"parametric_si"``).

    Parameters
    ----------
    incidence : array-like
        Daily incidence (cases per day, consecutive days). As in EpiEstim, the
        cases of the first day are taken as imported.
    t_start, t_end : array-like of int
        Start and end (1-based, inclusive) of each time window.
    mean_si, std_si : float
        Mean and standard deviation (days) of the serial interval.
    mean_prior, std_prior : float
        Mean and standard deviation of the gamma prior of Rt.

    Returns
    -------
    pandas.DataFrame
        One row per window: ``t_start``, ``t_end``, ``mean``, ``std`` and the
        quantiles ``q025`` ... ``q975`` of the posterior of Rt.
    """
    inc = np.nan_to_num(np.asarray(incidence, dtype=float))
    T = len(inc)
    t_start = np.atleast_1d(np.asarray(t_start, dtype=int))
    t_end = np.atleast_1d(np.asarray(t_end, dtype=int))
    if len(t_start) != len(t_end) or np.any(t_start < 2) or np.any(t_end > T) \
            or np.any(t_start > t_end):
        raise ValueError("t_start and t_end must be windows within 2..len(incidence)")

    local = inc.copy()
    local[0] = 0                                 # first day: imported
    si = np.zeros(T + 1)
    si[:T] = discr_si(np.arange(T), mean_si, std_si)
    final_mean_si = float(np.sum(si * np.arange(T + 1)))

    # overall infectivity: lambda_t = sum_{s=0}^{t-1} w_s I_{t-s}
    lam = np.full(T, np.nan)
    for t in range(1, T):                        # 0-based day t
        lam[t] = np.sum(si[: t + 1] * inc[t::-1])

    a_prior = (mean_prior / std_prior) ** 2
    b_prior = std_prior ** 2 / mean_prior
    rows = []
    for s, e in zip(t_start, t_end):
        if e <= final_mean_si:
            continue
        a = a_prior + local[s - 1: e].sum()
        b = 1 / (1 / b_prior + lam[s - 1: e].sum())
        row = {"t_start": int(s), "t_end": int(e), "mean": a * b,
               "std": math.sqrt(a) * b}
        for name, p in _QUANTILES.items():
            row[name] = stats.gamma.ppf(p, a, scale=b)
        rows.append(row)
    return pd.DataFrame(rows, columns=["t_start", "t_end", "mean", "std",
                                       *_QUANTILES])


def _num(s):
    return pd.to_numeric(pd.Series(s, dtype="string").str.strip(),
                         errors="coerce").astype("Int64")


def rt_dengue(dataset: pd.DataFrame, level: str = "country", place: str = "REP",
              cve_edo=None, cve_mpo=None, cve_loc=None, status_caso=(1, 2),
              ref_date=None, period_days: int = 365, window_days: int = 21,
              mean_si: float = 21, std_si: float = 6.1, plot: bool = True) -> dict:
    """Effective reproductive number (Rt) of dengue.

    Rt for the country, a state, a municipality or a locality, from the daily
    incidence of the SINAVE cases and a parametric serial interval (Cori et
    al. 2013), as ``rt_dengue()`` of densnv in R.

    Parameters
    ----------
    dataset : pandas.DataFrame
        SINAVE dengue dataset (e.g. ``densnv.read_sinave("DENGUE2_.txt")``).
        It needs ``FEC_INI_SIGNOS_SINT``, ``ESTATUS_CASO`` and, depending on
        ``level`` and ``place``, ``CVE_EDO_*``, ``CVE_MPO_*`` and
        ``CVE_LOC_*``.
    level : {"country", "state", "municipality", "locality"}
    place : {"REP", "RES"}
        Place used to select the cases: of report (default; the cases are
        geocoded with it) or of residence. In non-endemic areas the cases are
        reported where they were exposed, so the place of report is the place
        of transmission.
    cve_edo, cve_mpo, cve_loc : str, int or list
        INEGI keys of the state, municipality and locality. ``cve_mpo`` and
        ``cve_loc`` can be lists of the same length (pairs).
    status_caso : int or list of int
        1 probable, 2 confirmed.
    ref_date : str or date, optional
        Reference date; by default the latest date of onset of the cases.
    period_days : int
        Days, up to ``ref_date``, of the daily incidence (default 365).
    window_days : int
        Days, up to ``ref_date``, of the window whose Rt gives the phase of
        the outbreak (default 21).
    mean_si, std_si : float
        Mean and standard deviation (days) of the serial interval (21, 6.1).
    plot : bool
        Draw the Rt curve (matplotlib).

    Returns
    -------
    dict
        ``rt`` (Rt of the weekly sliding windows with at least 10 cases, or
        ``None``), ``rt_ventana`` (Rt of the window, its 95% credible interval
        and the phase: "En ascenso", "Estable", "En descenso" or
        "Sin estimar"), ``incidencia`` (daily incidence), ``area`` (name of
        the area) and ``figure`` (or ``None``).

    Notes
    -----
    The curve uses weekly sliding windows with at least 10 cases that start
    at least ``mean_si + 2 * std_si`` days after the first date of onset of
    the dataset; before it the infectiousness of the previous cases is
    unknown and Rt is overestimated. The last days are usually incomplete
    (delay of report), which lowers Rt at the end.
    """
    if level not in ("country", "state", "municipality", "locality"):
        raise ValueError("level must be 'country', 'state', 'municipality' or 'locality'")
    if place not in ("REP", "RES"):
        raise ValueError("place must be 'REP' or 'RES'")
    col_edo, col_mpo, col_loc = (f"CVE_{k}_{place}" for k in ("EDO", "MPO", "LOC"))
    need = ["FEC_INI_SIGNOS_SINT", "ESTATUS_CASO"] + {
        "country": [], "state": [col_edo], "municipality": [col_edo, col_mpo],
        "locality": [col_edo, col_mpo, col_loc]}[level]
    missing = [c for c in need if c not in dataset.columns]
    if missing:
        raise ValueError(f"the dataset does not have the columns: {', '.join(missing)}")
    as_list = lambda v: list(v) if isinstance(v, (list, tuple, set, range)) else [v]  # noqa: E731
    if level != "country" and cve_edo is None:
        raise ValueError("cve_edo is required")
    if level in ("municipality", "locality") and cve_mpo is None:
        raise ValueError("cve_mpo is required")
    if level == "locality":
        if cve_loc is None:
            raise ValueError("cve_loc is required")
        if len(as_list(cve_loc)) != len(as_list(cve_mpo)):
            raise ValueError("cve_mpo and cve_loc must have the same length")

    # Step 1. cases of the area
    onset_all = pd.to_datetime(dataset["FEC_INI_SIGNOS_SINT"], errors="coerce")
    datos_ini = onset_all.min()
    x = dataset.assign(onset=onset_all.dt.normalize())
    x = x[pd.to_numeric(x["ESTATUS_CASO"], errors="coerce").isin(as_list(status_caso))]
    if level != "country":
        x = x[_num(x[col_edo]).isin([int(v) for v in as_list(cve_edo)]).to_numpy()]
    if level == "municipality":
        x = x[_num(x[col_mpo]).isin([int(v) for v in as_list(cve_mpo)]).to_numpy()]
    if level == "locality":
        pairs = {(int(m), int(l)) for m, l in zip(as_list(cve_mpo), as_list(cve_loc))}
        keys = list(zip(_num(x[col_mpo]), _num(x[col_loc])))
        x = x[[k in pairs for k in keys]]
    x = x[x["onset"].notna()]

    def name(col):
        if col in x and len(x):
            return ", ".join(pd.unique(x[col].dropna().astype(str)))
        return ""
    area = {"country": "México", "state": name(f"DES_EDO_{place}"),
            "municipality": name(f"DES_MPO_{place}"),
            "locality": name(f"DES_LOC_{place}")}[level]
    if not area:
        area = "-".join(str(v) for v in as_list(cve_edo) + as_list(cve_mpo) + as_list(cve_loc)
                        if v is not None)

    # Step 2. daily incidence of the period
    if ref_date is None:
        if x.empty:
            raise ValueError("there are no cases in the selected area")
        ref_date = x["onset"].max()
    ref_date = pd.Timestamp(ref_date).normalize()
    per_ini = ref_date - pd.Timedelta(days=period_days - 1)
    win_ini = ref_date - pd.Timedelta(days=window_days - 1)
    dates = pd.date_range(per_ini, ref_date, freq="D")
    counts = x["onset"].value_counts()
    I = counts.reindex(dates, fill_value=0).to_numpy()
    incidencia = pd.DataFrame({"dates": dates, "I": I})

    # Step 3. Rt
    historia = math.ceil(mean_si + 2 * std_si)
    rt = rt_win = None
    if I.sum() > 0 and period_days > window_days:
        t_end = np.arange(8, period_days + 1)
        t_start = t_end - 6
        cum = np.concatenate([[0], np.cumsum(I)])
        casos_w = cum[t_end] - cum[t_start - 1]
        ok = (casos_w >= 10) & (dates[t_start - 1] >= datos_ini + pd.Timedelta(days=historia))
        if ok.any():
            rt = estimate_r(I, t_start[ok], t_end[ok], mean_si, std_si)
            rt.insert(2, "date_start", dates[rt["t_start"] - 1])
            rt.insert(3, "date_end", dates[rt["t_end"] - 1])
            if rt.empty:
                rt = None
        w = estimate_r(I, [period_days - window_days + 1], [period_days], mean_si, std_si)
        rt_win = None if w.empty else w.iloc[0]
    if win_ini < datos_ini + pd.Timedelta(days=historia):
        rt_win = None
    if rt_win is None:
        rt_med = rt_inf = rt_sup = np.nan
        fase = "Sin estimar"
    else:
        rt_med, rt_inf, rt_sup = rt_win["mean"], rt_win["q025"], rt_win["q975"]
        fase = "En ascenso" if rt_inf > 1 else "En descenso" if rt_sup < 1 else "Estable"
    rt_ventana = pd.DataFrame([{"inicio": win_ini, "fin": ref_date,
                                "casos": int(I[dates >= win_ini].sum()),
                                "rt": rt_med, "rt_inf": rt_inf, "rt_sup": rt_sup,
                                "fase": fase}])
    if I.sum() < 50:
        warnings.warn(f"fewer than 50 cases in the period ({I.sum()}): Rt is unreliable",
                      stacklevel=2)

    fig = _plot_rt(rt, area, win_ini, ref_date, rt_med, rt_inf, rt_sup, fase) \
        if plot and rt is not None else None
    print(f"{area}: {I.sum()} cases ({per_ini.date()} to {ref_date.date()}); "
          f"Rt {win_ini.date()} to {ref_date.date()} = {rt_med:.2f} "
          f"(95% CrI {rt_inf:.2f} to {rt_sup:.2f}), {fase.lower()}.")
    return {"rt": rt, "rt_ventana": rt_ventana, "incidencia": incidencia,
            "area": area, "figure": fig}


def _plot_rt(rt, area, win_ini, ref_date, rt_med, rt_inf, rt_sup, fase):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.fill_between(rt["date_end"], rt["q025"], rt["q975"], color="red", alpha=0.3,
                    lw=0, label="IC95%")
    ax.plot(rt["date_end"], rt["mean"], color="red", lw=1.5, label="Rt")
    ax.axhline(1, ls="--", color="gray")
    ax.axvline(win_ini, ls=":", color="blue")
    ax.set_ylabel("Número reproductivo efectivo (Rt)")
    ax.set_title(f"{area}. Rt de la ventana ({win_ini:%d/%m} a {ref_date:%d/%m/%Y}): "
                 f"{rt_med:.2f} (IC95% {rt_inf:.2f} a {rt_sup:.2f}), {fase.lower()}",
                 fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.autofmt_xdate()
    fig.tight_layout()
    return fig
