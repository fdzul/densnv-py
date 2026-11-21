"""Read the SINAVE vector-borne disease datasets."""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

#: columns kept by ``read(complete=False)`` (same as densnv in R)
ARBO_COLUMNS = [
    "VEC_ID", "FOL_ID", "IDE_EDA_ANO", "IDE_SEX", "DES_CAL", "IDE_CAL",
    "NUM_EXT", "NUM_INT", "IDE_COL", "IDE_CP",
    "CVE_LOC_RES", "DES_LOC_RES", "CVE_MPO_RES", "DES_MPO_RES", "DES_JUR_RES",
    "CVE_EDO_RES", "DES_EDO_RES", "CVE_EDO_REP", "DES_EDO_REP",
    "CVE_LOC_REP", "DES_LOC_REP", "CVE_MPO_REP", "DES_MPO_REP",
    "ESTATUS_CASO", "CVE_DIAG_PROBABLE", "DES_DIAG_PROBABLE", "DES_DIAG_FINAL",
    "FEC_INI_SIGNOS_SINT", "ANO", "SEM", "MANEJO", "DES_INS_UNIDAD",
    "DENGUE_SER_TRIPLEX", "FEC_INGRESO",
]

#: columns of ``read_dengue_dataset()`` (same as densnv in R)
DENGUE_COLUMNS = [
    "FOL_ID", "VEC_ID", "IDE_EDA_ANO", "IDE_SEX", "DES_CAL", "IDE_CAL",
    "NUM_EXT", "NUM_INT", "IDE_COL", "IDE_CP",
    "CVE_LOC_RES", "DES_LOC_RES", "CVE_MPO_RES", "DES_MPO_RES", "DES_JUR_RES",
    "CVE_EDO_RES", "DES_EDO_RES", "CVE_LOC_REP", "DES_LOC_REP", "CVE_MPO_REP",
    "DES_MPO_REP", "DES_JUR_REP", "CVE_EDO_REP", "DES_EDO_REP",
    "ESTATUS_CASO", "CVE_DIAG_PROBABLE", "DES_DIAG_PROBABLE", "DES_DIAG_FINAL",
    "FEC_INI_SIGNOS_SINT", "ANO", "SEM", "MANEJO", "DES_INS_UNIDAD",
    "DENGUE_SER_TRIPLEX", "FEC_INGRESO",
]

# numeric columns (as integers); the geographic keys stay as zero-padded text
_INT_COLUMNS = ["IDE_EDA_ANO", "IDE_SEX", "ESTATUS_CASO", "CVE_DIAG_PROBABLE",
                "ANO", "SEM", "MANEJO", "DENGUE_SER_TRIPLEX"]
_KEY_WIDTH = {"CVE_EDO": 2, "CVE_MPO": 3, "CVE_LOC": 4}
_DATE_COLUMNS = ["FEC_INI_SIGNOS_SINT", "FEC_INGRESO"]
_FOREIGN = ["OTROS PAISES", "OTROS PAISES DE LATINOAMERICA",
            "ESTADOS UNIDOS DE NORTEAMERICA"]


def _sniff_sep(path: Path, encoding: str) -> str:
    with open(path, encoding=encoding, errors="replace") as f:
        header = f.readline()
    counts = {s: header.count(s) for s in ("|", "\t", ",", ";")}
    return max(counts, key=counts.get)


def _pad(s: pd.Series, width: int) -> pd.Series:
    s = s.astype("string").str.strip()
    s = s.str.replace(r"\.0$", "", regex=True)
    return s.str.zfill(width)


def _parse_dates(s: pd.Series) -> pd.Series:
    s = s.astype("string").str.strip()
    iso = pd.to_datetime(s, format="%Y-%m-%d", errors="coerce")
    dmy = pd.to_datetime(s, format="%d/%m/%Y", errors="coerce")
    return iso.fillna(dmy)


def _tidy(x: pd.DataFrame) -> pd.DataFrame:
    # like data.table::fread(), remove the spaces around the values
    for c in x.columns:
        x[c] = x[c].str.strip()
    for c in _INT_COLUMNS:
        if c in x:
            x[c] = pd.to_numeric(x[c], errors="coerce").astype("Int64")
    for c in x.columns:
        for prefix, width in _KEY_WIDTH.items():
            if c.startswith(prefix + "_"):
                x[c] = _pad(x[c], width)
    for c in _DATE_COLUMNS:
        if c in x:
            x[c] = _parse_dates(x[c])
    return x


def read_sinave(path: str | Path, columns: list[str] | None = None,
                encoding: str = "latin-1") -> pd.DataFrame:
    """Read one SINAVE txt file.

    Parameters
    ----------
    path : str or Path
        The SINAVE txt file (e.g. ``"DENGUE2_.txt"``).
    columns : list of str, optional
        Columns to keep. All the columns by default.
    encoding : str
        Encoding of the file; SINAVE files are Latin-1.

    Returns
    -------
    pandas.DataFrame
        One row per case. The geographic keys (``CVE_EDO_*``, ``CVE_MPO_*``,
        ``CVE_LOC_*``) are zero-padded text; ``ESTATUS_CASO``, ``ANO``,
        ``SEM``, ``MANEJO`` and ``DENGUE_SER_TRIPLEX`` are integers; the dates
        are datetimes.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"the file does not exist: {path}")
    sep = _sniff_sep(path, encoding)
    usecols = None
    if columns is not None:
        wanted = set(columns)
        usecols = lambda c: c.strip() in wanted  # noqa: E731
    x = pd.read_csv(path, sep=sep, encoding=encoding, dtype=str,
                    quoting=csv.QUOTE_NONE, usecols=usecols,
                    keep_default_na=False, na_values=[""],
                    on_bad_lines="warn", low_memory=False)
    x.columns = [c.strip() for c in x.columns]
    return _tidy(x)


def read(path: str | Path, vbd: str, complete: bool = True) -> pd.DataFrame:
    """Read the SINAVE datasets of a vector-borne disease.

    Parameters
    ----------
    path : str or Path
        Folder with the SINAVE txt files.
    vbd : str
        Text that identifies the files of the disease, e.g. ``"DENGUE"``.
    complete : bool
        ``True`` keeps all the columns; ``False`` keeps the columns of the
        arboviruses (see ``ARBO_COLUMNS``).

    Returns
    -------
    pandas.DataFrame
    """
    files = sorted(p for p in Path(path).glob("*.txt") if vbd in p.name)
    if not files:
        raise FileNotFoundError(f"no txt file with '{vbd}' in {path}")
    cols = None if complete else ARBO_COLUMNS
    parts = [read_sinave(f, columns=cols) for f in files]
    return parts[0] if len(parts) == 1 else pd.concat(parts, ignore_index=True)


def read_dengue_dataset(path: str | Path, spatial_resolution: str,
                        status_caso=(1, 2), des_edo_res=None,
                        des_mpo_res=None) -> pd.DataFrame:
    """Read the dengue dataset of SINAVE for the country, a state or a municipality.

    Parameters
    ----------
    path : str or Path
        The dengue txt file (``"DENGUE2_.txt"``).
    spatial_resolution : {"country", "state", "municipality"}
    status_caso : int or list of int
        1 probable, 2 confirmed, 3 discarded.
    des_edo_res : str or list of str, optional
        State(s) of residence in uppercase, as in SINAVE (``"NUEVO LEON"``).
    des_mpo_res : str or list of str, optional
        Municipality(ies) of residence in uppercase.

    Returns
    -------
    pandas.DataFrame
        The cases of residence in Mexico (foreign residences removed).
    """
    if spatial_resolution not in ("country", "state", "municipality"):
        raise ValueError("spatial_resolution must be 'country', 'state' or 'municipality'")
    x = read_sinave(path, columns=DENGUE_COLUMNS)
    x = x[~x["DES_EDO_RES"].isin(_FOREIGN)]
    x = x[x["ESTATUS_CASO"].isin(_as_list(status_caso))]
    if spatial_resolution in ("state", "municipality"):
        x = x[x["DES_EDO_RES"].isin(_as_list(des_edo_res))]
    if spatial_resolution == "municipality":
        x = x[x["DES_MPO_RES"].isin(_as_list(des_mpo_res))]
    return x.reset_index(drop=True)


def _as_list(v):
    if v is None:
        return []
    return list(v) if isinstance(v, (list, tuple, set, range)) else [v]
