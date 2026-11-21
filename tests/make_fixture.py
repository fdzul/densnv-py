"""Build the synthetic SINAVE file used by the tests (no real cases)."""
from pathlib import Path

from densnv._read import DENGUE_COLUMNS

cols = DENGUE_COLUMNS + ["OTRA_COLUMNA"]

def row(**kw):
    base = dict(FOL_ID="1", VEC_ID="1", IDE_EDA_ANO="30", IDE_SEX="1", DES_CAL="CALLE",
                IDE_CAL="X", NUM_EXT="1", NUM_INT="", IDE_COL="CENTRO", IDE_CP="97000",
                CVE_LOC_RES="1", DES_LOC_RES="MÉRIDA", CVE_MPO_RES="50", DES_MPO_RES="MÉRIDA",
                DES_JUR_RES="MÉRIDA", CVE_EDO_RES="31", DES_EDO_RES="YUCATÁN",
                CVE_LOC_REP="1", DES_LOC_REP="MÉRIDA", CVE_MPO_REP="50", DES_MPO_REP="MÉRIDA",
                DES_JUR_REP="MÉRIDA", CVE_EDO_REP="31", DES_EDO_REP="YUCATÁN",
                ESTATUS_CASO="1", CVE_DIAG_PROBABLE="1", DES_DIAG_PROBABLE="DENGUE NO GRAVE",
                DES_DIAG_FINAL="DENGUE NO GRAVE", FEC_INI_SIGNOS_SINT="2026-03-02", ANO="2026",
                SEM="10", MANEJO="2", DES_INS_UNIDAD="SSA", DENGUE_SER_TRIPLEX="",
                FEC_INGRESO="2026-03-04", OTRA_COLUMNA="z")
    base.update(kw)
    return [base[c] for c in cols]

rows = [
    # Merida, 2026: 4 probable, 3 confirmed (DENV2, DENV2, DENV3), 2 discarded
    row(), row(), row(SEM="11"), row(SEM="11", MANEJO="1"),
    row(ESTATUS_CASO="2", DENGUE_SER_TRIPLEX="2", DES_DIAG_FINAL="DENGUE CON SIGNOS DE ALARMA"),
    row(ESTATUS_CASO="2", DENGUE_SER_TRIPLEX="2", SEM="11", MANEJO="1",
        DES_DIAG_FINAL="DENGUE GRAVE"),
    row(ESTATUS_CASO="2", DENGUE_SER_TRIPLEX="3", SEM="11"),
    row(ESTATUS_CASO="3"), row(ESTATUS_CASO="3", SEM="11"),
    # Kanasin, 2026: 1 confirmed DENV1
    row(CVE_MPO_REP="41", DES_MPO_REP="KANASÍN", CVE_MPO_RES="41", DES_MPO_RES="KANASÍN",
        ESTATUS_CASO="2", DENGUE_SER_TRIPLEX="1"),
    # Merida 2025 (other year), date in dd/mm/yyyy
    row(ANO="2025", FEC_INI_SIGNOS_SINT="15/07/2025"),
    # residence abroad, reported in Merida
    row(DES_EDO_RES="ESTADOS UNIDOS DE NORTEAMERICA", CVE_EDO_RES="99"),
    # residence in "other countries" written with a trailing space
    row(DES_EDO_RES="OTROS PAISES ", CVE_EDO_RES="35", ESTATUS_CASO="3", ANO="2025"),
    # missing locality of report
    row(CVE_LOC_REP="0"),
]

out = Path(__file__).parent / "data" / "DENGUE2_.txt"
out.parent.mkdir(exist_ok=True)
with open(out, "w", encoding="latin-1", newline="") as f:
    f.write("|".join(cols) + "\r\n")
    for r in rows:
        f.write("|".join(r) + "\r\n")
print(out)
