# Getting Started

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install "git+https://github.com/fdzul/densnv-py"
```

## Read the SINAVE files

`read()` reads every txt file of a folder whose name contains `vbd`:

```python
import densnv

x = densnv.read("~/sinave/2026", vbd="DENGUE")            # all the columns
x = densnv.read("~/sinave/2026", vbd="DENGUE", complete=False)  # arbovirus columns
```

The geographic keys are zero-padded text (`CVE_EDO_REP == "31"`,
`CVE_MPO_REP == "050"`, `CVE_LOC_REP == "0001"`), the status, year and week are
integers and the dates are datetimes. The names keep their accents
(`"YUCATÁN"`).

## Cases by municipality and locality

```python
mpo = densnv.cases_by_agem("~/sinave/2026", "DENGUE", status_case=[1, 2],
                           country=False, cve_edo="31")
mpo.plot(column="n", legend=True)

loc = densnv.cases_by_loc("~/sinave/2026/DENGUE2_.txt", "DENV", week=range(1, 11))
loc.sort_values("probable", ascending=False).head()
```

## Cases by week

```python
r = densnv.cases_week("~/sinave/2026", mpo="MERIDA", year=2026)
r["estatus"]        # probable, confirmed and estimated by week
r["figure"]         # matplotlib figure
```
