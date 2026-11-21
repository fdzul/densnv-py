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

## Effective reproductive number (Rt)

```python
dataset = densnv.read_sinave("~/sinave/2026/DENGUE2_.txt")

rt = densnv.rt_dengue(dataset, level="locality", place="REP",
                      cve_edo="31", cve_mpo="050", cve_loc="0001",
                      window_days=28)
rt["rt_ventana"]    # Rt of the last 28 days, 95% CrI and phase
rt["rt"]            # Rt curve (weekly sliding windows)
rt["figure"]        # matplotlib figure
```

## Active transmission foci and epidemiological scenario

```python
# geocoded cases: columns long, lat, ESTATUS_CASO, FEC_INI_SIGNOS_SINT, CVE_*_REP
res = densnv.mp_active_foci(z, cve_edo="31", locality="Merida",
                            temporal_scale=28, buffer_m=0,
                            method="st-dbscan", n_sim=999, seed=123,
                            print_map=True)
res["foci"]                      # active foci, Knox test by focus
res["blocks"]                    # blocks to treat
res["diagnostico"]["tabla"]      # scenario, criteria and action
res["map"]                       # interactive map (pip install folium)
```

## Graphs and maps

Each function returns the table behind the figure (`data`) and the figure.

```bash
pip install "densnv[plots] @ git+https://github.com/fdzul/densnv-py"   # plotly, squarify, folium
```

```python
ds = densnv.read_sinave("~/sinave/2026/DENGUE2_.txt")

densnv.plt_heatmap(ds, state="YUCATAN", year=2026, status=[1, 2])["figure"]   # plotly
densnv.plt_pyramid(ds, title="Casos confirmados")["figure"]                   # matplotlib
densnv.plt_serotype(ds, year=2026)["figure"]
densnv.plt_treemap(ds, country=True, year=2026)["figure"]
densnv.mp_bumpmap(ds, year=2026)["figure"]

loc = densnv.cases_by_loc("~/sinave/2026/DENGUE2_.txt", "DENV")
densnv.mp_cases_loc(loc, var="probable", top_n=10)["map"]                    # folium
densnv.mp_serotype("~/sinave/2026/DENGUE2_.txt", cve_edo="31")["map"]
```
