# densnv (Python)

**Exploratory analysis of the SINAVE vector-borne disease datasets of Mexico.**
Python version of the R package [densnv](https://github.com/fdzul/densnv), part
of the [dengueverse](https://fdzul.github.io/dengueverse/).

> Status: in development. Module 1 (reading and aggregating the cases) is
> ready; Rt, active transmission foci and the graphs come next.

## Installation

```bash
pip install "git+https://github.com/fdzul/densnv-py"
```

It installs [rgeomex](https://fdzul.github.io/rgeomex-py/) for the
municipalities and localities of Mexico.

## Example

```python
import densnv

# read the SINAVE files of a disease (folder with the txt files)
x = densnv.read("path/to/sinave", vbd="DENGUE")

# probable and confirmed cases of residence in a state
nl = densnv.read_dengue_dataset("DENGUE2_.txt", "state",
                                status_caso=[1, 2], des_edo_res="NUEVO LEON")

# cases by municipality (polygons) and by locality (points)
mpo = densnv.cases_by_agem("path/to/sinave", "DENGUE", status_case=[1, 2])
loc = densnv.cases_by_loc("DENGUE2_.txt", "DENV", week=range(1, 11))

# cases by week in a municipality, with the figure
r = densnv.cases_week("path/to/sinave", mpo="MERIDA", year=2026)
r["figure"].savefig("merida.png")
```

## Data

The SINAVE datasets are for the exclusive use of the staff of the State
Ministries of Health and of the federal level. The package reads them from your
computer; it never downloads or publishes them. The tests use a synthetic file.

## License

MIT. See [LICENSE](LICENSE).
