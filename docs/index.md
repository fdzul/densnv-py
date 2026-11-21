# densnv

**Exploratory analysis of the SINAVE vector-borne disease datasets of Mexico**,
in Python. It is the Python version of the R package
[densnv](https://fdzul.github.io/densnv/), part of the
[dengueverse](https://fdzul.github.io/dengueverse/).

!!! note "In development"
    Ready: reading and aggregating the cases (module 1) and the effective
    reproductive number, Rt (module 2). The active transmission foci and the
    graphs come next.

```python
import densnv

loc = densnv.cases_by_loc("DENGUE2_.txt", "DENV")
loc.explore(column="probable")
```

## Functions

1. **Read:** `read()`, `read_sinave()`, `read_dengue_dataset()`.
2. **Aggregate:** `cases_by_agem()` (municipalities), `cases_by_loc()`
   (localities), `cases_week()` (weeks).
3. **Rt:** `rt_dengue()` for the country, a state, a municipality or a
   locality (Cori et al. 2013, the same estimates as EpiEstim).

The results are the same as in R: every function is validated against densnv
with the real SINAVE dataset.

## Data

The SINAVE datasets are for the exclusive use of the staff of the State
Ministries of Health and of the federal level. The package reads them from your
computer; it never downloads or publishes them.
