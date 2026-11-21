# Español

**densnv** es la versión en Python del paquete de R
[densnv](https://fdzul.github.io/densnv/) para el análisis exploratorio de las
bases de datos de Enfermedades Transmitidas por Vector del SINAVE. Es parte del
[dengueverse](https://fdzul.github.io/dengueverse/es/).

!!! note "En desarrollo"
    Listos: la lectura y agregación de casos (módulo 1), el número
    reproductivo efectivo, Rt (módulo 2), y los focos activos de transmisión
    con el diagnóstico del escenario epidemiológico (módulo 3). Siguen las
    gráficas.

## Instalación

```bash
pip install "git+https://github.com/fdzul/densnv-py"
```

## Ejemplo

```python
import densnv

x   = densnv.read("~/sinave/2026", vbd="DENGUE")
mpo = densnv.cases_by_agem("~/sinave/2026", "DENGUE", status_case=[1, 2])
loc = densnv.cases_by_loc("~/sinave/2026/DENGUE2_.txt", "DENV")
r   = densnv.cases_week("~/sinave/2026", mpo="MERIDA", year=2026)

# Rt de la localidad de Mérida, ventana de 4 semanas
dataset = densnv.read_sinave("~/sinave/2026/DENGUE2_.txt")
rt = densnv.rt_dengue(dataset, level="locality", cve_edo="31",
                      cve_mpo="050", cve_loc="0001", window_days=28)
rt["rt_ventana"]

# focos activos y escenario epidemiológico (casos geocodificados z)
res = densnv.mp_active_foci(z, cve_edo="31", locality="Merida",
                            temporal_scale=28, buffer_m=0, method="st-dbscan",
                            n_sim=999, seed=123, print_map=True)
res["diagnostico"]["tabla"]
```

Las claves geográficas son texto con ceros a la izquierda (`"31"`, `"050"`,
`"0001"`) y los nombres conservan sus acentos.

## Datos

Las bases del SINAVE son para uso exclusivo del personal de las Secretarías de
Salud estatales y del nivel federal. El paquete las lee de tu computadora;
nunca las descarga ni las publica.
