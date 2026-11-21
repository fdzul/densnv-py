# Español

**densnv** es la versión en Python del paquete de R
[densnv](https://fdzul.github.io/densnv/) para el análisis exploratorio de las
bases de datos de Enfermedades Transmitidas por Vector del SINAVE. Es parte del
[dengueverse](https://fdzul.github.io/dengueverse/es/).

!!! note "En desarrollo"
    El módulo 1 (lectura y agregación de casos) está listo. Siguen el número
    reproductivo efectivo (Rt), los focos activos de transmisión y las gráficas.

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
```

Las claves geográficas son texto con ceros a la izquierda (`"31"`, `"050"`,
`"0001"`) y los nombres conservan sus acentos.

## Datos

Las bases del SINAVE son para uso exclusivo del personal de las Secretarías de
Salud estatales y del nivel federal. El paquete las lee de tu computadora;
nunca las descarga ni las publica.
