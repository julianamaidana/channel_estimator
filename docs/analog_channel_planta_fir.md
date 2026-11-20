# Planta FIR periódica de 16 fases

Issue asociado: #3 

Modelo en Python del canal equivalente del camino analógico del FI-DAC.

## Qué hace

Es un FIR cuyos coeficientes cambian según la fase de cada muestra. El patrón se repite cada 16 muestras:

```
r[n] = h[p,0]·x[n] + h[p,1]·x[n-1] + ... + h[p,L-1]·x[n-L+1],   con p = n mod 16
```

- `h` es una matriz de 16 filas por `L` columnas. La fila `p` es el kernel de la fase `p`.
- Para cada muestra de salida se usa el kernel de su fase, `p = n mod 16`.
- La primera muestra usa la fase 0, y antes del comienzo de la señal se toman ceros.
- Canal identidad: los 16 kernels iguales a `[1, 0, ..., 0]`.

El 16 sale del entrelazado: sobre la grilla de 72 GS/s, el patrón se repite cada 72 / 4,5 GHz = 16 muestras. En la teoría de Bonetti et al. (2023), es un sistema lineal periódicamente variante (PTV) en tiempo discreto (Eq. 36).

## Archivos

```text
modules/analog_channel/
├── python/
│   └── fir_periodico.py                  la planta
└── simulations/
    └── sim_validacion_fir_periodico.py   los tests
```

Hace falta Python 3 y NumPy (`pip install -r requirements.txt` desde la raíz del repositorio).

## Cómo se usa

Los ejemplos se corren desde `modules/analog_channel/python/`.

Con la señal completa:

```python
import numpy as np
from fir_periodico import planta_fir_periodica, N_FASES

h = np.zeros((N_FASES, 3))
h[:] = [0.6, 0.3, 0.1]          # el mismo kernel en las 16 fases

r = planta_fir_periodica([1, 0, 0, 0], h)
print(r)                         # [0.6 0.3 0.1 0. ]
```

Por bloques, como llega la señal en el sistema real. Entre un bloque y el siguiente se pasan dos cosas:

- `memoria`: las últimas `L-1` muestras del bloque anterior (al principio, ceros).
- `n`: el índice global de la próxima muestra (al principio, 0), para que la fase no vuelva a 0 en cada bloque.

```python
import numpy as np
from fir_periodico import planta_fir_periodica, planta_fir_periodica_bloque, N_FASES

h = np.random.default_rng(0).standard_normal((N_FASES, 3))   # L = 3
x = np.random.default_rng(1).standard_normal(20)

memoria = np.zeros(2)            # L-1 = 2 ceros
n = 0
r1, memoria, n = planta_fir_periodica_bloque(x[:10], h, memoria, n)
r2, memoria, n = planta_fir_periodica_bloque(x[10:], h, memoria, n)

r = np.concatenate([r1, r2])
print(np.allclose(r, planta_fir_periodica(x, h)))   # True
```

Después del primer bloque, `memoria` tiene `x[8]` y `x[9]`, y `n = 10`.

## Tests

Desde la raíz del repositorio:

```
python modules/analog_channel/simulations/sim_validacion_fir_periodico.py
```

Cada test imprime una tabla con la fase, la entrada, la salida y lo esperado. Si alguno falla, el programa se corta en ese test.

1. **Entrada nula:** con 16 kernels cualquiera, la salida es cero.
2. **Canal identidad:** la salida es igual a la entrada, con un impulso y con una secuencia aleatoria.
3. **Ganancia común:** con kernels `[2, 0, ..., 0]`, la salida es el doble de la entrada.
4. **FIR común:** con el kernel `[0.6, 0.3, 0.1]` en las 16 fases, un impulso devuelve el kernel y una secuencia aleatoria da lo mismo que `np.convolve`.
5. **Kernels distintos:** un impulso en cada una de las 16 fases. Cada muestra de salida usa el kernel de su fase, así que la salida recorre una diagonal de `h`: `r[m+k] = h[(m+k) mod 16, k]`.
6. **Entrada aleatoria:** con kernels distintos, comparada contra la convolución con cada kernel, tomando en cada muestra la de su fase.
7. **Bloques:** procesar en bloques da lo mismo que la señal completa, con bloques de 1, 2, 10, 16 y 103 muestras y `L` = 1, 3, 11 y 20. También muestra que, sin pasar la memoria y el índice, el resultado da mal.

## Pendiente

- Efectos del camino analógico: ancho de banda (cambiando los kernels) y offset (un valor por fase sumado a la salida).
- La implementación es un doble loop: clara, pero lenta para señales muy largas.

## Referencia

J. I. Bonetti, A. Galetto y M. R. Hueda, *Theory of Periodically Time-Variant Linear Systems*, arXiv:2305.07556, 2023.