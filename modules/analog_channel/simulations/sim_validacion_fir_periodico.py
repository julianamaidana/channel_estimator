"""
Tests de la planta FIR periodica de 16 fases.

Correr con:  python sim_validacion_fir_periodico.py

Cada test imprime una tabla con la fase, la entrada, la salida y lo esperado.
Si un test falla, el programa se corta en ese test.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))  # para importar fir_periodico.py desde el subdirectorio "python"

import numpy as np
from fir_periodico import planta_fir_periodica, planta_fir_periodica_bloque, N_FASES # type: ignore

N = 100   # largo de la senal de prueba
L = 11    # coeficientes por kernel

rng = np.random.default_rng(0)   # semilla fija: siempre salen los mismos numeros
np.set_printoptions(precision=3, suppress=True)


def mostrar(x, r, esperado, desde=0, hasta=12):
    """Imprime muestra por muestra: indice, fase, entrada, salida y valor esperado."""
    print("   n  fase     x[n]     r[n]  esperado")
    for n in range(desde, min(hasta, len(x))):
        print(f"{n:4d} {n % N_FASES:5d} {x[n]:8.3f} {r[n]:8.3f} {esperado[n]:9.3f}")


def kernels_iguales(kernel):
    """Matriz h con el mismo kernel en las 16 fases."""
    h = np.zeros((N_FASES, len(kernel)))
    h[:] = kernel
    return h


def procesar_en_bloques(x, h, tamano):
    """Procesa x en bloques de 'tamano' muestras y junta las salidas."""
    L = h.shape[1]
    memoria = np.zeros(L - 1)
    n = 0
    salidas = []
    for inicio in range(0, len(x), tamano):
        bloque = x[inicio:inicio + tamano]
        r_bloque, memoria, n = planta_fir_periodica_bloque(bloque, h, memoria, n)
        salidas.append(r_bloque)
    return np.concatenate(salidas)


# -----------------------------------------------------------------------------
# 1. Entrada nula
# -----------------------------------------------------------------------------
def test_1_entrada_nula():
    print("\n=== Test 1: entrada nula -> salida nula ===")
    h = rng.standard_normal((N_FASES, L))   # 16 kernels cualquiera
    x = np.zeros(N)
    r = planta_fir_periodica(x, h)

    mostrar(x, r, esperado=np.zeros(N), hasta=6)
    assert np.all(r == 0)
    print("OK")


# -----------------------------------------------------------------------------
# 2. Canal identidad
# -----------------------------------------------------------------------------
def test_2_canal_identidad():
    print("\n=== Test 2: canal identidad -> salida = entrada ===")
    kernel = np.zeros(L)
    kernel[0] = 1
    h = kernels_iguales(kernel)
    print("kernel en las 16 fases:", kernel)

    x_impulso = np.zeros(N)
    x_impulso[0] = 1
    x_aleatoria = rng.standard_normal(N)

    for x in (x_impulso, x_aleatoria):
        r = planta_fir_periodica(x, h)
        mostrar(x, r, esperado=x, hasta=6)
        assert np.allclose(r, x)
    print("OK")


# -----------------------------------------------------------------------------
# 3. Ganancia comun
# -----------------------------------------------------------------------------
def test_3_ganancia_comun():
    print("\n=== Test 3: ganancia comun -> salida = 2 x entrada ===")
    kernel = np.zeros(L)
    kernel[0] = 2
    h = kernels_iguales(kernel)
    print("kernel en las 16 fases:", kernel)

    x = rng.standard_normal(N)
    r = planta_fir_periodica(x, h)

    mostrar(x, r, esperado=2 * x, hasta=6)
    assert np.allclose(r, 2 * x)
    print("OK")


# -----------------------------------------------------------------------------
# 4. FIR comun con memoria
# -----------------------------------------------------------------------------
def test_4_fir_comun():
    print("\n=== Test 4: FIR comun, kernel [0.6, 0.3, 0.1] en las 16 fases ===")
    kernel = [0.6, 0.3, 0.1]
    h = kernels_iguales(kernel)

    # Con un impulso, la salida es el kernel
    x = np.array([1.0, 0, 0, 0])
    r = planta_fir_periodica(x, h)
    esperado = [0.6, 0.3, 0.1, 0.0]
    mostrar(x, r, esperado)
    assert np.allclose(r, esperado)

    # Con una secuencia aleatoria, da lo mismo que una convolucion comun
    x = rng.standard_normal(N)
    r = planta_fir_periodica(x, h)
    referencia = np.convolve(x, kernel)[:N]
    mostrar(x, r, referencia, hasta=6)
    assert np.allclose(r, referencia)
    print("OK")


# -----------------------------------------------------------------------------
# 5. Kernels distintos con memoria
# -----------------------------------------------------------------------------
def test_5_kernels_distintos():
    print("\n=== Test 5: kernels distintos, un impulso en cada fase ===")

    # Con un impulso en m, la salida es r[m+k] = h[(m+k) % 16, k]:
    # cada muestra de salida usa el kernel de SU fase (una diagonal de h).

    # Un caso impreso: L = 3, impulso en n = 16 (fase 0)
    h = rng.standard_normal((N_FASES, 3))
    x = np.zeros(20)
    x[16] = 1
    r = planta_fir_periodica(x, h)
    print("fila 0 de h:", h[0])
    for k in range(3):
        print(f"r[{16 + k}] = {r[16 + k]:6.3f}   h[{k},{k}] = {h[k, k]:6.3f}")
    print("La salida debería ser la diagonal de h.")

    # Todos los casos: impulsos en n = 16 ... 31 (las 16 fases), con varios L
    for L5 in (1, 3, 11, 20):
        h = rng.standard_normal((N_FASES, L5))
        for m in range(16, 32):
            x = np.zeros(m + L5 + 1)
            x[m] = 1
            r = planta_fir_periodica(x, h)

            esperado = np.zeros(len(x))
            for k in range(L5):
                esperado[m + k] = h[(m + k) % N_FASES, k]
            assert np.allclose(r, esperado)
    print("OK")


# -----------------------------------------------------------------------------
# 6. Kernels distintos con entrada aleatoria
# -----------------------------------------------------------------------------
def test_6_entrada_aleatoria():
    print("\n=== Test 6: kernels distintos y entrada aleatoria ===")
    h = rng.standard_normal((N_FASES, L))
    x = rng.standard_normal(N)
    r = planta_fir_periodica(x, h)

    # Referencia: la convolucion con cada uno de los 16 kernels,
    # y en cada muestra me quedo con la del kernel de su fase
    convoluciones = []
    for p in range(N_FASES):
        convoluciones.append(np.convolve(x, h[p])[:N])

    referencia = np.zeros(N)
    for n in range(N):
        referencia[n] = convoluciones[n % N_FASES][n]

    mostrar(x, r, referencia, desde=12, hasta=20)   # se ve la vuelta de la fase 15 a la 0
    assert np.allclose(r, referencia)
    print("OK")


# -----------------------------------------------------------------------------
# 7. Procesamiento por bloques
# -----------------------------------------------------------------------------
def test_7_bloques():
    print("\n=== Test 7: procesamiento por bloques = senal completa ===")

    # Ejemplo: L = 3, una senal de 20 muestras partida en dos bloques de 10
    h = rng.standard_normal((N_FASES, 3))
    x = rng.standard_normal(20)
    r_completa = planta_fir_periodica(x, h)

    memoria = np.zeros(2)   # L-1 = 2 ceros al principio
    n = 0
    r1, memoria, n = planta_fir_periodica_bloque(x[:10], h, memoria, n)
    print("Despues del bloque 1: memoria =", memoria, "(= x[8], x[9]), n =", n)
    r2, memoria, n = planta_fir_periodica_bloque(x[10:], h, memoria, n)
    r_bloques = np.concatenate([r1, r2])

    print("r[n]: en bloques, esperado: senal completa (el bloque 2 empieza en n = 10)")
    mostrar(x, r_bloques, r_completa, desde=7, hasta=14)
    assert np.allclose(r_bloques, r_completa)

    # Si el bloque 2 arranca sin la memoria ni el indice del bloque 1, da mal
    r2_mal, _, _ = planta_fir_periodica_bloque(x[10:], h, np.zeros(2), 0)
    print(f"Bloque 2 sin memoria ni indice: r[10] = {r2_mal[0]:.3f} en vez de {r_completa[10]:.3f}")
    assert not np.allclose(r2_mal, r_completa[10:])

    # Todos los casos: varios L y tamanos de bloque (103 no es multiplo de 16)
    for L7 in (1, 3, 11, 20):
        h = rng.standard_normal((N_FASES, L7))
        x = rng.standard_normal(103)
        r_completa = planta_fir_periodica(x, h)
        for tamano in (1, 2, 10, 16, 103):
            r_bloques = procesar_en_bloques(x, h, tamano)
            assert np.allclose(r_bloques, r_completa)
    print("OK")


if __name__ == "__main__":
    test_1_entrada_nula()
    test_2_canal_identidad()
    test_3_ganancia_comun()
    test_4_fir_comun()
    test_5_kernels_distintos()
    test_6_entrada_aleatoria()
    test_7_bloques()
    print("\nTodos los tests pasaron.")