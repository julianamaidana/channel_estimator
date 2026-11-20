"""
sim_no_idealidades_bw_offs.py — Validación del modelo de ancho de banda y offset (issue #4).

Prueba las funciones de no_idealidades.py. Cada test arma un caso con
resultado conocido y lo compara con un assert. Los parámetros (B, orden,
desajuste, semilla, L) están definidos al principio del archivo.

Uso:
    python sim_no_idealidades_bw_offs.py

Tests de offset
    O1. Entrada nula: la salida es el offset de su fase, u[n % 16].
    O2. Entrada aleatoria y kernels Butterworth: la diferencia entre la salida
        con offset y sin offset es u[n % 16]. El offset no toca la parte lineal.
    O3. random_offset: los 16 valores están en [-lsb, +lsb], y con la misma
        semilla sale el mismo u.

Tests de ancho de banda
    B1. Canal identidad (un 1 en la columna 0 de cada fila): la salida es igual
        a la entrada.
    B2. Con L = 200, la suma de los taps de cada fila da 1 (ganancia en DC de
        un pasabajos).
    B3. Con L = 200, en cada fase |H(B_fases[p])| = -3,01 dB ± 0,05: el corte
        cae en el ancho de banda sorteado.
    B4. Los B_fases están entre B·(1 - desajuste) y B·(1 + desajuste), y las
        filas de h son distintas entre sí.
    B5. Con desajuste = 0, la planta da lo mismo que np.convolve: sin desajuste
        es un filtro común.
    B6. Un corte por encima de Nyquist (30 GHz con ±25 % > 36 GHz) levanta
        ValueError.

Efecto de L (mediciones, sin assert)
    Las dos se hacen en el peor caso: el kernel que corta en B·(1 + desajuste),
    que es el de respuesta más larga. El L de los tests sale de largo_minimo()
    con ese mismo corte (L = 74 con ±15 %).
    L1. Energía perdida por truncamiento para cada L del barrido. Umbral:
        -44 dB (ruido de un DAC con ENOB de 7 bits).
    L2. Diferencia máxima entre |H| del kernel truncado y |H| del kernel largo
        (L = 400), hasta 30 GHz.
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import matplotlib.pyplot as plt

import numpy as np
from fir_periodico import N_FASES # type: ignore
from no_idealidades import planta_con_offset, kernel_butterworth, random_offset, random_butterworth, largo_minimo # type: ignore

import scipy.signal as signal

FS = 72e9
BITS = 8
BITS_ENOB = 7
RANGO = 1.0     # rango de la senal de entrada (y del offset)
B = 30e9
orden = 4
semilla_o = 0
semilla_x = 1
semilla_bw = 2
desajuste = 0.15   # desajuste relativo del ancho de banda (0.15 = ±15 %)
UMBRAL_DB = -(6.02 * BITS_ENOB + 1.76)                   # -44 dB: ruido de un DAC con ENOB de 7 bits
L = largo_minimo(B * (1 + desajuste), orden, UMBRAL_DB)  # peor caso: la fase que corta más alto (L = 74 con ±15 %)
L_BARRIDO = [5, 11, 21, 31, 41, 51, 61, 71, 81, 91, 101]
    
def test_O1():
    #-------------------------------------------------------------
    # Test con entrada nula
    #-------------------------------------------------------------
    print("\n=== Test O1: entrada nula -> salida = offset de su fase ===")
    lsb = RANGO / (2 ** BITS)  # tamaño del LSB
    u = random_offset(lsb, semilla_o)
    print(f"Offset por fase (lsb={lsb:.3f}):\n {u} \n")
    x = np.zeros(100)
    h = np.zeros((N_FASES, L))
    r = planta_con_offset(x, h, u)
    print(f"Salida con offset:\n {r[0:16]}")
    print("--------------------------------")
    r_esperado = u[np.arange(len(x)) % N_FASES]
    assert np.allclose(r, r_esperado)

def test_O2():
    #-------------------------------------------------------------
    # Test con entrada aleatoria y kernels Butterworth
    #-------------------------------------------------------------
    print("\n=== Test O2: entrada aleatoria, kernels Butterworth ===")
    lsb = RANGO / (2 ** BITS)  # tamaño del LSB
    u = random_offset(lsb, semilla_o)
    print(f"Offset por fase (lsb={lsb:.3f}):\n {u} \n")
    x = np.random.default_rng(semilla_x).uniform(-RANGO/2, RANGO/2, 100)
    h = np.zeros((N_FASES, L))
    for p in range(N_FASES):
        h[p, :] = kernel_butterworth(B, orden, L)
    r_sin_offset = planta_con_offset(x, h, np.zeros(N_FASES))
    r_con_offset = planta_con_offset(x, h, u)
    print(f"Salida sin offset:\n {r_sin_offset[0:16]}\n")
    print(f"Salida con offset:\n {r_con_offset[0:16]}\n")
    diferencia = r_con_offset - r_sin_offset
    r_esperado = u[np.arange(len(x)) % N_FASES]
    assert np.allclose(diferencia, r_esperado)
    print("--------------------------------")

def test_O3():
    #-------------------------------------------------------------
    # Test random_offset: todos los valores en [-lsb, +lsb] y con la misma semilla sale el mismo u
    #-------------------------------------------------------------
    print("\n=== Test O3: random_offset ===")
    lsb = RANGO / (2 ** BITS)  # tamaño del LSB
    u1 = random_offset(lsb, semilla_o)
    u2 = random_offset(lsb, semilla_o)
    print(f"Offset por fase (lsb={lsb:.3f}):\n {u1} \n")
    assert np.all(u1 >= -lsb) and np.all(u1 <= lsb)
    assert np.allclose(u1, u2)
    print("--------------------------------")

def test_B1():
    #-------------------------------------------------------------
    # Test con canal identidad (un 1 en la columna 0 de cada fila)
    #-------------------------------------------------------------
    print("\n=== Test B1: canal identidad ===")
    h = np.zeros((N_FASES, L))
    h[:, 0] = 1.0
    x = np.random.default_rng(semilla_x).uniform(-RANGO/2, RANGO/2, 100)
    r = planta_con_offset(x, h, np.zeros(N_FASES))
    print(f"Entrada:\n {x[0:16]}\n")
    print(f"Salida:\n {r[0:16]}\n")
    assert np.allclose(r, x)
    print("--------------------------------")

def test_B2():
    #-------------------------------------------------------------
    # Test con L = 200: la suma de los taps de cada fila da 1
    #-------------------------------------------------------------
    print("\n=== Test B2: suma de los taps de cada fila da 1 ===")
    h, B_fases = random_butterworth(B, orden, 200, desajuste, semilla_bw)
    suma_taps = np.sum(h, axis=1)
    print(f"Suma de los taps por fase:\n {suma_taps}\n")
    assert np.allclose(suma_taps, np.ones(N_FASES))
    print("--------------------------------")

def test_B3():
    #-------------------------------------------------------------
    # Test con L = 200: en cada fase |H(B_fases[p])| = -3,01 dB ± 0,05
    #-------------------------------------------------------------
    print("\n=== Test B3: corte en B_fases ===")
    h, B_fases = random_butterworth(B, orden, 200, desajuste, semilla_bw)
    for p in range(N_FASES):
        w, H = signal.freqz(h[p, :], 1, worN=[B_fases[p]], fs=FS)
        modulo = abs(H[0])
        modulo_db = 20 * np.log10(modulo)
        print(f"Fase {p:2d}: B={B_fases[p]/1e9:.3f} GHz  |H(B)|={modulo:.4f}  ->  {modulo_db:.2f} dB")
        assert abs(modulo_db - (-3.01)) < 0.05, f"Fase {p}: el corte no cae en B"
    print("--------------------------------")
    
def test_B4():
    #-------------------------------------------------------------
    # Test: B_fases entre B·(1-d) y B·(1+d), filas distintas entre sí
    #-------------------------------------------------------------
    print("\n=== Test B4: B_fases entre B·(1-d) y B·(1+d), filas distintas ===")
    h, B_fases = random_butterworth(B, orden, 32, desajuste, semilla_bw)
    print(f"Ancho de banda sorteado por fase (B={B/1e9:.1f} GHz, desajuste={desajuste*100:.0f} %):\n {B_fases/1e9} GHz")
    assert np.all(B_fases >= B * (1 - desajuste)) and np.all(B_fases <= B * (1 + desajuste))
    for p in range(N_FASES):
        for q in range(p + 1, N_FASES):
            assert not np.allclose(h[p, :], h[q, :]), f"Fase {p} y {q}: kernels iguales"
    print("--------------------------------")

def test_B5():
    #-------------------------------------------------------------
    # Test: desajuste = 0 -> planta == np.convolve(x, h[0])[:len(x)]
    #-------------------------------------------------------------
    print("\n=== Test B5: desajuste = 0 -> planta == np.convolve(x, h[0])[:len(x)] ===")
    h, B_fases = random_butterworth(B, orden, 32, 0.0, semilla_bw)
    x = np.random.default_rng(semilla_x).uniform(-RANGO/2, RANGO/2, 100)
    r_planta = planta_con_offset(x, h, np.zeros(N_FASES))
    r_convolve = np.convolve(x, h[0, :])[:len(x)]
    assert np.allclose(r_planta, r_convolve)
    print("OK")
    print("--------------------------------")

def test_B6():
    #-------------------------------------------------------------
    # Test: corte > Nyquist (30 GHz con ±25 % > 36 GHz) levanta ValueError
    #-------------------------------------------------------------
    print("\n=== Test B6: corte > Nyquist levanta ValueError ===")
    try:
        h, B_fases = random_butterworth(30e9, 4, 32, 0.25, semilla_bw)
    except ValueError as e:
        print(f"Se levantó ValueError como se esperaba: {e}")
    else:
        raise AssertionError("No se levantó ValueError cuando el corte es mayor a Nyquist.")
    print("--------------------------------")

def test_L1():
    #-------------------------------------------------------------
    # Test: energía perdida en el peor caso para cada L, con umbral de -44 dB
    #-------------------------------------------------------------
    print("\n=== Test L1: energía perdida por truncamiento (peor caso) ===")
    B_peor = B * (1 + desajuste)                      # la fase que corta más alto es la de respuesta más larga
    h_larga = kernel_butterworth(B_peor, orden, 400)
    E_total = np.sum(h_larga ** 2)
    perdidas_db = []
    for l_test in L_BARRIDO:
        E_perdida = np.sum(h_larga[l_test:] ** 2)     # energía de la cola que se tira
        E_perdida_db = 10 * np.log10(E_perdida / E_total)
        perdidas_db.append(E_perdida_db)
        cumple = "cumple" if E_perdida_db < UMBRAL_DB else "no cumple"
        print(f"L={l_test:3d}: energía perdida = {E_perdida_db:7.2f} dB  ({cumple})")
    print(f"\nL mínimo con corte en B·(1 + desajuste) = {B_peor/1e9:.1f} GHz: L = {L}")
    print("--------------------------------")
    return perdidas_db

def test_L2():
    #-------------------------------------------------------------
    # Test: diferencia máxima entre |H| del FIR cortado y |H| del largo (L = 400), hasta 30 GHz
    #-------------------------------------------------------------
    print("\n=== Test L2: diferencia máxima entre |H| del FIR cortado y |H| del largo (peor caso) ===")
    B_peor = B * (1 + desajuste)
    h_larga = kernel_butterworth(B_peor, orden, 400)
    f = np.linspace(0, 30e9, 1000)   # frecuencias donde se compara: de 0 a 30 GHz

    # |H| del kernel largo, en dB (se calcula una sola vez)
    w, H = signal.freqz(h_larga, 1, worN=f, fs=FS)
    H_larga_db = 20 * np.log10(np.abs(H))

    diferencias_db = []
    for l_test in L_BARRIDO:
        w, H = signal.freqz(h_larga[:l_test], 1, worN=f, fs=FS)   # kernel cortado en l_test
        H_cortada_db = 20 * np.log10(np.abs(H))
        diferencia_db = np.max(np.abs(H_cortada_db - H_larga_db))
        diferencias_db.append(diferencia_db)
        print(f"L={l_test:3d}: diferencia máxima de |H| = {diferencia_db:.3f} dB")
    print("--------------------------------")
    return diferencias_db

if __name__ == "__main__":
    # Tests offset
    #   O1. entrada nula -> la salida es u[n % 16]
    test_O1()
    #   O2. x aleatoria, h Butterworth: r_con_offset - r_sin_offset == u[n % 16]
    test_O2()
    #   O3. random_offset: todos los valores en [-lsb, +lsb] y con la misma semilla sale el mismo u
    test_O3()

    # Tests BW
    #   B1. canal identidad (un 1 en la columna 0 de cada fila): salida == entrada
    test_B1()
    #   B2. random_butterworth con L = 200: la suma de los taps de cada fila da 1
    test_B2()
    #   B3. random_butterworth con L = 200: en cada fase |H(B_fases[p])| = -3,01 dB ± 0,05
    test_B3()
    #   B4. random_butterworth con L = 32: B_fases entre B·(1-d) y B·(1+d), filas distintas entre sí
    test_B4()
    #   B5. random_butterworth con desajuste = 0: planta == np.convolve(x, h[0])[:len(x)]
    test_B5()
    #   B6. random_butterworth(30e9, 4, 32, 0.25, semilla_bw): tiene que dar ValueError (corte > 36 GHz)
    test_B6()

    # Efecto de L (mediciones, sin assert), en el peor caso: corte en B·(1 + desajuste)
    #   L1. energía perdida para cada L, con umbral de -44 dB
    perdidas_db = test_L1()
    #   L2. diferencia máxima entre |H| del FIR cortado y |H| del largo (L = 400), hasta 30 GHz
    test_L2()

    print("\n Se pasaron todos los tests de offset y ancho de banda, y se midió el efecto de L.")

    # Gráfico de L1: energía perdida en el peor caso contra L
    plt.plot(L_BARRIDO, perdidas_db, "o-")
    plt.axhline(UMBRAL_DB, linestyle="--")
    plt.axvline(L, linestyle=":")
    plt.grid(True)
    plt.title(f"Energía perdida por truncamiento (corte en {B * (1 + desajuste) / 1e9:.1f} GHz, L = {L})")
    plt.xlabel("Cantidad de coeficientes por kernel")
    plt.ylabel("Energía perdida [dB]")
    plt.show()
