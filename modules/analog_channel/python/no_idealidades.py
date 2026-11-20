"""
no_idealidades.py — No idealidades del camino analógico: ancho de banda y offset (issue #4).

Ancho de banda
    Cada fase p tiene su propio filtro Butterworth, con un ancho de banda B_p
    sorteado alrededor del nominal B. El kernel de la fase es su respuesta al
    impulso, truncada a L coeficientes:

        h[p, :] = butterworth(B_p)    truncado a L coeficientes

Largo L
    Lo fija la fase que corta más alto, porque es la de respuesta más larga.
    largo_minimo() busca el L más chico que, para ese corte, pierde menos
    energía que un umbral (-44 dB: el ruido de un DAC con ENOB de 7 bits).
    Con B = 30 GHz y ±15 %, el peor corte es 34,5 GHz y da L = 74.

Offset
    No entra en h, porque con entrada nula los kernels dan salida nula. Va
    aparte: un valor fijo por fase, sumado a la salida de la planta:

        r[n] = sum_k h[p, k]·x[n-k] + u[p],    con p = n mod 16

Ganancia y skew
    Los modela Misael, también como kernels de 16 x L. Para combinar
    los efectos, se convoluciona fila por fila y se vuelve a truncar a L:

        h_total[p, :] = h_bw[p, :] * h_otro[p, :]

Referencias: Bonetti 2026 (Sec. II, Ec. 17, y Sec. V-A), Galetto 2021 (Sec. V).
"""

import numpy as np
from fir_periodico import planta_fir_periodica, N_FASES

from scipy import signal

FS = 72e9

# ===========================================================================
# Offset por fase
# ===========================================================================

def planta_con_offset(x, h, u):
    """
    Planta de 16 fases con offset:
        r[n] = planta_fir_periodica(x, h)[n] + u[n mod 16]

    Entradas
        x : señal de entrada
        h : matriz de kernels, 16 x L
        u : offset de cada fase, 16 valores
    Salida
        r : señal de salida, mismo largo que x
    """
    r = planta_fir_periodica(x, h)
    for n in range(len(r)):
        p = n % N_FASES
        r[n] += u[p]
    return r

def random_offset(lsb, semilla):
    """
    Sortea un offset aleatorio para cada fase, entre -lsb y +lsb

    Entradas
        lsb : tamaño del LSB (en la misma unidad que el offset)
        semilla : semilla para el generador de números aleatorios
    Salida
        u : vector de 16 offsets, uno por fase
    """
    rng = np.random.default_rng(semilla)
    u = rng.uniform(-lsb, lsb, N_FASES)
    return u

# ===========================================================================
# BW
# ===========================================================================

def kernel_butterworth(B, orden, L):
    """
    Kernel de Butterworth (truncado a L coeficientes).

    Entradas
        B : ancho de banda (Hz)
        orden : orden del filtro
        L : cantidad de coeficientes del kernel
    Salida
        h : kernel de Butterworth, largo L
    """
    b, a = signal.butter(orden, B, fs=FS)
    # Obtener la respuesta al impulso del filtro digital
    impulso = np.zeros(L)
    impulso[0] = 1
    h = signal.lfilter(b, a, impulso)
    return h


def random_butterworth(B, orden, L, desajuste, semilla):
    """
    Arma un kernel de Butterworth distinto para cada fase. El ancho de banda
    de cada fase se sortea alrededor del nominal B, entre B·(1 - desajuste)
    y B·(1 + desajuste).

    Entradas
        B         : ancho de banda nominal [Hz]
        orden     : orden del filtro
        L         : cantidad de coeficientes por kernel
        desajuste : desajuste máximo, relativo (0.15 = ±15 %)
        semilla   : semilla del generador de números aleatorios
    Salidas
        h       : matriz de kernels, 16 x L
        B_fases : ancho de banda sorteado para cada fase [Hz], 16 valores
    """
    if B * (1 + desajuste) >= FS / 2:
        raise ValueError("B * (1 + desajuste) tiene que ser menor que FS/2")

    rng = np.random.default_rng(semilla)
    B_fases = B * (1 + rng.uniform(-desajuste, desajuste, N_FASES))

    h = np.zeros((N_FASES, L))
    for p in range(N_FASES):
        h[p, :] = kernel_butterworth(B_fases[p], orden, L)

    return h, B_fases


def largo_minimo(B_corte, orden, umbral_db):
    """
    L más chico tal que, al truncar el kernel de corte B_corte a L coeficientes,
    la energía que se tira (la cola) queda por debajo de umbral_db respecto
    de la energía total.

    Para elegir el L del modelo se usa el peor caso, B_corte = B·(1 + desajuste):
    la fase que corta más alto es la de respuesta más larga.

    Entradas
        B_corte   : ancho de banda del kernel [Hz]
        orden     : orden del filtro
        umbral_db : energía perdida máxima [dB], por ejemplo -44
    Salida
        L : cantidad de coeficientes
    """
    h_larga = kernel_butterworth(B_corte, orden, 2000)   # largo de sobra: la cola ya es despreciable
    E_total = np.sum(h_larga ** 2)

    for L in range(1, 2000):
        E_perdida = np.sum(h_larga[L:] ** 2)              # energía de la cola que se tira
        if 10 * np.log10(E_perdida / E_total) < umbral_db:
            return L

    raise ValueError("Ningún L menor que 2000 cumple el umbral")
