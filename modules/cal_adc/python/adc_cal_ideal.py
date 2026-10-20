"""
ADC de calibración ideal (sin jitter ni otros errores).

Dispara cada M*N + 1 muestras a 72 GHz (Galetto 2021, ec. 14): con M = 320 son
385 muestras, o sea 5.35 ns. La señal de entrada viene en una grilla L veces más
fina (576 GHz con L = 8), así que en esa grilla el paso es 385*8 = 3080 puntos.
"""
import numpy as np
import matplotlib.pyplot as plt

FS = 72e9   # frecuencia de muestreo del sistema [Hz]
N = 1.2     # muestras por símbolo (72 GHz / 60 GBd)
L = 8       # sobremuestreo de la grilla fina


def adc_cal_ideal(r, M=320, L=8):
    """
    r: señal del FI-DAC en la grilla fina [V]
    M: factor de submuestreo (320 o 640)
    L: sobremuestreo con el que viene r

    Devuelve y_cal (valor de la señal en cada disparo) y n_m (índice a 72 GHz de
    cada disparo, para saber contra qué muestra de la referencia comparar).
    """
    paso = round(M * N) + 1     # 385; el +1 hace que recorra todas las fases del canal
    i = np.arange(0, len(r), paso * L)
    y_cal = r[i]
    n_m = i // L
    return y_cal, n_m