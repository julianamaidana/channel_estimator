"""
ADC de calibración con jitter.

Dispara en los mismos instantes que adc_cal_ideal, pero cada disparo se corre un
tiempo delta ~ N(0, sigma^2), independiente entre disparos. Con Taylor de primer
orden:

    y_m = r(t_m) + delta_m * r'(t_m)

La derivada se aproxima con la diferencia hacia adelante en la grilla fina.
El front-end (ganancia, offset, ancho de banda) no está incluido; si se agrega,
va antes de adc_cal.
"""
import numpy as np

from adc_cal_ideal import adc_cal_ideal, FS


def adc_cal(x, sigma_jitter=0.0, M=320, L=8, semilla=None):
    """
    x: señal a la entrada del ADC, en la grilla fina [V]
    sigma_jitter: jitter rms [s]
    semilla: para repetir la misma realización del jitter

    Devuelve y_cal y n_m, igual que adc_cal_ideal. Con sigma_jitter = 0 da lo
    mismo que el ideal.
    """
    y_ideal, n_m = adc_cal_ideal(x, M=M, L=L)
    i = n_m * L

    # la pendiente usa el punto siguiente: si el último disparo cae en el último
    # punto de x, se descarta
    ok = i + 1 < len(x)
    y_ideal, n_m, i = y_ideal[ok], n_m[ok], i[ok]

    rng = np.random.default_rng(semilla)
    delta = sigma_jitter * rng.standard_normal(len(i))   # [s]

    dt = 1 / (L * FS)                                    # 1.74 ps con L = 8
    pendiente = (x[i + 1] - x[i]) / dt                   # [V/s]

    y_cal = y_ideal + delta * pendiente
    return y_cal, n_m