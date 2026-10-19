"""
cal_adc.py — Modelo comportamental del ADC de calibración (Cal ADC).

Cadena del modelo (en el mismo orden que en el hardware):

    r(t) -> [ pasabajos RC, ancho_banda ] -> T&H (disparo cada MN+1 muestras)
         -> [ × ganancia ] -> [ + offset ] -> [ recorte ±V_FS ] -> y_cal[m]

Los instantes de disparo siguen Galetto et al. (2021), ec. (14):
    t_m = m(MN)Ts + mTs  ->  un disparo cada MN + 1 muestras a FS.

Solo usa numpy.
"""
import numpy as np

# ===========================================================================
# Sistema
# ===========================================================================
FS = 72e9               # frecuencia de muestreo del sistema [Hz]
N  = 1.2                # muestras por símbolo (72 GHz / 60 GBd)
M  = 320                # submuestreo (320 o 640)
L  = 8                  # puntos de la grilla fina por cada muestra a 72 GHz


def filtro_rc(r, ancho_banda, L=8):
    """
    Pasabajos de 1er orden (polo RC) aplicado a la señal analógica en la grilla
    fina. Discretización bilineal con prewarping: el -3 dB cae exactamente en
    ancho_banda. Es causal, como el filtro real.

        y[n] = b·(x[n] + x[n-1]) - a·y[n-1]
    """
    fs_fina = L * FS
    if ancho_banda >= fs_fina / 2:
        raise ValueError("ANCHO_BANDA debe ser menor que la mitad de la grilla fina")
    K = np.tan(np.pi * ancho_banda / fs_fina)
    b = K / (1 + K)
    a = (K - 1) / (K + 1)

    y = np.empty_like(r, dtype=float)
    x_ant = 0.0                         # el filtro arranca descargado
    y_ant = 0.0
    for n in range(len(r)):
        y_ant = b * (r[n] + x_ant) - a * y_ant
        x_ant = r[n]
        y[n] = y_ant
    return y


def respuesta_rc(f, ancho_banda, L=8):
    """Respuesta en frecuencia H(f) del filtro_rc (para comparar con la teoría)."""
    fs_fina = L * FS
    K = np.tan(np.pi * ancho_banda / fs_fina)
    b = K / (1 + K)
    a = (K - 1) / (K + 1)
    z1 = np.exp(-2j * np.pi * f / fs_fina)          # z^-1
    return b * (1 + z1) / (1 + a * z1)


def adc_cal(r, M=320, L=8, ganancia=1.0, offset=0.0, ancho_banda=None, V_FS=None):
    """
    Entradas
        r           : señal del FI-DAC en la grilla fina de 576 GHz [V]
        M           : submuestreo (320 o 640)
        L           : puntos de la grilla fina por cada muestra a 72 GHz
        ganancia    : ganancia del ADC (1 = ideal)
        offset      : offset del ADC [V] (0 = ideal)
        ancho_banda : ancho de banda -3 dB del front-end [Hz] (None = ideal)
        V_FS        : rango de entrada ±V_FS [V] (None = sin recorte)
    Salidas
        y_cal : el valor entregado por el ADC en cada disparo [V]
        n_m   : índice a 72 GHz de cada disparo
    """
    # 1) Ancho de banda: filtra la señal analógica ANTES de muestrear
    if ancho_banda is not None:
        r = filtro_rc(r, ancho_banda, L)

    # 2) Disparos: uno cada MN + 1 muestras a 72 GHz (Galetto, ec. 14)
    paso = round(M * N) + 1
    i = np.arange(0, len(r), paso * L)
    n_m = i // L

    # 3) Track & hold: valor de la señal en el instante del disparo
    v = r[i]

    # 4) Ganancia y offset
    y_cal = ganancia * v + offset

    # 5) Rango de entrada finito
    if V_FS is not None:
        y_cal = np.clip(y_cal, -V_FS, V_FS)

    return y_cal, n_m
