"""
sim_cal_adc.py — Compara la señal "analógica" con la salida del Cal ADC.

Aplica a un tono de 10 GHz el modelo del Cal ADC con los efectos que se
configuran abajo (GANANCIA, OFFSET, ANCHO_BANDA, V_FS), verifica la amplitud,
fase y offset estimados contra la teoría y grafica:
    Panel 1: señal analógica + señal vista por el T&H + muestras del ADC.
    Panel 2: salida del ADC vs. valor exacto en cada disparo, y el error.

Uso (desde la raíz del repositorio):
    python modules/cal_adc/simulations/sim_cal_adc.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Modelo del Cal ADC (modules/cal_adc/python/cal_adc.py)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from cal_adc import FS, N, M, L, adc_cal, filtro_rc, respuesta_rc  # noqa: E402

# ===========================================================================
# Efectos del ADC  (valores ideales: GANANCIA = 1, OFFSET = 0, ANCHO_BANDA = None)
# ===========================================================================
GANANCIA    = 1.05      # ganancia del ADC (1.05 = +5 % de error de ganancia)
OFFSET      = 0.05      # offset del ADC [V]
ANCHO_BANDA = 15e9      # ancho de banda -3 dB del front-end (T&H) [Hz]; None = infinito
V_FS        = None      # rango de entrada ±V_FS [V]; None = sin recorte

PUNTOS = 200000         # cantidad de puntos de la grilla fina (≈ 347 ns a 576 GHz)

if __name__ == "__main__":
    f = 10e9                                   # tono de 10 GHz
    A = 1.0                                    # amplitud [V]

    t = np.arange(PUNTOS) / (L * FS)
    r = A * np.sin(2 * np.pi * f * t)

    y_cal, n_m = adc_cal(r, M=M, L=L, ganancia=GANANCIA, offset=OFFSET,
                         ancho_banda=ANCHO_BANDA, V_FS=V_FS)

    t_m = n_m / FS
    exacto = A * np.sin(2 * np.pi * f * t_m)   # lo que entregaría un ADC ideal
    err = y_cal - exacto
    fs_adc = FS / (round(M * N) + 1)

    # Respuesta esperada del front-end a la frecuencia del tono
    if ANCHO_BANDA is not None:
        H = respuesta_rc(f, ANCHO_BANDA, L)
        H_analog = 1 / (1 + 1j * f / ANCHO_BANDA)
        r_bw = filtro_rc(r, ANCHO_BANDA, L)    # solo para el gráfico
        tau = 1 / (2 * np.pi * ANCHO_BANDA)
    else:
        H = H_analog = 1.0 + 0j
        r_bw = r
        tau = 0.0

    # Verificación: ajuste y = a·sin + b·cos + c (frecuencia conocida),
    # descartando el transitorio inicial del filtro (primeros 10·tau)
    ok = t_m > 10 * tau
    X = np.column_stack([np.sin(2 * np.pi * f * t_m[ok]),
                         np.cos(2 * np.pi * f * t_m[ok]),
                         np.ones(ok.sum())])
    (c_s, c_c, os_est), *_ = np.linalg.lstsq(X, y_cal[ok], rcond=None)
    amp_est = np.hypot(c_s, c_c) / A
    fase_est = np.arctan2(c_c, c_s)

    print(f"El ADC se quedó con {len(y_cal)} muestras (fs_ADC = {fs_adc/1e6:.2f} MHz)")
    print(f"{'':22s}{'esperado':>12s}{'estimado':>12s}")
    print(f"{'amplitud G·|H|':22s}{GANANCIA*abs(H):12.4f}{amp_est:12.4f}")
    print(f"{'fase ∠H [rad]':22s}{np.angle(H):+12.4f}{fase_est:+12.4f}")
    print(f"{'offset [mV]':22s}{OFFSET*1e3:12.2f}{os_est*1e3:12.2f}")
    if ANCHO_BANDA is not None:
        print(f"Polo RC analógico ideal: |H| = {abs(H_analog):.4f}, ∠H = {np.angle(H_analog):+.4f} rad")
    print(f"Error máximo vs. ADC ideal: {np.max(np.abs(err[ok])):.4f} V  |  MSE: {np.mean(err[ok]**2):.2e} V²")
    if V_FS is not None:
        print(f"Muestras recortadas en +V_FS: {np.sum(y_cal >= V_FS)}  |  en −V_FS: {np.sum(y_cal <= -V_FS)}")

    # ---------------- Gráficos ----------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7),
                                   gridspec_kw={"height_ratios": [3, 2]})

    ax1.plot(t * 1e9, r, lw=1, color="C0",
             label=f"señal analógica r (seno, {f/1e9:g} GHz)")
    if ANCHO_BANDA is not None:
        ax1.plot(t * 1e9, r_bw, lw=1, color="C1",
                 label=f"señal vista por el T&H (BW = {ANCHO_BANDA/1e9:g} GHz)")
    ax1.step(np.append(t_m, t[-1]) * 1e9, np.append(y_cal, y_cal[-1]), where="post",
             lw=1.5, color="C3", alpha=0.6, label="salida ADC (hold)")
    ax1.plot(t_m * 1e9, y_cal, "o", color="C3",
             label=f"muestras ADC (G = {GANANCIA:g}, offset = {OFFSET*1e3:g} mV, fs = {fs_adc/1e6:.2f} MHz)")
    ax1.set_xlabel("t [ns]"); ax1.set_ylabel("amplitud [V]")
    ax1.set_title("Señal analógica vs. salida del ADC de calibración")
    ax1.grid(alpha=0.3); ax1.legend(loc="upper right", fontsize=8)

    m = np.arange(len(y_cal))
    ax2.stem(m, exacto, linefmt="C0-", markerfmt="C0s", basefmt=" ",
             label="valor exacto r(t_m)")
    ax2.plot(m, y_cal, "x", color="C3", ms=9, mew=2, label="salida ADC y_cal[m]")
    ax2.plot(m, err, "^", color="C2", label="error y_cal[m] − r(t_m)")
    ax2.set_xlabel("índice de disparo del ADC m"); ax2.set_ylabel("amplitud [V]")
    ax2.grid(alpha=0.3); ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    plt.show()
