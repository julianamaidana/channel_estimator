"""
SNR por jitter del ADC cal en función de sigma: simulación vs. teoría.

Para cada sigma se muestrea la misma señal con y sin jitter; el error es la
diferencia. Teoría para un tono de frecuencia f (Kester, MT-007):

    SNR = -20*log10(2*pi*f*sigma)

Para la señal de prueba se usa f_rms en lugar de f.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from modules.cal_adc.python.adc_cal_ideal import FS
from modules.cal_adc.python.adc_cal_efecto_jitter import adc_cal

L = 8
DISPAROS = 2000
LARGO = DISPAROS * 385 * L      # 385 = paso con M = 320
SIGMAS = np.array([25, 50, 100, 200, 400]) * 1e-15

QPSK = (-1.0, 1.0)              # un riel (I o Q); cada riel va por su propio FI-DAC
PAM4 = (-3.0, -1.0, 1.0, 3.0)


def snr_jitter(r, sigma):
    """SNR [dB] debido solo al jitter."""
    ideal, _ = adc_cal(r, sigma_jitter=0, L=L)
    con_jitter, _ = adc_cal(r, sigma_jitter=sigma, L=L, semilla=1)
    error = con_jitter - ideal
    return 10 * np.log10(np.mean(ideal**2) / np.mean(error**2))


def f_rms(r):
    """f_rms = rms(r') / (2*pi*rms(r)). Para un tono da su frecuencia."""
    pendiente = np.diff(r) * (L * FS)
    return np.sqrt(np.mean(pendiente**2)) / (2 * np.pi * np.sqrt(np.mean(r**2)))


def senal_prueba(largo, niveles=QPSK, semilla=0):
    """
    Símbolos al azar a 60 GBd interpolados a la grilla fina (banda limitada a
    30 GHz). Se usa en lugar de la salida del FI-DAC hasta tener ese modelo.
    """
    rng = np.random.default_rng(semilla)
    simbolos = rng.choice(niveles, largo * 5 // (6 * L) + 40)
    return signal.resample_poly(simbolos, 6 * L, 5)[:largo]   # 9.6 puntos por símbolo


if __name__ == "__main__":
    t = np.arange(LARGO) / (L * FS)
    senales = {
        "tono de 10 GHz": np.sin(2 * np.pi * 10e9 * t),
        "riel de QPSK de prueba (60 GBd)": senal_prueba(LARGO),
        "tono de 30 GHz (borde de banda)": np.sin(2 * np.pi * 30e9 * t),
    }
    del t

    snr_medido, teoria = {}, {}
    for nombre, r in senales.items():
        fe = f_rms(r)
        snr_medido[nombre] = np.array([snr_jitter(r, s) for s in SIGMAS])
        teoria[nombre] = -20 * np.log10(2 * np.pi * fe * SIGMAS)
        print(f"\n{nombre}   (f_rms = {fe / 1e9:.1f} GHz)")
        print("   sigma    SNR modelo   SNR teoría   señal/error")
        for s, snr, te in zip(SIGMAS, snr_medido[nombre], teoria[nombre]):
            print(f"  {s * 1e15:4.0f} fs   {snr:6.1f} dB    {te:6.1f} dB    {10 ** (snr / 20):11.0f}")

    colores = ["#2a78d6", "#eb6834", "#1baf7a"]
    x = SIGMAS * 1e15
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for nombre, color in zip(senales, colores):
        ax.plot(x, teoria[nombre], color=color, lw=2, label=nombre)
        ax.plot(x, snr_medido[nombre], "o", color=color, ms=8, mfc="white", mew=2)
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:.0f}" for v in x])
    ax.set_xlabel("jitter rms σ [fs]")
    ax.set_ylabel("SNR por jitter [dB]")
    ax.set_title("SNR por jitter en el ADC cal: puntos = modelo, líneas = teoría", fontsize=10, loc="left")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    plt.show()
