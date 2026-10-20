"""
Figuras para la presentación: muestreo del ADC cal y efecto del jitter.

    python presentacion_jitter.py         # las 4 figuras
    python presentacion_jitter.py 3 4     # solo la 3 y la 4

Se guardan en docs/img/cal_adc/ del repo. La señal de prueba es un riel de QPSK
(senal_prueba, de estudio_jitter.py).
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(RAIZ))
from modules.cal_adc.python.adc_cal_ideal import adc_cal_ideal, FS, N
from modules.cal_adc.python.adc_cal_efecto_jitter import adc_cal
from modules.cal_adc.simulations.sim_estudio_jitter import senal_prueba, snr_jitter, f_rms

MOSTRAR = True              # False: guarda las figuras sin abrir ventanas
CARPETA = Path(__file__).resolve().parents[3] / "docs" / "img" / "cal_adc"
L = 8
M = 320
PASO = round(M * N) + 1     # 385
T_DISPARO = PASO / FS       # 5.35 ns
DISPAROS = 2000
P = 16                      # período del canal, en muestras

AZUL, NARANJA, VERDE, GRIS = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.titlelocation": "left",
                     "axes.grid": True, "grid.alpha": 0.3, "axes.axisbelow": True,
                     "axes.spines.top": False, "axes.spines.right": False})


def tiempo_ns(i):
    return i / (L * FS) * 1e9


def pendiente_en(r, i):
    # misma aproximación que en adc_cal [V/s]
    return (r[i + 1] - r[i]) * (L * FS)


def figura_1(r):
    """Qué hace el ADC cal: una muestra cada 385."""
    y_cal, n_m = adc_cal_ideal(r, M=M, L=L)

    fig = plt.figure(figsize=(13, 7))
    grilla = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])

    # arriba: los primeros 4 disparos
    ax = fig.add_subplot(grilla[0, :])
    fin = 3 * PASO * L + 60 * L
    ax.plot(tiempo_ns(np.arange(fin)), r[:fin], color=GRIS, lw=0.5,
            label="señal de prueba: un riel de QPSK (grilla fina, 576 GHz)")
    ax.plot(n_m[:4] / FS * 1e9, y_cal[:4], "o", color=NARANJA, ms=11, label="lo que toma el ADC cal")
    ax.set_title(f"El ADC cal toma un punto cada {PASO} muestras = {T_DISPARO * 1e9:.2f} ns "
                 f"(ADC a {FS / PASO / 1e6:.0f} MHz)")
    ax.set_xlabel("tiempo [ns]")
    ax.set_ylabel("tensión [V]")
    lo, hi = r[:fin].min(), r[:fin].max()
    ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.45 * (hi - lo))   # lugar para la leyenda
    ax.legend(loc="upper right", fontsize=9, ncol=2)

    # abajo: zoom en los disparos 1, 2 y 3
    for col, m in enumerate([1, 2, 3]):
        ax = fig.add_subplot(grilla[1, col])
        i0 = n_m[m] * L
        tramo = np.arange(i0 - 60, i0 + 61)   # +-104 ps
        ax.plot(tiempo_ns(tramo), r[tramo], ".-", color=AZUL, ms=3, lw=1.2, label="señal (puntos de la grilla fina)")
        ax.plot(n_m[m] / FS * 1e9, y_cal[m], "o", color=NARANJA, ms=11, label="muestra del ADC")
        ax.set_title(f"disparo {m}: n_m = {n_m[m]}, fase {n_m[m] % P}", fontsize=10)
        ax.set_xlabel("tiempo [ns]")
        if col == 0:
            ax.set_ylabel("tensión [V]")
            ax.legend(loc="best", fontsize=8)
    fig.tight_layout()

    print(f"paso: {PASO} muestras = {T_DISPARO * 1e9:.2f} ns ({FS / PASO / 1e6:.0f} MHz), "
          f"{PASO / N:.0f} símbolos entre disparos")
    return fig


def figura_2(r):
    """El +1: cada disparo avanza una fase del canal."""
    _, n_m = adc_cal_ideal(r, M=M, L=L)
    fase = n_m % P
    fase_sin = (np.arange(len(n_m)) * round(M * N)) % P   # mismo ADC con paso 384

    fig, ejes = plt.subplots(1, 2, figsize=(13, 4.8))
    m = np.arange(48)
    ejes[0].plot(m, fase[:48], "o-", color=AZUL, ms=6, lw=1, label=f"con el +1 (paso {PASO})")
    ejes[0].plot(m, fase_sin[:48], "s", color=NARANJA, ms=6, label=f"sin el +1 (paso {PASO - 1})")
    ejes[0].set_title("Cada disparo avanza una fase del canal")
    ejes[0].set_xlabel("número de disparo")
    ejes[0].set_ylabel(f"fase del canal (n_m % {P})")
    ejes[0].set_yticks(range(0, P, 3))
    ejes[0].set_ylim(-1, P + 3)
    ejes[0].legend(loc="upper left", ncol=2, fontsize=9)

    ejes[1].bar(range(P), np.bincount(fase, minlength=P), color=AZUL)
    ejes[1].set_title(f"Las {P} fases reciben la misma cantidad de disparos")
    ejes[1].set_xlabel(f"fase del canal (n_m % {P})")
    ejes[1].set_ylabel("cantidad de disparos")
    ejes[1].set_xticks(range(P))
    fig.tight_layout()

    print(f"{PASO} = {PASO // P}*{P} + {PASO % P}")
    print(f"fases visitadas: {len(np.unique(fase))} con el +1, {len(np.unique(fase_sin))} sin el +1")
    return fig


def figura_3(r):
    """Error de jitter = delta * pendiente, en un punto plano y en uno empinado."""
    _, n_m = adc_cal_ideal(r, M=M, L=L)
    i = n_m[1:400] * L
    s = pendiente_en(r, i)
    elegidos = [i[np.argmin(np.abs(s))], i[np.argmax(np.abs(s))]]
    titulos = ["Señal casi plana: el error es chico", "Señal empinada: el mismo δ da un error grande"]
    sigma_vis = 2e-12    # exagerado para que se vea (la referencia es 100 fs)
    rng = np.random.default_rng(3)

    fig, ejes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, i0, titulo in zip(ejes, elegidos, titulos):
        tramo = np.arange(i0 - 40, i0 + 41)   # +-69 ps
        t_ps = (tramo - i0) / (L * FS) * 1e12
        pend = pendiente_en(r, i0)
        delta = sigma_vis * rng.standard_normal(80)
        ax.plot(t_ps, r[tramo], color=AZUL, lw=2, label="señal")
        tt = np.array([-4, 4]) * sigma_vis
        ax.plot(tt * 1e12, r[i0] + pend * tt, "--", color=GRIS, lw=1.2, label="pendiente (recta tangente)")
        ax.plot(delta * 1e12, r[i0] + pend * delta, "o", color=NARANJA, ms=5, alpha=0.6,
                label="lo que mide con jitter: r(t) + δ·pendiente")
        ax.plot(0, r[i0], "o", ms=11, mfc="white", mec=AZUL, mew=2.5, label="lo que debería medir")
        dv = sigma_vis * abs(pend)
        ax.annotate("", xy=(12, r[i0] + dv), xytext=(12, r[i0] - dv),
                    arrowprops=dict(arrowstyle="|-|,widthA=0.4,widthB=0.4", color=NARANJA, lw=1.5))
        ax.text(0.98, 0.96, f"barra naranja: ±σ·|pendiente| = ±{dv * 1e3:.1f} mV", transform=ax.transAxes,
                ha="right", va="top", fontsize=9, color=NARANJA)
        ax.set_title(titulo)
        ax.set_xlabel("tiempo desde el disparo [ps]")
    ejes[0].set_ylabel("tensión [V]")
    ejes[0].legend(loc="best", fontsize=8)
    fig.suptitle("Jitter: el disparo se corre δ en el tiempo → error de tensión = δ × pendiente   "
                 "(σ exagerado a 2 ps para que se vea; el valor de referencia es 100 fs)", fontsize=11, x=0.01, ha="left")
    fig.tight_layout()

    ideal, _ = adc_cal(r, 0, M=M, L=L)
    con, _ = adc_cal(r, 100e-15, M=M, L=L, semilla=1)
    err = con - ideal
    print(f"sigma = 100 fs: error rms = {np.sqrt(np.mean(err ** 2)) * 1e3:.1f} mV, "
          f"señal rms = {np.sqrt(np.mean(ideal ** 2)):.2f} V")
    return fig


def figura_4(r):
    """SNR por jitter en función de sigma: modelo (puntos) y teoría (líneas)."""
    sigmas = np.array([25, 50, 100, 200, 400]) * 1e-15
    t = np.arange(len(r)) / (L * FS)
    senales = {"tono de 10 GHz": np.sin(2 * np.pi * 10e9 * t),
               "riel de QPSK de prueba (60 GBd)": r,
               "tono de 30 GHz (borde de banda)": np.sin(2 * np.pi * 30e9 * t)}

    fig, ax = plt.subplots(figsize=(10, 5.5))
    con_100 = {}
    print("SNR con sigma = 100 fs:")
    for (nombre, x), color in zip(senales.items(), [AZUL, NARANJA, VERDE]):
        fe = f_rms(x)
        medido = [snr_jitter(x, s) for s in sigmas]
        teoria = -20 * np.log10(2 * np.pi * fe * sigmas)
        con_100[nombre] = medido[2]   # sigmas[2] = 100 fs
        ax.plot(sigmas * 1e15, teoria, color=color, lw=2, label=f"{nombre}  (f_rms = {fe / 1e9:.1f} GHz)")
        ax.plot(sigmas * 1e15, medido, "o", color=color, ms=9, mfc="white", mew=2)
        print(f"  {nombre:32s} modelo {medido[2]:5.1f} dB, teoría {teoria[2]:5.1f} dB, "
              f"señal/error = {10 ** (medido[2] / 20):.0f}")
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks(sigmas * 1e15)
    ax.set_xticklabels([f"{v:.0f}" for v in sigmas * 1e15])
    ax.set_xlabel("jitter rms σ [fs]")
    ax.set_ylabel("SNR por jitter [dB]")
    ax.set_title("Cuánto afecta el jitter: SNR = −20·log10(2π·f·σ)   (puntos = modelo, líneas = teoría)")
    ax.legend(fontsize=9)

    # referencia de 100 fs
    ax.axvline(100, color=GRIS, lw=1, ls="--", zorder=0)
    q = con_100["riel de QPSK de prueba (60 GBd)"]
    ax.text(0.02, 0.04,
            "Con σ = 100 fs (línea punteada):\n"
            f"   tono de 10 GHz: {con_100['tono de 10 GHz']:.1f} dB\n"
            f"   riel de QPSK: {q:.1f} dB → la señal es unas {10 ** (q / 20):.0f} veces el error\n"
            f"   tono de 30 GHz: {con_100['tono de 30 GHz (borde de banda)']:.1f} dB",
            transform=ax.transAxes, fontsize=9, va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GRIS, lw=0.8))
    fig.tight_layout()
    return fig


# número de figura: (nombre del archivo, función que la arma)
FIGURAS = {
    1: ("submuestreo", figura_1),
    2: ("recorrido_fases", figura_2),
    3: ("error_jitter", figura_3),
    4: ("snr_vs_jitter", figura_4),
}

if __name__ == "__main__":
    CARPETA.mkdir(parents=True, exist_ok=True)
    pedidas = [int(a) for a in sys.argv[1:]] or list(FIGURAS)
    r = senal_prueba(DISPAROS * PASO * L)
    for k in pedidas:
        nombre, armar = FIGURAS[k]
        print(f"\n--- figura {k}: {nombre}.png ---")
        fig = armar(r)
        fig.savefig(CARPETA / f"{nombre}.png", dpi=150)
        if MOSTRAR:
            plt.show()
        plt.close(fig)