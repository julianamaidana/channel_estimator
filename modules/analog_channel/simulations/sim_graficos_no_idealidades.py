"""
sim_graficos_no_idealidades.py — Gráficos del modelo de ancho de banda y offset (issue #4).

Usa las funciones de no_idealidades.py con los mismos parámetros y semillas que
sim_no_idealidades_bw_offs.py, así las figuras coinciden con lo validado.
Cada figura se guarda en PNG (presentación) y PDF (reporte LaTeX) en la carpeta
figuras/, al lado de este script. Esa carpeta no se sube al repositorio.

Uso:
    python sim_graficos_no_idealidades.py

Figuras
    F1. Magnitud de los 16 kernels: el desajuste de ancho de banda entre fases.
    F2. Fase de los 16 kernels y su diferencia respecto del nominal.
    F3. Taps de la fase más lenta y la más rápida: dónde se corta la cola.
    F4. Efecto de L: energía perdida y error de |H| en la peor fase.
    F5. Offset: salida con entrada nula, un patrón que se repite cada 16 muestras.
    F6. Espectro con un tono a la entrada: espurias en ±f0 + k·4,5 GHz.
    F7. Misma señal por el canal nominal y por el desajustado, y el error.
    F8. SNR del canal respecto del nominal contra el desajuste.
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

from fir_periodico import planta_fir_periodica, N_FASES  # type: ignore
from no_idealidades import planta_con_offset, kernel_butterworth, random_offset, random_butterworth, largo_minimo  # type: ignore

# ---------------------------------------------------------------------------
# Parámetros (los mismos que en sim_no_idealidades_bw_offs.py)
# ---------------------------------------------------------------------------
FS = 72e9
BITS = 8
BITS_ENOB = 7
RANGO = 1.0
B = 30e9
orden = 4
desajuste = 0.15
semilla_o = 0
semilla_x = 1
semilla_bw = 2

LSB = RANGO / (2 ** BITS)
UMBRAL_DB = -(6.02 * BITS_ENOB + 1.76)   # ≈ -44 dB

# El L lo fija la fase que corta más alto (la más cerca de Nyquist). Con un desajuste de ±d,
# el peor caso posible es una fase que corte justo en el borde, B·(1 + d): no depende del sorteo.
B_PEOR = B * (1 + desajuste)
L_ELEGIDO = largo_minimo(B_PEOR, orden, UMBRAL_DB)   # peor caso: se recalcula si cambiás B o el desajuste

MOSTRAR = True                            # False: solo guarda las figuras, no abre ventanas
CARPETA = Path(__file__).resolve().parent / "figuras"   # no se sube al repositorio

# ---------------------------------------------------------------------------
# Estilo: fondo claro, grilla tenue, líneas finas
# ---------------------------------------------------------------------------
AZUL = "#2a78d6"        # serie 1
NARANJA = "#eb6834"     # serie 2
GRIS = "#52514e"        # nominal, umbrales y referencias
MAPA_B = plt.cm.Blues   # color de cada fase según su ancho de banda (más oscuro = más ancho)

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb",
    "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#b8b7b1",
    "axes.grid": True,
    "grid.color": "#e4e3df",
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.frameon": False,
    "lines.linewidth": 1.6,
    "font.size": 10,
})


def guardar(fig, nombre):
    """Guarda la figura en PNG y PDF dentro de figuras/."""
    CARPETA.mkdir(exist_ok=True)
    fig.savefig(CARPETA / f"{nombre}.png", dpi=200, bbox_inches="tight")
    fig.savefig(CARPETA / f"{nombre}.pdf", bbox_inches="tight")
    print(f"Guardada: figuras/{nombre}.png y .pdf")


def color_de_fase(B_fases, p):
    """Color de la fase p según su ancho de banda: del celeste (más angosto) al azul oscuro (más ancho)."""
    t = (B_fases[p] - B_fases.min()) / (B_fases.max() - B_fases.min())
    return MAPA_B(0.35 + 0.6 * t)


def barra_de_color(fig, ax, B_fases):
    """Barra de color que traduce el tono de cada curva a su ancho de banda."""
    norma = plt.Normalize(B_fases.min() / 1e9, B_fases.max() / 1e9)
    mapa = plt.cm.ScalarMappable(norm=norma, cmap=plt.cm.colors.ListedColormap(MAPA_B(np.linspace(0.35, 0.95, 256))))
    barra = fig.colorbar(mapa, ax=ax, pad=0.02)
    barra.set_label("Ancho de banda de la fase [GHz]")


# ===========================================================================
# F1. Magnitud de los 16 kernels
# ===========================================================================
def grafico_magnitud():
    h, B_fases = random_butterworth(B, orden, L_ELEGIDO, desajuste, semilla_bw)
    h_nominal = kernel_butterworth(B, orden, L_ELEGIDO)
    f = np.linspace(0, FS / 2, 2000)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax in (ax1, ax2):
        for p in range(N_FASES):
            w, H = signal.freqz(h[p, :], 1, worN=f, fs=FS)
            ax.plot(f / 1e9, 20 * np.log10(np.abs(H)), color=color_de_fase(B_fases, p), linewidth=1.2)
        w, H = signal.freqz(h_nominal, 1, worN=f, fs=FS)
        ax.plot(f / 1e9, 20 * np.log10(np.abs(H)), color=GRIS, linestyle="--", linewidth=1.8, label=f"Nominal ({B/1e9:.0f} GHz)")
        ax.axhline(-3, color=GRIS, linewidth=0.8, linestyle=":")
        ax.set_xlabel("Frecuencia [GHz]")

    ax1.set_xlim(0, 36); ax1.set_ylim(-45, 3)
    ax1.set_ylabel("|H| [dB]")
    ax1.set_title("Las 16 fases, de 0 a 36 GHz")
    ax1.legend(loc="lower left")
    ax1.text(1, -2.2, "−3 dB", color=GRIS)

    ax2.set_xlim(20, 36); ax2.set_ylim(-12, 1)
    ax2.set_title("Zoom en la zona de corte")
    barra_de_color(fig, ax2, B_fases)

    fig.suptitle(f"Respuesta en magnitud de los kernels (B0 = {B/1e9:.0f} GHz, desajuste ±{desajuste*100:.0f} %, L = {L_ELEGIDO})", y=1.02)
    guardar(fig, "F1_magnitud_16_fases")


# ===========================================================================
# F2. Fase de los 16 kernels
# ===========================================================================
def grafico_fase():
    h, B_fases = random_butterworth(B, orden, L_ELEGIDO, desajuste, semilla_bw)
    h_nominal = kernel_butterworth(B, orden, L_ELEGIDO)
    f = np.linspace(0, B, 1500)                       # hasta 30 GHz: la banda de la señal

    w, H_nom = signal.freqz(h_nominal, 1, worN=f, fs=FS)
    fase_nom = np.degrees(np.unwrap(np.angle(H_nom)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    for p in range(N_FASES):
        w, H = signal.freqz(h[p, :], 1, worN=f, fs=FS)
        fase = np.degrees(np.unwrap(np.angle(H)))
        ax1.plot(f / 1e9, fase, color=color_de_fase(B_fases, p), linewidth=1.2)
        ax2.plot(f / 1e9, fase - fase_nom, color=color_de_fase(B_fases, p), linewidth=1.2)
    ax1.plot(f / 1e9, fase_nom, color=GRIS, linestyle="--", linewidth=1.8, label="Nominal")

    ax1.set_title("Fase de las 16 fases")
    ax1.set_xlabel("Frecuencia [GHz]"); ax1.set_ylabel("Fase [°]")
    ax1.legend(loc="lower left")
    ax2.axhline(0, color=GRIS, linewidth=0.8)
    ax2.set_title("Diferencia respecto del nominal")
    ax2.set_xlabel("Frecuencia [GHz]"); ax2.set_ylabel("Fase − fase nominal [°]")
    barra_de_color(fig, ax2, B_fases)

    fig.suptitle("Respuesta en fase de los kernels, hasta 30 GHz", y=1.02)
    guardar(fig, "F2_fase_16_fases")


# ===========================================================================
# F3. Taps de la fase más lenta y la más rápida
# ===========================================================================
def grafico_taps():
    h_larga, B_fases = random_butterworth(B, orden, 200, desajuste, semilla_bw)
    p_lenta = np.argmin(B_fases)                      # la que corta más bajo
    p_rapida = np.argmax(B_fases)                     # la que corta más alto (más cerca de Nyquist)
    k = np.arange(200)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))

    # Izquierda: los primeros taps, en lineal
    n_ver = 25
    ax1.plot(k[:n_ver], h_larga[p_lenta, :n_ver], "o-", color=AZUL, markersize=5,
             label=f"Fase {p_lenta} ({B_fases[p_lenta]/1e9:.1f} GHz)")
    ax1.plot(k[:n_ver], h_larga[p_rapida, :n_ver], "o-", color=NARANJA, markersize=5,
             label=f"Fase {p_rapida} ({B_fases[p_rapida]/1e9:.1f} GHz)")
    ax1.axhline(0, color=GRIS, linewidth=0.8)
    ax1.set_title("Primeros taps")
    ax1.set_xlabel("Índice del tap k"); ax1.set_ylabel("h[k]")
    ax1.legend()

    # Derecha: el tamaño de cada tap en dB, para ver la cola
    n_ver = 120
    ax2.plot(k[:n_ver], 20 * np.log10(np.abs(h_larga[p_lenta, :n_ver]) + 1e-12), color=AZUL)
    ax2.plot(k[:n_ver], 20 * np.log10(np.abs(h_larga[p_rapida, :n_ver]) + 1e-12), color=NARANJA)
    ax2.axvline(L_ELEGIDO, color=GRIS, linestyle="--", linewidth=1.2)
    ax2.axvspan(L_ELEGIDO, n_ver, color="#e4e3df", alpha=0.6)
    ax2.text(L_ELEGIDO + 2, -15, f"Cola descartada\n(k ≥ L = {L_ELEGIDO})", color=GRIS)
    ax2.set_ylim(-160, 5); ax2.set_xlim(0, n_ver)
    ax2.set_title("Tamaño de los taps: la cola que se corta")
    ax2.set_xlabel("Índice del tap k"); ax2.set_ylabel("20·log10|h[k]| [dB]")

    fig.suptitle("Respuesta al impulso: cuanto más cerca de Nyquist corta, más larga es la cola", y=1.02)
    guardar(fig, "F3_taps")


# ===========================================================================
# F4. Efecto de L
# ===========================================================================
def grafico_efecto_L():
    # Tres cortes: el borde inferior, el nominal y el borde superior (el peor caso)
    cortes = [
        (B * (1 - desajuste), AZUL, "-"),
        (B, GRIS, "--"),
        (B_PEOR, NARANJA, "-"),
    ]
    valores_de_L = list(range(5, 122, 2))
    f = np.linspace(0, B, 800)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))

    # Izquierda: energía perdida contra L para cada corte
    for B_corte, color, estilo in cortes:
        h_larga = kernel_butterworth(B_corte, orden, 2000)
        E_total = np.sum(h_larga ** 2)
        energia_db = [10 * np.log10(np.sum(h_larga[L_prueba:] ** 2) / E_total) for L_prueba in valores_de_L]
        L_necesario = largo_minimo(B_corte, orden, UMBRAL_DB)
        ax1.plot(valores_de_L, energia_db, linestyle=estilo, color=color,
                 label=f"Corte en {B_corte/1e9:.1f} GHz: L = {L_necesario}")
    ax1.axhline(UMBRAL_DB, color=GRIS, linestyle=":", linewidth=1.2)
    ax1.text(valores_de_L[-1], UMBRAL_DB + 2, f"Umbral ENOB {BITS_ENOB} bits ({UMBRAL_DB:.0f} dB)", color=GRIS, ha="right")
    ax1.set_ylim(-120, 0)
    ax1.set_title("Energía perdida por truncamiento")
    ax1.set_xlabel("L (coeficientes por kernel)"); ax1.set_ylabel("Energía perdida [dB]")
    ax1.legend(loc="upper right")

    # Derecha: error de |H| del peor caso contra L
    h_peor = kernel_butterworth(B_PEOR, orden, 2000)
    w, H = signal.freqz(h_peor, 1, worN=f, fs=FS)
    H_larga_db = 20 * np.log10(np.abs(H))
    error_db = []
    for L_prueba in valores_de_L:
        w, H = signal.freqz(h_peor[:L_prueba], 1, worN=f, fs=FS)
        error_db.append(np.max(np.abs(20 * np.log10(np.abs(H)) - H_larga_db)))
    ax2.semilogy(valores_de_L, error_db, color=NARANJA)
    ax2.axvline(L_ELEGIDO, color=GRIS, linestyle=":", linewidth=1.2)
    ax2.text(L_ELEGIDO + 2, max(error_db) / 3, f"L = {L_ELEGIDO}", color=GRIS)
    ax2.set_title(f"Error máximo de |H| hasta 30 GHz (corte en {B_PEOR/1e9:.1f} GHz)")
    ax2.set_xlabel("L (coeficientes por kernel)"); ax2.set_ylabel("|H truncado| − |H largo| [dB]")

    fig.suptitle(f"El L lo fija la fase que corta más alto (B0 = {B/1e9:.0f} GHz, desajuste ±{desajuste*100:.0f} %)", y=1.02)
    guardar(fig, "F4_efecto_L")


# ===========================================================================
# F5. Offset con entrada nula
# ===========================================================================
def grafico_offset():
    u = random_offset(LSB, semilla_o)
    h = np.zeros((N_FASES, L_ELEGIDO))
    x = np.zeros(3 * N_FASES)                         # tres períodos
    r = planta_con_offset(x, h, u)
    n = np.arange(len(r))

    fig, ax = plt.subplots(figsize=(12, 3.8))
    ax.vlines(n, 0, r / LSB, color=AZUL, linewidth=1.6)
    ax.plot(n, r / LSB, "o", color=AZUL, markersize=6)
    for borde in (N_FASES, 2 * N_FASES):
        ax.axvline(borde - 0.5, color=GRIS, linestyle="--", linewidth=1)
    ax.axhline(1, color=GRIS, linestyle=":", linewidth=1)
    ax.axhline(-1, color=GRIS, linestyle=":", linewidth=1)
    ax.axhline(0, color=GRIS, linewidth=0.8)
    ax.text(len(r) - 0.5, 1.05, "+1 LSB", color=GRIS, ha="right", va="bottom")
    ax.text(len(r) - 0.5, -1.05, "−1 LSB", color=GRIS, ha="right", va="top")
    for i in range(3):
        ax.text(i * N_FASES + N_FASES / 2 - 0.5, 1.35, f"período {i + 1}", color=GRIS, ha="center")
    ax.set_ylim(-1.5, 1.6)
    ax.set_xlabel("Muestra n"); ax.set_ylabel("Salida r[n] [LSB]")
    ax.set_title("Offset por fase: con entrada nula, la salida es u[n mod 16] y se repite cada 16 muestras")
    guardar(fig, "F5_offset_entrada_nula")


# ===========================================================================
# F6. Espectro con un tono a la entrada
# ===========================================================================
def grafico_espectro_tono():
    # Muestreo coherente: el tono y sus copias caen justo en un bin de la FFT,
    # así se ven como líneas limpias sin ventana.
    N_FFT = 16 * 1024
    bin_tono = 1138
    f0 = bin_tono * FS / N_FFT                        # ≈ 5 GHz
    transitorio = 200                                 # muestras que se descartan al principio
    n = np.arange(N_FFT + transitorio)
    x = 0.5 * np.cos(2 * np.pi * f0 / FS * n)

    h_nominal = np.tile(kernel_butterworth(B, orden, L_ELEGIDO), (N_FASES, 1))
    h_desajuste, B_fases = random_butterworth(B, orden, L_ELEGIDO, desajuste, semilla_bw)
    u = random_offset(LSB, semilla_o)

    f = np.fft.rfftfreq(N_FFT, 1 / FS)
    casos = [
        ("Canal nominal (16 fases iguales, sin offset)", h_nominal, np.zeros(N_FASES), GRIS),
        (f"Canal con desajuste ±{desajuste*100:.0f} % y offset ±1 LSB", h_desajuste, u, AZUL),
    ]

    fig, ejes = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    for ax, (titulo, h, u_caso, color) in zip(ejes, casos):
        r = planta_con_offset(x, h, u_caso)[transitorio:]
        R = np.abs(np.fft.rfft(r))
        R_dbc = 20 * np.log10(R / R[bin_tono] + 1e-12)      # relativo al tono: dBc
        hay_tono = R_dbc > -110                              # bins con un tono o una espuria
        ax.vlines(f[hay_tono] / 1e9, -125, R_dbc[hay_tono], color=color, linewidth=1.4)
        ax.plot(f[hay_tono] / 1e9, R_dbc[hay_tono], "o", color=color, markersize=3)

        # Dónde se esperan las espurias: ±f0 + k·4,5 GHz (copias del tono) y k·4,5 GHz (offset)
        for k in range(0, 9):
            ax.axvline(k * FS / N_FASES / 1e9, color="#e4e3df", linewidth=1, zorder=0)

        # Marcar la espuria más grande
        R_sin_tono = R_dbc.copy(); R_sin_tono[bin_tono] = -999
        k_max = np.argmax(R_sin_tono)
        if R_sin_tono[k_max] > -100:
            ax.annotate(f"{R_sin_tono[k_max]:.1f} dBc en {f[k_max]/1e9:.1f} GHz",
                        xy=(f[k_max] / 1e9, R_sin_tono[k_max]), xytext=(f[k_max] / 1e9 + 1.5, R_sin_tono[k_max] + 12),
                        color=GRIS, arrowprops=dict(arrowstyle="->", color=GRIS))
        ax.set_ylim(-125, 10)
        ax.set_ylabel("Espectro [dBc]")
        ax.set_title(titulo)

    ejes[0].text(f0 / 1e9 + 0.4, 0, f"tono f0 = {f0/1e9:.2f} GHz", color=GRIS, va="top")
    ejes[1].set_xlabel("Frecuencia [GHz]  (líneas verticales: múltiplos de fs/16 = 4,5 GHz)")
    ejes[1].set_xlim(0, 36)
    fig.suptitle("Un tono a la entrada: el canal periódico genera copias en ±f0 + k·4,5 GHz y el offset, líneas en k·4,5 GHz", y=1.0)
    guardar(fig, "F6_espectro_tono")


# ===========================================================================
# Señal de prueba de banda limitada (para F7 y F8)
# ===========================================================================
def senal_banda_limitada(n_muestras, f_max=25e9):
    """Ruido uniforme filtrado hasta f_max: una señal de banda ancha dentro de la banda útil."""
    blanco = np.random.default_rng(semilla_x).uniform(-RANGO / 2, RANGO / 2, n_muestras + 200)
    filtro = signal.firwin(201, f_max, fs=FS)
    return np.convolve(blanco, filtro)[200:200 + n_muestras]


# ===========================================================================
# F7. Señal por el canal nominal y por el desajustado
# ===========================================================================
def grafico_senal_tiempo():
    x = senal_banda_limitada(2000)
    h_nominal = np.tile(kernel_butterworth(B, orden, L_ELEGIDO), (N_FASES, 1))
    h_desajuste, B_fases = random_butterworth(B, orden, L_ELEGIDO, desajuste, semilla_bw)
    u = random_offset(LSB, semilla_o)

    r_nominal = planta_fir_periodica(x, h_nominal)
    r_desajuste = planta_con_offset(x, h_desajuste, u)
    error = r_desajuste - r_nominal
    estable = slice(L_ELEGIDO, None)                  # sin el transitorio del principio
    snr_db = 10 * np.log10(np.mean(r_nominal[estable] ** 2) / np.mean(error[estable] ** 2))

    ventana = np.arange(400, 464)                     # 64 muestras = 4 períodos de 16
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax1.plot(ventana, r_nominal[ventana], "o-", color=GRIS, markersize=4, label="Canal nominal")
    ax1.plot(ventana, r_desajuste[ventana], "o-", color=AZUL, markersize=4, label=f"Canal con desajuste ±{desajuste*100:.0f} % y offset")
    ax1.set_ylabel("Salida r[n]")
    ax1.legend(loc="upper right")
    ax1.set_title("La misma señal (banda limitada a 25 GHz) por los dos canales")

    ax2.vlines(ventana, 0, error[ventana], color=NARANJA, linewidth=1.6)
    ax2.plot(ventana, error[ventana], "o", color=NARANJA, markersize=4)
    ax2.axhline(0, color=GRIS, linewidth=0.8)
    for borde in range(ventana[0], ventana[-1] + 2, N_FASES):
        for ax in (ax1, ax2):
            ax.axvline(borde - 0.5, color="#e4e3df", linewidth=1, zorder=0)
    ax2.set_xlabel("Muestra n  (líneas verticales: cada 16 muestras)")
    ax2.set_ylabel("Error")
    ax2.set_title(f"Diferencia entre los dos: SNR respecto del nominal = {snr_db:.1f} dB")
    guardar(fig, "F7_senal_nominal_vs_desajuste")


# ===========================================================================
# F8. SNR del canal respecto del nominal contra el desajuste
# ===========================================================================
def grafico_snr_vs_desajuste():
    x = senal_banda_limitada(4000)
    h_nominal = np.tile(kernel_butterworth(B, orden, L_ELEGIDO), (N_FASES, 1))
    r_nominal = planta_fir_periodica(x, h_nominal)
    estable = slice(L_ELEGIDO, None)
    P_senal = np.mean(r_nominal[estable] ** 2)

    desajustes = [0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15]
    sorteos = range(10)                               # promedio sobre 10 sorteos distintos
    snr_media, snr_min, snr_max = [], [], []
    for d in desajustes:
        snr = []
        for s in sorteos:
            h, B_fases = random_butterworth(B, orden, L_ELEGIDO, d, s)
            error = planta_fir_periodica(x, h) - r_nominal
            snr.append(10 * np.log10(P_senal / np.mean(error[estable] ** 2)))
        snr_media.append(np.mean(snr)); snr_min.append(np.min(snr)); snr_max.append(np.max(snr))
        print(f"desajuste ±{d*100:4.1f} %: SNR media {np.mean(snr):5.1f} dB")

    # Offset solo: el error es directamente u[n mod 16]
    snr_offset = np.mean([10 * np.log10(P_senal / np.mean(random_offset(LSB, s) ** 2)) for s in sorteos])
    snr_enob = 6.02 * BITS_ENOB + 1.76

    d_pct = np.array(desajustes) * 100
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.fill_between(d_pct, snr_min, snr_max, color=AZUL, alpha=0.15, linewidth=0, label="Rango entre sorteos")
    ax.plot(d_pct, snr_media, "o-", color=AZUL, markersize=6, label="Desajuste de ancho de banda (media)")
    ax.axhline(snr_offset, color=NARANJA, linestyle="--", linewidth=1.4, label=f"Solo offset ±1 LSB ({snr_offset:.0f} dB)")
    ax.axhline(snr_enob, color=GRIS, linestyle=":", linewidth=1.4, label=f"Techo del DAC, ENOB {BITS_ENOB} bits ({snr_enob:.0f} dB)")
    ax.set_xlabel("Desajuste de ancho de banda [± %]")
    ax.set_ylabel("SNR respecto del nominal [dB]")
    ax.set_ylim(0, 50)
    ax.legend(loc="lower left")
    ax.set_title("Cuánto se aleja el canal del nominal (señal limitada a 25 GHz)")
    guardar(fig, "F8_snr_vs_desajuste")


if __name__ == "__main__":
    grafico_magnitud()
    grafico_fase()
    grafico_taps()
    grafico_efecto_L()
    grafico_offset()
    grafico_espectro_tono()
    grafico_senal_tiempo()
    grafico_snr_vs_desajuste()
    if MOSTRAR:
        plt.show()
