"""
sim_nmse_barridos.py — Barridos de NMSE del Cal ADC.

Varía un parámetro por vez (ganancia, offset, ancho de banda), dejando los
demás en su valor ideal, y compara el NMSE simulado con la teoría:
    ganancia:        NMSE = 20·log10|G − 1|
    offset:          NMSE = 10·log10(V_os² / P_x)
    ancho de banda:  NMSE = 20·log10( (f/BW) / sqrt(1 + (f/BW)²) )

Genera tres figuras con el mismo eje vertical para comparar la importancia
de cada efecto.

Uso (desde la raíz del repositorio):
    python modules/cal_adc/simulations/sim_nmse_barridos.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Modelo del Cal ADC (modules/cal_adc/python/cal_adc.py)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from cal_adc import FS, M, L, adc_cal  # noqa: E402

# ===========================================================================
# Barridos de NMSE: se varía UN parámetro por vez (los otros quedan ideales)
#
#   NMSE = mean( (y_cal[m] − r(t_m))² ) / mean( r(t_m)² )      [dB]
#
# Los tres gráficos usan el mismo eje vertical para poder comparar la
# importancia de cada efecto.
# ===========================================================================
A_BARRIDO      = 0.8            # amplitud del tono [V] (deja margen hasta V_FS)
F_BARRIDO      = 10e9           # frecuencia del tono [Hz]
V_FS_BARRIDO   = 1.0            # rango de entrada ±V_FS [V]
PUNTOS_BARRIDO = 2_000_000      # ≈ 650 disparos del ADC

GANANCIAS    = [0.9, 1.1, 1.2, 1.3]                 # [-]
OFFSETS      = [0.025, 0.05, 0.1, 0.2, 0.25]        # [V]
ANCHOS_BANDA = [50e9, 30e9, 15e9, 10e9, 5e9]        # [Hz]

NMSE_LIM = (-40, 0)             # eje vertical común [dB]


def nmse_db(y, ref):
    """Error cuadrático medio normalizado a la potencia de la referencia [dB]."""
    return 10 * np.log10(np.mean((y - ref) ** 2) / np.mean(ref ** 2))


def medir_nmse(r, t_fina, ganancia=1.0, offset=0.0, ancho_banda=None):
    """Pasa r por el ADC con un solo efecto activo y devuelve (NMSE [dB], recortes)."""
    y, n_m = adc_cal(r, M=M, L=L, ganancia=ganancia, offset=offset,
                     ancho_banda=ancho_banda, V_FS=V_FS_BARRIDO)
    t_m = n_m / FS
    ref = A_BARRIDO * np.sin(2 * np.pi * F_BARRIDO * t_m)
    # descarta el transitorio inicial del filtro (primeros 10·tau)
    tau = 0.0 if ancho_banda is None else 1 / (2 * np.pi * ancho_banda)
    ok = t_m > 10 * tau
    recortes = int(np.sum(np.abs(y) >= V_FS_BARRIDO))
    return nmse_db(y[ok], ref[ok]), recortes


if __name__ == "__main__":
    t_b = np.arange(PUNTOS_BARRIDO) / (L * FS)
    r_b = A_BARRIDO * np.sin(2 * np.pi * F_BARRIDO * t_b)
    P_x = A_BARRIDO ** 2 / 2                         # potencia del seno

    print("\n=========== Barridos de NMSE ===========")
    print(f"Tono {F_BARRIDO/1e9:g} GHz, A = {A_BARRIDO} V, V_FS = ±{V_FS_BARRIDO} V")

    # ---------------- Ganancia ----------------
    res_g = [medir_nmse(r_b, t_b, ganancia=g) for g in GANANCIAS]
    print("\nGanancia:   G      NMSE sim   NMSE teo   recortes")
    for g, (nm, rc) in zip(GANANCIAS, res_g):
        print(f"          {g:5.2f}   {nm:8.2f}   {20*np.log10(abs(g-1)):8.2f}   {rc:5d}")

    g_teo = np.linspace(0.7, 1.4, 701)
    g_teo = g_teo[np.abs(g_teo - 1) > 1e-6]
    G_clip = V_FS_BARRIDO / A_BARRIDO                # a partir de acá el pico recorta

    fig_g, ax = plt.subplots(figsize=(8, 5))
    ax.axvspan(G_clip, g_teo.max(), color="C3", alpha=0.08, label=f"zona de recorte (G·A > V_FS)")
    ax.plot(g_teo, 20 * np.log10(np.abs(g_teo - 1)), color="C0", lw=1.5,
            label="teoría: 20·log₁₀|G − 1|")
    ax.plot(GANANCIAS, [nm for nm, _ in res_g], "o", color="C3", ms=8, label="simulación")
    ax.axvline(1.0, color="gray", lw=0.8, ls=":")
    ax.set_xlabel("ganancia G"); ax.set_ylabel("NMSE [dB]")
    ax.set_ylim(*NMSE_LIM); ax.set_xlim(g_teo.min(), g_teo.max())
    ax.set_title("NMSE vs. ganancia del ADC")
    ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=9)
    fig_g.tight_layout()

    # ---------------- Offset ----------------
    res_o = [medir_nmse(r_b, t_b, offset=o) for o in OFFSETS]
    print("\nOffset:   V_os [mV]   NMSE sim   NMSE teo   recortes")
    for o, (nm, rc) in zip(OFFSETS, res_o):
        print(f"          {o*1e3:7.1f}   {nm:8.2f}   {10*np.log10(o**2/P_x):8.2f}   {rc:5d}")

    o_teo = np.linspace(0.005, 0.3, 300)
    o_clip = V_FS_BARRIDO - A_BARRIDO                # a partir de acá el pico positivo recorta

    fig_o, ax = plt.subplots(figsize=(8, 5))
    ax.axvspan(o_clip * 1e3, o_teo.max() * 1e3, color="C3", alpha=0.08,
               label="zona de recorte (A + V_os > V_FS)")
    ax.plot(o_teo * 1e3, 10 * np.log10(o_teo ** 2 / P_x), color="C0", lw=1.5,
            label="teoría: 10·log₁₀(V_os² / P_x)")
    ax.plot(np.array(OFFSETS) * 1e3, [nm for nm, _ in res_o], "o", color="C3", ms=8,
            label="simulación")
    ax.set_xlabel("offset V_os [mV]"); ax.set_ylabel("NMSE [dB]")
    ax.set_ylim(*NMSE_LIM); ax.set_xlim(0, o_teo.max() * 1e3)
    ax.set_title("NMSE vs. offset del ADC")
    ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=9)
    fig_o.tight_layout()

    # ---------------- Ancho de banda ----------------
    res_bw = [medir_nmse(r_b, t_b, ancho_banda=bw) for bw in ANCHOS_BANDA]
    x_sim = F_BARRIDO / np.array(ANCHOS_BANDA)       # f / BW
    print("\nAncho de banda:  BW [GHz]   f/BW   atenuación   NMSE sim   NMSE teo")
    for bw, x, (nm, _) in zip(ANCHOS_BANDA, x_sim, res_bw):
        print(f"                 {bw/1e9:6.0f}   {x:5.2f}   {-10*np.log10(1+x*x):8.2f} dB"
              f"   {nm:8.2f}   {20*np.log10(x/np.sqrt(1+x*x)):8.2f}")

    x_teo = np.logspace(-2, np.log10(3), 400)
    fig_b, ax = plt.subplots(figsize=(8, 5))
    ax.semilogx(x_teo, 20 * np.log10(x_teo / np.sqrt(1 + x_teo ** 2)), color="C0", lw=1.5,
                label="teoría NMSE: 20·log₁₀|H − 1|")
    ax.semilogx(x_teo, -10 * np.log10(1 + x_teo ** 2), color="gray", lw=1.2, ls="--",
                label="atenuación |H| [dB] (referencia)")
    ax.semilogx(x_sim, [nm for nm, _ in res_bw], "o", color="C3", ms=8, label="simulación")
    for bw, x, (nm, _) in zip(ANCHOS_BANDA, x_sim, res_bw):
        ax.annotate(f"{bw/1e9:g} GHz", (x, nm), textcoords="offset points",
                    xytext=(6, -12), fontsize=8, color="C3")
    ax.set_xlabel("f / BW   (frecuencia del tono / ancho de banda del T&H)")
    ax.set_ylabel("NMSE [dB]")
    ax.set_ylim(*NMSE_LIM); ax.set_xlim(x_teo.min(), x_teo.max())
    ax.set_title(f"NMSE vs. ancho de banda del ADC (tono de {F_BARRIDO/1e9:g} GHz)")
    ax.grid(alpha=0.3, which="both"); ax.legend(loc="lower right", fontsize=9)
    fig_b.tight_layout()

    plt.show()
