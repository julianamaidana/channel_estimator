"""
Chequea que el ADC cal pase por todas las fases del canal (fase = n_m % P) con
M = 320 y M = 640. P = 16: cada sub-DAC vuelve a disparar cada 72 / 4.5 = 16
muestras.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from modules.cal_adc.python.adc_cal_ideal import adc_cal_ideal, N

L = 8
P = 16  
DISPAROS = 320     # 20 vueltas de 16 fases
AZUL, NARANJA = "#2a78d6", "#eb6834"

for M in (320, 640):
    r = np.zeros(DISPAROS * (round(M * N) + 1) * L)   # la señal no importa, solo los n_m
    _, n_m = adc_cal_ideal(r, M=M, L=L)
    fase = n_m % P
    avanza_1 = np.all(np.diff(fase) % P == 1)
    fases_vistas = len(np.unique(fase))
    print(f"M = {M}: paso {n_m[1] - n_m[0]}, {len(n_m)} disparos, "
          f"avanza de a una fase: {avanza_1}, fases visitadas: {fases_vistas}/{P}")
    assert avanza_1 and fases_vistas == P, "el ADC no recorre todas las fases"
print("OK")

# gráfico con M = 320, comparado con el mismo ADC sin el +1 (paso 384)
M = 320
r = np.zeros(DISPAROS * (round(M * N) + 1) * L)
_, n_m = adc_cal_ideal(r, M=M, L=L)
fase = n_m % P
fase_sin_mas_1 = (np.arange(DISPAROS) * round(M * N)) % P

fig, ejes = plt.subplots(1, 2, figsize=(11, 4.2))

m = np.arange(48)
ejes[0].plot(m, fase[:48], "o-", color=AZUL, ms=5, lw=1, label="con el +1 (paso 385)")
ejes[0].plot(m, fase_sin_mas_1[:48], "s", color=NARANJA, ms=5, label="sin el +1 (paso 384)")
ejes[0].set_title("Cada disparo avanza una fase", fontsize=10, loc="left")
ejes[0].set_xlabel("número de disparo")
ejes[0].set_ylabel(f"fase del canal (n_m % {P})")
ejes[0].set_yticks(range(0, P, 3))
ejes[0].set_ylim(-1, P + 3)   # lugar para la leyenda
ejes[0].legend(loc="upper left", fontsize=8, ncol=2)

ejes[1].bar(range(P), np.bincount(fase, minlength=P), color=AZUL)
ejes[1].set_title(f"Las {P} fases, parejas ({DISPAROS} disparos)", fontsize=10, loc="left")
ejes[1].set_xlabel(f"fase del canal (n_m % {P})")
ejes[1].set_ylabel("cantidad de disparos")
ejes[1].set_xticks(range(P))
ejes[1].set_yticks(range(0, DISPAROS // P + 1, 5))

for ax in ejes:
    ax.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
fig.tight_layout()
plt.show()