import numpy as np

N_FASES = 16


def planta_fir_periodica(x, h):
    """Canal equivalente: FIR periodico de 16 fases (procesa la senal completa)."""
    x = np.asarray(x, dtype=float)
    h = np.asarray(h, dtype=float)

    if h.shape[0] != N_FASES:
        raise ValueError(f"h debe tener {N_FASES} filas, tiene {h.shape[0]}")

    N = len(x)
    L = h.shape[1]
    r = np.zeros(N)

    for n in range(N):
        p = n % N_FASES  # Calcular la fase p
        acc = 0.0  # Inicializar acumulador
        for k in range(L):
            if n - k >= 0:
                acc += h[p, k] * x[n - k]  # Acumular el resultado
        r[n] = acc

    return r


def planta_fir_periodica_bloque(x_bloque, h, memoria, n_inicio):
    """Procesa un bloque de la senal.

    memoria:  ultimas L-1 muestras del bloque anterior (al principio, ceros)
    n_inicio: indice global de la primera muestra del bloque (al principio, 0)

    Devuelve la salida del bloque, y la memoria y el indice para el bloque siguiente.
    """
    x_bloque = np.asarray(x_bloque, dtype=float)
    h = np.asarray(h, dtype=float)

    if h.shape[0] != N_FASES:
        raise ValueError(f"h debe tener {N_FASES} filas, tiene {h.shape[0]}")

    L = h.shape[1]

    # Pongo la memoria delante del bloque, asi x[n-k] siempre existe
    x_ext = np.concatenate([memoria, x_bloque])

    r_bloque = np.zeros(len(x_bloque))
    for i in range(len(x_bloque)):
        n = n_inicio + i  # Indice global: la fase no se reinicia en cada bloque
        p = n % N_FASES
        acc = 0.0
        for k in range(L):
            acc += h[p, k] * x_ext[i + L - 1 - k]  # Es x[n-k]
        r_bloque[i] = acc

    memoria_nueva = x_ext[len(x_bloque):]  # Las ultimas L-1 muestras
    n_siguiente = n_inicio + len(x_bloque)

    return r_bloque, memoria_nueva, n_siguiente