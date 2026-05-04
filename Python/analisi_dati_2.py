import numpy as np

g = 9.80655

# ===================
# 1. CARICAMENTO DATI
# ===================
df = np.loadtxt('../dati_pendolo_semplice_2.csv', delimiter = ',')

lunghezze = df[0] / 100
omega = df[1]
err_frequenze = df[2]

fps = 60
ris = 1 / fps

periodo = 2 * np.pi / omega
err_periodo = 2 * np.pi / omega**2 * err_frequenze
periodo_teo = 2 * np.pi * np.sqrt(lunghezze / g)
