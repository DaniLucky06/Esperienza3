import numpy as np
import matplotlib.pyplot as plt

# ===================
# 1. CARICAMENTO DATI
# ===================
df = np.loadtxt('../dati_pendolo_semplice_2.csv', delimiter=',')

lunghezze = df[0] / 100
omega = df[1]
dev_std_omega = df[2]

g = 9.80655

fps = 60
ris = 1 / fps
err_ris = ris / np.sqrt(12)

periodi = 2 * np.pi / omega
dev_std_periodo = 2 * np.pi / omega**2 * dev_std_omega
periodo_teo = 2 * np.pi * np.sqrt(lunghezze / g)

err_periodi = np.sqrt(err_ris**2 + dev_std_periodo**2)

ris_lunghezze = 0.001
err_lunghezze = (ris_lunghezze / np.sqrt(12)) * np.ones_like(lunghezze)

# ==========
# 2. GRAFICO
# ==========
plt.figure(figsize=(10, 6))

plt.errorbar(lunghezze, periodi, 
             xerr=err_lunghezze, yerr=err_periodi, 
             fmt='o', color='blue', ecolor='red', capsize=3, 
             markersize=5, label='Dati sperimentali')

indici_ordinati = np.argsort(lunghezze)
plt.plot(lunghezze[indici_ordinati], periodo_teo[indici_ordinati], 
         '--', color='gray', label=r'Andamento teorico $\mathcal{T} \propto \ell^{1/2}$')

plt.title('Periodo di oscillazione in funzione della lunghezza')
plt.xlabel(r'Lunghezza $\ell$ [m]')
plt.ylabel(r'Periodo $\mathcal{T}$ [s]')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()

plt.tight_layout()
plt.show()

# lunghezze_log = np.log(lunghezze)
# periodi_log = np.log(periodi)

plt.errorbar(lunghezze, periodi, 
             xerr=err_lunghezze, yerr=err_periodi, 
             fmt='o', color='blue', ecolor='red', capsize=3, 
             markersize=5, label='Dati sperimentali')

plt.plot(lunghezze[indici_ordinati], periodo_teo[indici_ordinati], '--', color='gray', label=r'Andamento teorico $\mathcal{T} \propto \ell^{1/2}$')

plt.yscale('log')
plt.xscale('log')


plt.title('Periodo di oscillazione in scala logaritmica in funzione della lunghezza')
plt.xlabel(r'Lunghezza $\ell$ [m]')
plt.ylabel(r'Periodo $\mathcal{T}$ [s]')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()

plt.tight_layout()
plt.show()
