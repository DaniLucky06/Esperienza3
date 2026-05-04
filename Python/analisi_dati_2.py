import numpy as np
import matplotlib.pyplot as plt
import os

dir = os.path.dirname(os.path.realpath(__file__))

# ===================
# 1. CARICAMENTO DATI
# ===================
df = np.loadtxt(os.path.join(dir, '../dati_pendolo_semplice_2.csv'), delimiter=',')

lunghezze = df[0] / 100
omega = df[1]
dev_std_omega = df[2]

g = 9.80655

fps = 60
ris = 1 / fps
err_ris = ris / np.sqrt(12)

periodi = 2 * np.pi / omega
dev_std_periodo = 2 * np.pi / omega**2 * dev_std_omega
periodi_teo = 2 * np.pi * np.sqrt(lunghezze / g)

err_periodi = np.sqrt(err_ris**2 + dev_std_periodo**2)

ris_lunghezze = 0.001
err_lunghezze = (ris_lunghezze / np.sqrt(12)) * np.ones_like(lunghezze)

err_periodi_teo = np.pi / np.sqrt(lunghezze * g)

w = 1 / err_periodi**2

# ======================
# 2. REGRESSIONE LINEARE
# ======================
lunghezze_log = np.log(lunghezze)
periodi_log = np.log(periodi)
err_lunghezze_log = 1 / lunghezze * err_lunghezze
err_periodi_log = 1 / periodi * err_periodi

S   = np.sum(w)
Sx  = np.sum(w * lunghezze_log)
Sy  = np.sum(w * periodi_log)
Sxx = np.sum(w * lunghezze_log**2)
Sxy = np.sum(w * lunghezze_log * periodi_log)

Delta = S * Sxx - Sx**2

A = (Sxx * Sy - Sx * Sxy) / Delta
B = (S * Sxy - Sx * Sy) / Delta
sigma_A = np.sqrt(Sxx / Delta)
sigma_B = np.sqrt(S / Delta)
B1 = B + 2 * sigma_B

while abs(B - B1) > sigma_B:
    sig_y_i = np.sqrt(err_periodi_teo**2 + (B * err_lunghezze_log)**2);

    W = 1 / (sig_y_i**2);
    S_W   = sum(W);
    S_XW  = sum(lunghezze_log * W);
    S_YW  = sum(periodi_log * W);
    S_XXW = sum(lunghezze_log**2 * W);
    S_XYW = sum(lunghezze_log * periodi_log * W);

    D_W = S_W * S_XXW - S_XW**2;

    sigma_B = np.sqrt(S_W / D_W);
    B1 = B;
    B = (1 / D_W) * (S_W * S_XYW - S_XW * S_YW);

A = (1 / D_W) * (S_XXW * S_YW - S_XW * S_XYW);

sigma_A = np.sqrt(S_XXW / D_W);
sigma_B = np.sqrt(S_W / D_W);

# ==========
# 3. GRAFICI
# ==========
plt.figure(figsize=(10, 6))

plt.errorbar(lunghezze, periodi, 
             xerr=err_lunghezze, yerr=err_periodi, 
             fmt='o', color='blue', ecolor='red', capsize=3, 
             markersize=5, label='Dati sperimentali')

indici_ordinati = np.argsort(lunghezze)
plt.plot(lunghezze[indici_ordinati], periodi_teo[indici_ordinati], 
         '--', color='black', label=r'Andamento teorico $\mathcal{T} \propto \ell^{1/2}$')

plt.plot()

plt.title('Periodo di oscillazione in funzione della lunghezza')
plt.xlabel(r'Lunghezza $\ell$ [m]')
plt.ylabel(r'Periodo $\mathcal{T}$ [s]')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()

plt.tight_layout()
plt.show()

x_range = np.linspace(lunghezze_log.min(), lunghezze_log.max(), 100)
y_fit = A + B * x_range

plt.plot(x_range, y_fit, color='red', alpha=0.7,
         label=f'Fit lineare (B = {B:.2e} s/g)')

plt.errorbar(lunghezze_log, periodi_log, 
             xerr=err_lunghezze_log, yerr=err_periodi_teo, 
             fmt='o', color='blue', ecolor='red', capsize=3, 
             markersize=5, label='Dati sperimentali')

periodi_teo_log = np.log(periodi_teo)
plt.plot(lunghezze_log[indici_ordinati], periodi_teo_log[indici_ordinati], '--', color='black', label=r'Andamento teorico $\mathcal{T} \propto \ell^{1/2}$')

plt.title('Periodo di oscillazione in scala logaritmica in funzione della lunghezza')
plt.xlabel(r'Lunghezza $\ell$ [m]')
plt.ylabel(r'Periodo $\mathcal{T}$ [s]')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()

plt.tight_layout()
plt.show()
