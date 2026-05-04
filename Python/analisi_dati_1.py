import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ===================
# 1. CARICAMENTO DATI
# ===================
df = np.loadtxt('../dati_pendolo_semplice_1.csv', delimiter=',')

x = df[:, 0]
y = df[:, 1]
err_ris_x = 0.1/1000;
sig_y = df[:, 2]
sig_x = err_ris_x/np.sqrt(12);
w = 1 / (sig_y**2)


# =============================================
# 2. MODELLO 1: PERIODO COSTANTE (Media Pesata)
# =============================================
media_pesata = np.sum(y * w) / np.sum(w)
incertezza_media = 1 / np.sqrt(np.sum(w))

print("=== MODELLO 1: PERIODO COSTANTE ===")
print(f"Media pesata del periodo: {media_pesata:.5f} ± {incertezza_media:.5f} s")

chi2_costante = np.sum(((y - media_pesata) / sig_y)**2)
dof_costante = len(df) - 1
p_value_costante = 1 - stats.chi2.cdf(chi2_costante, dof_costante)

print(f"Chi-quadrato osservato: {chi2_costante:.3f}")
print(f"Gradi di libertà (dof): {dof_costante}")
print(f"P-value: {p_value_costante:.4f}")
if p_value_costante > 0.05:
    print("Esito: Il modello a periodo costante è ACCETTABILE.\n")
else:
    print("Esito: Il modello a periodo costante NON è accettabile.\n")

# ===============================================
# 3. MODELLO 2: REGRESSIONE LINEARE (T = A + B*m)
# ===============================================
S   = np.sum(w)
Sx  = np.sum(w * x)
Sy  = np.sum(w * y)
Sxx = np.sum(w * x**2)
Sxy = np.sum(w * x * y)

Delta = S * Sxx - Sx**2

A = (Sxx * Sy - Sx * Sxy) / Delta
B = (S * Sxy - Sx * Sy) / Delta
sigma_A = np.sqrt(Sxx / Delta)
sigma_B = np.sqrt(S / Delta)
B1 = B + 2 * sigma_B

while abs(B - B1) > sigma_B:
    sig_y_i = np.sqrt(sig_y**2 + (B * sig_x)**2);

    W = 1 / (sig_y_i**2);
    S_W   = sum(W);
    S_XW  = sum(x * W);
    S_YW  = sum(y * W);
    S_XXW = sum(x**2 * W);
    S_XYW = sum(x * y * W);

    D_W = S_W * S_XXW - S_XW**2;

    sigma_B = np.sqrt(S_W / D_W);
    B1 = B;
    B = (1 / D_W) * (S_W * S_XYW - S_XW * S_YW);

A = (1 / D_W) * (S_XXW * S_YW - S_XW * S_XYW);

sigma_A = np.sqrt(S_XXW / D_W);
sigma_B = np.sqrt(S_W / D_W);

print("=== MODELLO 2: REGRESSIONE LINEARE ===")
print(f"Intercetta A: {A:.6f} ± {sigma_A:.6f} s")
print(f"Pendenza B:   {B:.8e} ± {sigma_B:.8e} s/g")

t_student = abs(B) / sigma_B
print(f"Compatibilità della pendenza con 0: {t_student:.2f} deviazioni standard")
if t_student < 3:
    print("Esito: La pendenza è compatibile con zero (T indipendente dalla massa).")
    limite_sup = abs(B) + 3 * sigma_B
    print(f"Limite superiore |B| + 3σ_B: {limite_sup:.2e} s/g")
else:
    print("Esito: La pendenza NON è statisticamente compatibile con zero.")

y_model_lineare = A + B * x
chi2_lineare = np.sum(((y - y_model_lineare) / sig_y)**2)
dof_lineare = len(df) - 2
p_value_lineare = 1 - stats.chi2.cdf(chi2_lineare, dof_lineare)

print(f"\nChi-quadrato modello lineare: {chi2_lineare:.3f}")
print(f"Gradi di libertà (dof): {dof_lineare}")
print(f"P-value: {p_value_lineare:.4f}")
if p_value_lineare > 0.05:
    print("Esito: Il modello lineare è ACCETTABILE.\n")
else:
    print("Esito: Il modello lineare NON è accettabile.\n")

# =================
# 4. GRAFICO FINALE
# =================
plt.figure(figsize=(10, 6))

plt.errorbar(x, y, yerr=sig_y, fmt='o', capsize=5, color='black',
             label='Dati sperimentali')

plt.axhline(y=float(media_pesata), color='green', linestyle='--', alpha=0.7,
            label=f'Media pesata ({media_pesata:.4f} s)')

x_range = np.linspace(x.min(), x.max(), 100)
y_fit = A + B * x_range
plt.plot(x_range, y_fit, color='red', alpha=0.7,
         label=f'Fit lineare (B = {B:.2e} s/g)')

plt.xlabel('Massa [g]')
plt.ylabel('Periodo [s]')
plt.title('Dipendenza del periodo dalla massa del pendolo')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
