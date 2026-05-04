import numpy as np

# ===================
# 1. CARICAMENTO DATI
# ===================
df = np.loadtxt('../dati_pendolo_semplice_2.csv', delimiter = ',')

lunghezze = df[0] / 100
omega = df[1]
dev_std_frequenze = df[2]

g = 9.80655
fps = 60
ris = 1 / fps
err_ris = ris / np.sqrt(12)

periodo = 2 * np.pi / omega
err_periodo = 2 * np.pi / omega**2 * dev_std_frequenze
periodo_teo = 2 * np.pi * np.sqrt(lunghezze / g)

# ===============================================
# 2. REGRESSIONE LINEARE (T = A + B*m)
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


