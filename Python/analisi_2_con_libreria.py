import sys, os
dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(dir, "analisi-dati-python"))
from analysis import *
import matplotlib.pyplot as plt
import numpy as np

df = np.loadtxt(os.path.join(dir, '../dati_pendolo_semplice_2.csv'), delimiter=',')

chi2_valido = 1.5

df = df[:, np.argsort(df[0])]
lunghezze_pre_log = df[0] / 100 - 0.0005
lunghezze = np.log(lunghezze_pre_log)
omega = df[1]

periodi = 2 * np.pi / omega
dev_std_omega = df[2]

dev_std_periodi = dev_std_omega / omega

periodi = np.log(periodi)

err_filo = 0.001 / np.sqrt(12)
err_baricentro = 0.001 / np.sqrt(12)
err_lunghezze = np.sqrt(err_filo**2 + err_baricentro**2)
dev_std_lunghezze = np.ones(len(lunghezze)) * err_lunghezze / (lunghezze_pre_log)

measurements = WeightedMeasurements(lunghezze, periodi, dev_std_lunghezze, dev_std_periodi)
reg_data, xy_data = weightedLinReg(measurements)
a = reg_data.a
b = reg_data.b

measurements.y = measurements.y - (a + measurements.x * b)
xy_data.y = xy_data.y - (a + xy_data.x * b)

endIndex = None
reg_data_abs, _ = weightedLinReg(WeightedMeasurements(lunghezze[0:endIndex], periodi[0:endIndex], dev_std_lunghezze[0:endIndex], dev_std_periodi[0:endIndex]))

index = 0
x = np.linspace(lunghezze[0], lunghezze[-1], 10)
print(reg_data.chi2, reg_data.nu)
for i in range(1, len(periodi) - 5):
    reg_data, _ = weightedLinReg(WeightedMeasurements(lunghezze[i:endIndex], periodi[i:endIndex], dev_std_lunghezze[i:endIndex], dev_std_periodi[i:endIndex]))
    valido = reg_data.chi2 / reg_data.nu < chi2_valido
    style = 'b' if valido else 'r'
    print(reg_data.chi2, reg_data.nu, reg_data.chi2 / reg_data.nu, "Valido:", valido)
    if valido:
        index = i
        break

dev_std_effettiva = np.sqrt(dev_std_periodi**2 + (reg_data.b * dev_std_lunghezze)**2)

### --- SEZIONE PENDOLO REALE --- ###
I_cm = 0.000115216532
m1 = 0.2003
m2 = 0.0844
m_pir = 0
M = m1 + m2
g = 9.80655

def T(L):
    I_pir = m_pir * (L - 0.0194)**2
    num = I_cm + M * L**2 + I_pir
    den = M * L * g
    return 2 * np.pi * np.sqrt(num / den)

periodi_reali = np.log(T(lunghezze_pre_log - 0.0005))
### --- FINE SEZIONE PENDOLO REALE --- ###

plt.figure()
plt.hlines(0, x[0], x[-1], 'b')

plt.errorbar(lunghezze[index:-1], periodi[index:-1] - (reg_data.a + reg_data.b * lunghezze[index:-1]), yerr=dev_std_effettiva[index:-1], fmt='gx')
plt.errorbar(lunghezze[0:index], periodi[0:index] - (reg_data.a + reg_data.b * lunghezze[0:index]), yerr=dev_std_effettiva[0:index], fmt='rx')
plt.plot(lunghezze, periodi_reali - (reg_data.a + reg_data.b * lunghezze))
plt.show()

a = reg_data.a
a_err = reg_data.a_err
a_real = np.exp(a)
g_experiment = (2 * np.pi / a_real) ** 2

err_a_real = np.exp(a) * a_err
err_g_experiment = 8 * np.pi ** 2 / a_real ** 3 * err_a_real
print(g_experiment, err_g_experiment)