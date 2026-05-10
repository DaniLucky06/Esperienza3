import sys, os
dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(dir, "analisi-dati-python"))
from analysis import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

df = np.loadtxt(os.path.join(dir, '../dati_pendolo_semplice_2.csv'), delimiter=',')

fps = 60
ris = 0.07 / fps
err_ris_t = ris / np.sqrt(12)

def f(t, a, b, c, d, e):
    return a + b * np.exp(c * t) * np.sin(d * t + e)

filePaths = [os.path.join(dir, f) for f in os.listdir(os.path.join(dir, "../videoData")) if os.path.isfile(os.path.join(dir, f))]
omegas = []
for i, fPath in filePaths:
    data = np.loadtxt(fPath)
    t = data[:, 0]
    y = data[:, 1]
    a, b, c, d, e = opt.curve_fit(f, t, y)
    err_y = err_ris_t * np.abs(b * d * np.cos(d * t + e))
    a, b, c, d, e = opt.curve_fit(f, t, y, sigma=err_y, absolute_sigma=True)
    omegas.append(d)
print(omegas)

df = df[:, np.argsort(df[0])]
lunghezze = np.log(df[0] / 100)
omega = df[1]

periodi = 2 * np.pi / omega
dev_std_omega = df[2]
dev_std_periodi = 2 * np.pi / omega**2 * dev_std_omega
dev_std_periodi = np.abs(1 / periodi) * np.sqrt(dev_std_periodi**2 + err_ris_t**2)
periodi = np.log(periodi)

err_lunghezze = 0.001 / np.sqrt(12)
dev_std_lunghezze = np.ones(len(lunghezze)) * err_lunghezze

measurements = WeightedMeasurements(lunghezze, periodi, dev_std_lunghezze, dev_std_periodi)
reg_data, xy_data = weightedLinReg(measurements)
a = reg_data.a
b = reg_data.b

# plotGraphs([measurements, xy_data], GraphVisuals(['rx','b']))

measurements.y = measurements.y - (a + measurements.x * b)
xy_data.y = xy_data.y - (a + xy_data.x * b)


# plt.figure()
endIndex = -1
reg_data_abs, _ = weightedLinReg(WeightedMeasurements(lunghezze[0:endIndex], periodi[0:endIndex], dev_std_lunghezze[0:endIndex], dev_std_periodi[0:endIndex]))

x = np.linspace(lunghezze[0], lunghezze[-1], 10)
# plt.plot(x, reg_data.a + reg_data.b * x)
print(reg_data.chi2, reg_data.nu)
for i in range(1, len(periodi) - 5):
    reg_data, _ = weightedLinReg(WeightedMeasurements(lunghezze[i:endIndex], periodi[i:endIndex], dev_std_lunghezze[i:endIndex], dev_std_periodi[i:endIndex]))
    valido = reg_data.chi2 / reg_data.nu < 1.2
    style = 'b' if valido else 'r'
    print(reg_data.chi2, reg_data.nu, "Valido: ", valido)
    if valido:
        # plt.plot(x, reg_data.a + reg_data.b * x - (reg_data_abs.a + reg_data_abs.b * x), style)
        index = i
        break
# plt.errorbar(lunghezze, periodi - (reg_data_abs.a + reg_data_abs.b * lunghezze), dev_std_periodi, dev_std_lunghezze, 'kx')
# plt.show()

plt.figure()
plt.hlines(0, x[0], x[-1], 'b')
plt.errorbar(lunghezze[index:-1], periodi[index:-1] - (reg_data.a + reg_data.b * lunghezze[index:-1]), dev_std_periodi[index:-1], dev_std_lunghezze[index:-1], 'gx')
plt.errorbar(lunghezze[0:index], periodi[0:index] - (reg_data.a + reg_data.b * lunghezze[0:index]), dev_std_periodi[0:index], dev_std_lunghezze[0:index], 'rx')
plt.show()
