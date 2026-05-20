import sys, os
dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(dir, "analisi-dati-python"))
from analysis import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

df = np.loadtxt(os.path.join(dir, '../dati_pendolo_semplice_2.csv'), delimiter=',')

fps = 60
ris = 1 / fps
err_ris_t = ris / np.sqrt(12)
err_y_intrinsic = 0.0001 # for zero division

chi2_valido = 1.2

def f(t, a, b, c, d, e):
    return a + b * np.exp(c * t) * np.sin(d * t + e)
video_dir = os.path.join(dir, "../videoData")
filePaths = [os.path.join(video_dir, f) for f in os.listdir(video_dir) if os.path.isfile(os.path.join(video_dir, f))]
omegas = []
errs = []

for i, fPath in enumerate(filePaths):
    data = np.loadtxt(fPath)
    t = data[:, 0]
    y = data[:, 1]

    a_guess = np.mean(y)
    b_guess = (np.max(y) - np.min(y)) / 2.
    c_guess = -.1

    # stima la pulsazione
    zero_intercepts = np.where(np.diff(np.sign(y - a_guess)))[0] # trova le intersezioni con lo zero
    period_guess = 2 * np.mean(np.diff(t[zero_intercepts])) # media del tempo fra intersezioni
    d_guess = 2 * np.pi / period_guess

    p0 = [a_guess, b_guess, c_guess, d_guess, 0.0]

    popt, _ = opt.curve_fit(f, t, y, p0)
    a, b, c, d, e = popt

    dy_dt = b * np.exp(c * t) * (c * np.sin(d * t + e) + d * np.cos(d * t + e))
    err_y = np.sqrt(err_y_intrinsic**2 + (err_ris_t * dy_dt)**2)

    popt, pcov  = opt.curve_fit(f, t, y, p0=popt, sigma=err_y, absolute_sigma=True)

    omegas.append(popt[3])
    err = np.sqrt(np.diag(pcov))[3]
    errs.append(err)

    """
    plt.figure()
    plt.plot(t, y)
    plt.plot(t, f(t, a, b, c, d, e))
    plt.title(fPath)
    plt.show()
    """

omegas = np.array(omegas)
errs = np.array(errs)


df = df[:, np.argsort(df[0])]
lunghezze = np.log(df[0] / 100)
omega = df[1]

indices = np.argsort(omegas)
for i, j in enumerate(indices[::-1]):
    print(f"Omegas[j]: {omegas[j]}, omega[i]: {omega[i]}, diff: {omegas[j] - omega[i]}")

# omega = omegas[indices[::-1]]

periodi = 2 * np.pi / omega
dev_std_omega = df[2]

print((dev_std_omega - errs[indices[::-1]]) / dev_std_omega)
dev_std_omega = errs[indices[::-1]]
dev_std_periodi = 2 * np.pi / omega**2 * dev_std_omega
dev_std_periodi = np.abs(1 / periodi) * np.sqrt(dev_std_periodi**2 + err_ris_t**2)

dev_std_periodi = dev_std_omega / omega


periodi = np.log(periodi)

err_filo = 0.001 / np.sqrt(12)
err_baricentro = 0.001 / np.sqrt(12)
err_lunghezze = np.sqrt(err_filo**2 + err_baricentro**2)
dev_std_lunghezze = np.ones(len(lunghezze)) * err_lunghezze / (df[0] / 100)

measurements = WeightedMeasurements(lunghezze, periodi, dev_std_lunghezze, dev_std_periodi)
reg_data, xy_data = weightedLinReg(measurements)
a = reg_data.a
b = reg_data.b

# plotGraphs([measurements, xy_data], GraphVisuals(['rx','b']))

measurements.y = measurements.y - (a + measurements.x * b)
xy_data.y = xy_data.y - (a + xy_data.x * b)


# plt.figure()
endIndex = None
reg_data_abs, _ = weightedLinReg(WeightedMeasurements(lunghezze[0:endIndex], periodi[0:endIndex], dev_std_lunghezze[0:endIndex], dev_std_periodi[0:endIndex]))

index = 0
x = np.linspace(lunghezze[0], lunghezze[-1], 10)
# plt.plot(x, reg_data.a + reg_data.b * x)
print(reg_data.chi2, reg_data.nu)
for i in range(1, len(periodi) - 5):
    reg_data, _ = weightedLinReg(WeightedMeasurements(lunghezze[i:endIndex], periodi[i:endIndex], dev_std_lunghezze[i:endIndex], dev_std_periodi[i:endIndex]))
    valido = reg_data.chi2 / reg_data.nu < chi2_valido
    style = 'b' if valido else 'r'
    print(reg_data.chi2, reg_data.nu, reg_data.chi2 / reg_data.nu, "Valido:", valido)
    if valido:
        # plt.plot(x, reg_data.a + reg_data.b * x - (reg_data_abs.a + reg_data_abs.b * x), style)
        index = i
        break
# plt.errorbar(lunghezze, periodi - (reg_data_abs.a + reg_data_abs.b * lunghezze), dev_std_periodi, dev_std_lunghezze, 'kx')
# plt.show()

dev_std_effettiva = np.sqrt(dev_std_periodi**2 + (reg_data.b * dev_std_lunghezze)**2)
print(dev_std_effettiva / (reg_data.b * dev_std_lunghezze))

plt.figure()
plt.hlines(0, x[0], x[-1], 'b')
plt.errorbar(lunghezze, np.log(2 * np.pi / omegas[indices[::-1]]) - (reg_data.a + reg_data.b * lunghezze), yerr=dev_std_effettiva, fmt='bx')

plt.errorbar(lunghezze[index:-1], periodi[index:-1] - (reg_data.a + reg_data.b * lunghezze[index:-1]), yerr=dev_std_effettiva[index:-1], fmt='gx')
plt.errorbar(lunghezze[0:index], periodi[0:index] - (reg_data.a + reg_data.b * lunghezze[0:index]), yerr=dev_std_effettiva[0:index], fmt='rx')
plt.show()
