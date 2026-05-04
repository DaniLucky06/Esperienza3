import numpy as np

N = int(input("Numero di masse: "))
file = np.loadtxt("masse.csv", delimiter=",")

massa_pirulo = 1.2 # g
h_pirulo = 6.5 # mm
l_filo = 850 # mm
l_voluta = 859.5233598409543 # mm

masse = np.empty(N)
h = np.empty(N)
h_tot = 0
for i in range(N):
    indice_massa = int(input(f"Massa {i + 1}: ")) - 1
    masse[i] = file[0,indice_massa]
    h_i = file[1,indice_massa]
    h[i] = h_tot + h_i/2
    h_tot += h_i

h_mean = (sum(masse * h) + massa_pirulo * (h_tot + h_pirulo/2)) / (sum(masse) + massa_pirulo)
h_top = h_tot - h_mean
h_voluta = l_voluta - h_top - h_pirulo
print("Altezza dalla cima (no tappo): ", h_top)
print("Altezza cima-perno Giorno 2", 100 - h_top)
print("Altezza perno-cima della pila: ", h_voluta + h_pirulo)
print("Altezza perno-pirulo: ", h_voluta)