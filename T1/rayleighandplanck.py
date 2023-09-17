import numpy as np
import matplotlib.pyplot as plt

h = 6.626e-34
c = 3.0e+8
k_b = 1.38e-23
def planck_f(freq, T):
    a = 2.0*h*(freq**3)
    b = h*freq/(k_b*T)
    intensity =  a/( (c**2 * (np.exp(b) - 1.0) ))
    return intensity
def rayleigh(freq,T):
    a = 2*k_b*T*freq**2
    b = c**2
    return a/b

wavelengths = np.arange(1e-9,3e-6,1e-9) 
freq = np.arange(1e13,3e17,1e13) 

intensity10000 = planck_f(freq, 10000)
rayleigh10000 = rayleigh(freq, 10000)
difference = (rayleigh10000[0:250]-intensity10000[0:250])*1e7
Meansquarederror = np.power(difference,2)
Meansquarederror = np.sum(Meansquarederror)
Meansquarederror = (1/len(rayleigh10000))* Meansquarederror
print(Meansquarederror)
plt.plot(freq*1e-12, rayleigh10000,label = '$Rayleigh$', color='#a832a6') # 7000K black line
plt.plot(freq*1e-12, intensity10000,label = '$Planck$', color='blue') # 7000K black line
plt.legend()
plt.xlabel('frecuencia (THz)')
plt.ylabel('densidad de energia espectral($W*sr^{-1}m^{3}$)')
plt.ylim(0,0.2e-6)
plt.xlim(0,2500)
# show the plot
plt.show()
