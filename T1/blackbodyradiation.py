import numpy as np
import matplotlib.pyplot as plt

h = 6.626e-34
c = 3.0e+8
k_b = 1.38e-23
def energy_dist(lam,T):
    a = 2*h*c**2
    exponent = h*c/(lam*k_b*T)
    intensity = a/((lam**5)*(np.exp(exponent)-1.0))
    return intensity

wavelengths = np.arange(1e-9, 3e-6, 1e-9) 

intensity3000 = energy_dist(wavelengths, 3000.)
intensity5000 = energy_dist(wavelengths, 5000.)
intensity8000 = energy_dist(wavelengths, 8000.)
intensity10000 = energy_dist(wavelengths, 10000.)


plt.plot(wavelengths*1e9,intensity3000,label = '$T=3000K$',color ='red') 
# plot intensity4000 versus wavelength in nm as a red line
plt.plot(wavelengths*1e9, intensity5000,label = '$T=5000K$',color='#0ddeda') # 5000K green line
plt.plot(wavelengths*1e9, intensity8000,label = '$T=8000K$',color='#ed980e') # 6000K blue line
plt.plot(wavelengths*1e9, intensity10000,label = '$T=10000K$', color='#a832a6') # 7000K black line
plt.legend()
plt.xlabel('longitud de onda (nm)')
plt.ylabel('densidad de energia espectral($W*sr^{-1}m^{3}$)')
# show the plot
plt.show()
