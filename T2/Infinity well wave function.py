import numpy as np
import matplotlib.pyplot as plt

hb=1.054*10**(-34)

def psi(x,t,a,n,m):
    
    Cn=8*15**(1/2)/(n*np.pi)**3
    ogp=(2/a)**(1/2)*np.sin(n*np.pi*x/a)
    ex=np.exp(-1j* n**2 *hb*np.pi**2 *t / (2*m* a**2))
    
    psi=Cn*ogp*ex
    
    if n % 2==0: 
        return 0
    else:
        return psi

def probability_distribution(psi):
    psi_conjugate = np.conjugate(psi)
    return np.real(psi*psi_conjugate)

#ps=(psi(5,2,1,3,1))
#print(ps)
#print(probability_distribution(ps))

# Constants
a = 1
t=5
t_values = np.linspace(0, t, 500) 
x_values = np.linspace(0, a, 500) 
m = 1.0


n_values_str = input("Enter values of n (comma-separated): ")
n_values = [int(n) for n in n_values_str.split(",")]

for n in n_values:
    psi_values = np.zeros((len(t_values),len(x_values)), dtype=complex)
    prob_dist_values = np.zeros((len(t_values),len(x_values)))

    for i, t in enumerate(t_values):
        psi_values[i,:] =psi(x_values,t,a,n,m)
        prob_dist_values[i,:]=probability_distribution(psi_values[i, :])

    plt.figure(figsize=(12, 6))
    plt.subplot(121)
    plt.title(f'Wave Function (n={n})')
    plt.plot(x_values, np.real(psi_values[-1, :]), label='Real')
    #plt.plot(x_values, np.imag(psi_values[-1, :]), label='Imaginary')
    plt.xlabel('Position (x)')
    plt.ylabel('Wave Function')
    plt.legend()
    plt.subplot(122)
    plt.title(f'Probability Distribution (n={n})')
    plt.plot(x_values, prob_dist_values[-1, :])
    plt.xlabel('Position (x)')
    plt.ylabel('Probability')
    plt.tight_layout()
    plt.show()
