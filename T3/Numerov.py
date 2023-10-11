import numpy as np
import matplotlib.pyplot as plt
#b)

hbar=1
m=1
omega=1
# Despejando Rn=Vn-En tal que y''-yRn=S(x)=0
def Rn(En, Vn):
    return (Vn-En)

#c)
#-5<x<5 con n=1000
x_values=np.linspace(-5,5,num=1000)

#d)
V_values=[]
for i in x_values:
    V_values.append(0.5* i**2)
N=1000
    
def numerov(V_values, E):
    V_values=np.array(V_values)
    psi=np.zeros(N)
    psi[0] = 0.0
    psi[1] = 1e-5

    h=1/(N-1)
    h2=h**2
    Rn=200*(V_values-E)
    #Numerov
    for i in range(1, len(x_values)-1):
        psi[i+1]=(2*psi[i]*(1+(5/12)*h2*Rn[i])-psi[i-1]*(1-(1/12)*h2*Rn[i-1]))/(1-(1/12)*h2*Rn[i+1])
        
    return psi/max(psi)


def eigenval(V_values,E, N):
    nlist = [0,1,2,3,4,5]
    eig = []
    V_values = np.array(V_values)
    for item in nlist:
        E = item
        Rn = 200*(V_values-E)
        P1 = numerov(V_values,E)[-1]
        dE = 0.001
        a=True
        while a:
            Rn = 200*(V_values - E)
            P2 = numerov(V_values,E)[-1]
            if P1*P2<0:
                dE = -dE/2
                if (abs(dE)<1e-12):
                    eig.append(E)
                    a=False
            P1 = P2
            E = E + dE
    return eig

energy_levels= eigenval(V_values,0,N)
print(energy_levels)
#f)

for i,energy_level in enumerate(energy_levels):
    psi=(numerov(V_values,energy_level))
    print(f"Eigenvalue{i}: {energy_level}")
    rounded_energy=round(energy_level,1)
    plt.plot(x_values,psi,label=f'E_{i}: {rounded_energy}')

plt.legend()
plt.title("Potencial oscilador armonico cuantico")
plt.show()

#i)
def Geigenval(V_values,E, N):
    nlist =[0,-1,-2,-3,-4]#[-5.5,-6,-7,-8,-9]
    nlistd=[]
    for i in nlist:
        nlistd.append(i+E)
    eig = []
    V_values = np.array(V_values)
    for item in nlistd:
        E = item
        Rn = 200*(V_values-E)
        P1 = numerov(V_values,E)[-1]
        dE = 0.001
        a=True
        while a:
            Rn = 200*(V_values - E)
            P2 = numerov(V_values,E)[-1]
            if P1*P2<0:
                dE = -dE/2
                if (abs(dE)<1e-12):
                    eig.append(E)
                    a=False
            P1 = P2
            E = E - dE
    return eig


GV_values=[]
for i in x_values:
    GV_values.append(-10*np.exp(-i**2 /20))
    
Genergy_levels=Geigenval(GV_values,-5.5,N)
print("Potencial Gaussiano: ")
print(Genergy_levels)
for i,energy_level in enumerate(Genergy_levels):
    psi=(numerov(GV_values,energy_level))
    print(f"Eigenvalue{i}: {energy_level}")
    rounded_energy=round(energy_level,2)
    plt.plot(x_values,psi,label=f'E_{i}: {rounded_energy}')

plt.legend()
plt.title("Potencial Gaussiano")
plt.plot()

#j)
#x_values=np.linspace(-2,2,num=1000)
def Reigenval(V_values,E, N):
    nlist =[0,-1]#[-5.5,-6,-7,-8,-9]
    eig = []
    V_values = np.array(V_values)
    for item in nlist:
        E=item
        Rn=200*(V_values-E)
        P1=numerov(V_values,E)[-1]
        dE=0.001
        a=True
        while a:
            Rn=200*(V_values-E)
            P2=numerov(V_values,E)[-1]
            if P1*P2<0:
                dE=-dE/2
                if (abs(dE)<1e-12):
                    eig.append(E)
                    a=False
            P1=P2
            E=E-dE
    return eig


RV_values=[]
for i in x_values:
    RV_values.append(-4/((1+i**2)**2))
    
Renergy_levels=Reigenval(RV_values,0,N)
print("Potencial racional: ")
print(Renergy_levels)












