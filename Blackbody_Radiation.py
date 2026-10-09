import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c,h,k

def planck_wavelength(lemda,T):
    a=2*h*c**2
    b=h*c/(lemda*k*T)
    B=a/(lemda**5*(np.exp(b)-1))
    return B 

N=1000
x=1e-7 
y=3e-6
wavelengths=np.linspace(x,y,N)  
temperatures=[3000,4000,5000]

plt.figure(figsize=(8,5))

for T in temperatures:
    B=planck_wavelength(wavelengths,T)
    plt.plot(wavelengths*1e9,B,label=f'T={T}K')

plt.title('BLACKBODY RADIATION CURVE')
plt.xlabel('Wavelength(nm)')   
plt.ylabel('Spectral Radiance [W/(m**2.sr.m)]')
plt.legend()
plt.grid(True)
plt.show()