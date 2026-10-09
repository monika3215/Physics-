import numpy as np
import matplotlib.pyplot as plt
# constants
hbar=1
m=1
w=1

def VPot(x):
  return 0.5*m*w**2*x**2

# grid points
N=100
# upper and lower limit
a=-5
b=5
x=np.linspace(a,b,N)
dx=x[1]-x[0]
# creating Kinetic energy matrix
T=np.zeros((N-2,N-2))
# diagonal elements
for i in range(N-2):
    for j in range(N-2):
        if i==j:
         T[i,j]=hbar**2/(m*dx**2) 
# lower diagonal
    if i>0: 
         T[i,i-1]=-hbar**2/(2*m*dx**2)
# N-3 upper diagonal
    if i<N-3:
         T[i,i+1]=-hbar**2/(2*m*dx**2)
print(T)

# creating potential energy matrix
V=np.zeros((N-2,N-2))
for i in range (N-2):
  for j in range (N-2):
      if i==j:
       V[i,j]=VPot(x[i+1])
      
print(V)
H=T+V
print(H)
E,Wavefunction=np.linalg.eigh(H)
print("Eigen Values:",E[:3])

# plotting ground state wave function
#  np.zeros create empty array for wave function
plt.figure(figsize=(8, 5))

for n in range(3):
    psi=np.zeros(N)
    psi[1:-1]=Wavefunction[:,n]
    if psi[np.argmax(np.abs(psi))]<0:
        psi=-psi
        # psi=psi/normalisation factor
    psi=psi/np.sqrt(np.sum(np.abs(psi)**2))*dx
    plt.plot(x,psi,label=f"n={n}, E={E[n]:.3f}")


plt.xlabel("x")
plt.ylabel("psi")
plt.title("First Three Harmonic Oscillator Wavefunctions ")
plt.legend()
plt.grid()
plt.show()