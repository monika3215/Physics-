import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
hbar=1
u=1
e=1
k=1

# 1. Computational Domain and Grid
N = 1000  # Number of grid points
a = 0.01
b = 20 # Maximum radial distance (in Bohr radii)
r = np.linspace( a,b,N)  # Avoid r=0 to prevent division by zero
dr = r[1] - r[0]  # Step size

# 2. Quantum Number and Potential
l=0
def Vpot():
  return (l * (l + 1))*hbar*hbar / (2*u*r**2) - e*e*k /r  # Effective Coulomb potential

# 3. Kinetic Energy Matrix (Finite Difference Second Derivative)
T=np.zeros((N-2,N-2))
# diagonal elements
for i in range(N-2):
    for j in range(N-2):
        if i==j:
         T[i,j]=hbar**2/(u*dr**2) 
# lower diagonal
    if i>0: 
         j=i-1
         T[i,j]=-hbar**2/(2*u*dr**2 )
# N-3 upper diagonal
    if i<N-3:
         j=i+1
         T[i,j]=-hbar**2/(2*u*dr**2)
print(T)

V=np.zeros((N-2,N-2))
for i in range (N-2):
  for j in range (N-2):
      if i==j:
        V[i,j]=Vpot()[i+1]
print(V)

H= T+ V

# Fix boundary conditions: u(0) = 0 and u(r_max) = 0 (handled by restricted interior or zero padding)

# 4. Solve the Eigenvalue Problem
E,Radial_Wavefunction = eigh(H)
print("Lowest three energies:")
print(E[0:3])

# 5. Ground-state radial wavefunction
u0 = Radial_Wavefunction[:, 2]
# Convert discrete normalization to continuum normalization
u0 /= np.sqrt(np.sum(np.abs(u0)**2) * dr)
# Add boundary conditions
u_full = np.zeros(N)
u_full[1:-1] = u0

# Plot
plt.figure(figsize=(8, 5))
plt.plot(r, u_full, label=f"$E_0 = {E[0]:.6f}$ Ha")

plt.title("Hydrogen Atom Ground-State Radial Wavefunction")
plt.xlabel("$r$ (Bohr radii)")
plt.ylabel("$u(r)$")
plt.legend()
plt.grid(True)
plt.show()