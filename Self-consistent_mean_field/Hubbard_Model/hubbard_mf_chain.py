import numpy as np
import matplotlib.pyplot as plt
import sys

t=1. #Hopping
U=1.0 #Interaction

#Lattice Parameters
L=256#Lattice side length
Ntotal=L/2
#Ncell=2*(Ntotal/4)
#Ne=Ntotal/Ncell
a=1.

Deltak=2*np.pi/L
M=int(2*np.pi/(a*Deltak))
print M

E=np.array([])
for Nup in range(Ntotal+1):
	Ndown=Ntotal-Nup
	nup=float(Nup)/float(L)
	ndown=float(Ndown)/float(L)
	kx=np.linspace(-np.pi/a,np.pi/a,M)
	Eup=-2.*t*np.cos(kx*a)+U*ndown
	Eup=np.sort(Eup)
	Edown=-2.*t*np.cos(kx*a)+U*nup
	Edown=np.sort(Edown)
	E=np.append(E,((Eup[:Nup].sum()+Edown[:Ndown].sum())/L)-U*nup*ndown)

plt.plot(E)
plt.show()
