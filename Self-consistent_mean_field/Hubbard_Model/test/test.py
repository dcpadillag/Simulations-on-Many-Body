import numpy as np
from scipy.linalg import eig
import matplotlib.pyplot as plt
import sys

t=float(sys.argv[1]) #Hopping
U=float(sys.argv[2]) #Interaction

#Lattice Parameters
L=20#Lattice side length
N=L*L
Ne=float(sys.argv[3])
Ntotal=int(Ne*N)
Ncell=(L-1)*(L-1)

a=1.

Deltak=2*np.pi/L
M=int(2*np.pi/(a*Deltak))
print M

kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)
Kx, Ky = np.meshgrid(kx, ky)

Eup=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)
Edown=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)

Nup=Ntotal/2
Ndown=Ntotal-Nup
Eup=np.sort(Eup.reshape((1,M*M)))
Edown=np.sort(Edown.reshape((1,M*M)))
E=(Eup[:Nup].sum()+Edown[:Ndown].sum())/N
print E
a=1.

b1=np.array([1.,1.])
b2=np.array([1.,-1.])

kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)
Kx, Ky = np.meshgrid(kx, ky)
def gamma(k):
	return -(1+np.exp(-1j*np.dot(k,b1+b2))+np.exp(-1j*np.dot(k,b1))+np.exp(-1j*np.dot(k,b2)))

Eup1=np.zeros((M,M))
Edown1=np.zeros((M,M))
Eup2=np.zeros((M,M))
Edown2=np.zeros((M,M))
Eup=np.zeros((M,M))
for i in range(M):
	for j in range(M):
		k=np.array([kx[i],ky[j]])

		Hup=np.array([[0.,t*gamma(k)],[t*gamma(k).conjugate(),0.]])
		vup, wup=eig(Hup)
		Eup1[i][j]=(vup[0]).real
		Eup2[i][j]=(vup[1]).real
		vup=np.sort(vup.real)
		Hdown=np.array([[0.,t*gamma(k)],[t*gamma(k).conjugate(),0.]])
		vdown, wdown=eig(Hup)
		Edown1[i][j]=(vdown[0]).real
		Edown2[i][j]=(vdown[1]).real

Nup1=Nup/2
Nup2=Nup/2
Ndown1=Ndown/2
Ndown2=Ndown/2

Eup1=Eup1.reshape((1,M*M))
Eup2=Eup2.reshape((1,M*M))
Edown1=Edown1.reshape((1,M*M))
Edown2=Edown2.reshape((1,M*M))
#Eup1=np.sort(Eup1)
#Edown1=np.sort(Edown1)
#Eup2=np.sort(Eup2)
#Edown2=np.sort(Edown2)
Eup=np.hstack((Eup1,Eup2))
Eup=np.sort(Eup)
print Eup[0][:Nup].sum()/N
#Edown=np.hstack((Eup1,Eup2))








