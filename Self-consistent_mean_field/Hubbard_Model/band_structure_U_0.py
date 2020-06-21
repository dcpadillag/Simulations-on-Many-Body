import numpy as np
from scipy.linalg import eig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

t=1. #Hopping
beta=300.

#Lattice Parameters
L=100
Ntotal=L*L
Ncell=2*(Ntotal/4)
#Ne=Ntotal/Ncell
Ne=1.6
a=1.
b=2*np.pi/a
a1=a*np.array([1,0])
a2=a*np.array([0,1])

b1=np.array([1.,1.])
b2=np.array([1.,-1.])

Deltak=2*np.pi/L
M=int(2*np.pi/(a*Deltak))

kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)

def gamma(k):
	return -(1+np.exp(-1j*np.dot(k,b1+b2))+np.exp(-1j*np.dot(k,b1))+np.exp(-1j*np.dot(k,b2)))

def E0(k):
    return 2.*np.cos(np.dot(k,a1))+2.*np.cos(np.dot(k,a2))

def DOS(E,E1):
    sigma=1.5*dE
    suma=0.
    delta =((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-E1)/sigma)**2)).sum()
#    delta+=((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-E2)/sigma)**2)).sum()
    suma=suma+delta
    return suma


def F(mu,D,E):
    integral=0.
    dE=(max(E)-min(E))/len(E)
    for i in range(len(E)):
        integral+=(D[i]/(np.exp(beta*(E[i]-mu))+1))*dE
    return integral-Ne

def biseccion(mu1,mu2,D,E):
    F1=F(mu1,D,E)
    F2=F(mu2,D,E)
    if (F2*F1<0):
        mu_m=0.5*(mu1+mu2)
        Fb=F(mu_m,D,E)
        while(abs(Fb)>10.E-6):
            if (Fb>0):
                mu2=mu_m
                mu_m=0.5*(mu1+mu2)
            else:
                mu1=mu_m
                mu_m=0.5*(mu1+mu2)
            Fb=F(mu_m,D,E)
    else:
        print ("Wrong chemical potential interval")
	mu_m=0

    return mu_m

#########################################################################

E1=np.zeros([M,M])
E2=np.zeros([M,M])
for i in range(M):
	for j in  range(M):
		k=np.array([kx[i],ky[j]])
		E1[i][j]=-t*E0(k)
#		E1[i][j]=-t*np.sqrt(((gamma(k)*gamma(k).conjugate()).real))
#		E2[i][j]=t*np.sqrt(((gamma(k)*gamma(k).conjugate()).real))

cp = plt.contourf(kx, ky, E1,cmap=plt.cm.bone)
plt.colorbar(cp)
plt.xlabel("$k_{x}$",fontsize=18)
plt.ylabel("$k_{y}$",fontsize=18)
plt.savefig("square/brioullin_square.png")	
plt.show()
#Density of State
E=np.linspace(-5.*t,5.*t,M)
dE=(max(E)-min(E))/M
D=np.array([])
for i in range(len(E)):
    D=np.append(D,DOS(E[i],E1))

D=D/Ncell

plt.plot(D,E/t)
plt.show()

