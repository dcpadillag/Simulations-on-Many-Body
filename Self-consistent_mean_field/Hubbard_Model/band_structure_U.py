import numpy as np
from scipy.linalg import eig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

t=float(sys.argv[1]) #Hopping
U=float(sys.argv[2]) #Interaction

#Lattice Parameters
L=int(sys.argv[3])
N=L*L
Ncell=2*(N/4)
Ne=float(sys.argv[4])
a=1.
b=2*np.pi/a
a1=a*np.array([1,0])
a2=a*np.array([0,1])

b1=np.array([1.,1.])
b2=np.array([1.,-1.])
Deltak=2*np.pi/L
M=int(2*np.pi/(a*Deltak))
print M

#State densities
#density=np.load("density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+"_L_"+str(L)+".npy")
density=np.load("square/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".npy")
Nup1=density[0]
Ndown1=density[1]
Nup2=density[2]
Ndown2=density[3]

#Nup1=0.62
#Ndown1=0.18
#Nup2=0.18
#Ndown2=0.62

print Nup1,Ndown1,Nup2,Ndown2

kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)

#def gamma(k):
#    return -(2.*np.cos(np.dot(k,a1))+2.*np.cos(np.dot(k,a2)))

def gamma(k):
	return -(1+np.exp(-1j*np.dot(k,b1+b2))+np.exp(-1j*np.dot(k,b1))+np.exp(-1j*np.dot(k,b2)))

def diag_up(k,Nup1,Nup2,Ndown1,Ndown2):
    Hup=np.array([[U*Ndown1,t*gamma(k)],[t*gamma(k).conjugate(), U*Ndown2]])
    Eup , Psiup= eig(Hup)
    return Eup, Psiup

def diag_down(k,Nup1,Nup2,Ndown1,Ndown2):
    Hdown=np.array([[U*Nup1,t*gamma(k)],[t*gamma(k).conjugate(), U*Nup2]])
    Edown , Psidown= eig(Hdown)
    return Edown, Psidown

def State(kx,ky,Nup1,Nup2,Ndown1,Ndown2):
	M=len(kx)
	Eup1=np.zeros([M,M])
	Eup2=np.zeros([M,M])
	Edown1=np.zeros([M,M])
	Edown2=np.zeros([M,M])
	
	EU=-U*(Nup1*Ndown1+Nup2*Ndown2)
	for i in range(M):
		for j in range(M):
			k=np.array([kx[i],ky[j]])
			vup, wup = diag_up(k,Nup1,Nup2,Ndown1,Ndown2)
			vdown, wdown = diag_down(k,Nup1,Nup2,Ndown1,Ndown2)
			Eup1[i][j]=(vup[0]).real
			Eup2[i][j]=(vup[1]).real
			Edown1[i][j]=(vdown[0]).real
			Edown2[i][j]=(vdown[1]).real
	return Eup1, Eup2, Edown1, Edown2
#	return Eup1+EU, Eup2+EU, Edown1+EU, Edown2+EU


def DOS(E,E01,E02,dE):
    sigma=1.5*dE
    suma=0.
    delta =((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-E01)/sigma)**2)).sum()
    delta+=((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-E02)/sigma)**2)).sum()
    suma=suma+delta
    return suma

#########################################################################

#Diagonalization
Eup1, Eup2, Edown1, Edown2 = State(kx,ky,Nup1,Nup2,Ndown1,Ndown2)

E=np.linspace(-5.*t,15.*t,M)
dE=(max(E)-min(E))/M
Dup=np.array([])
Ddown=np.array([])
for i in range(len(E)):
     Dup=np.append(Dup,DOS(E[i],Eup1,Eup2,dE))
     Ddown=np.append(Ddown,DOS(E[i],Edown1,Edown2,dE))

Dup=Dup/(4*N)
Ddown=Ddown/(4*N)

plt.plot(Dup,E/t,color='red',label='$Espin$ $1/2$')
plt.plot(Ddown,E/t,"--",color='blue',label='$Espin$ $-1/2$')

plt.ylim([-5.,16.0])
plt.legend(loc='best')
plt.xlabel("$DOS$ $por$ $celda$ $unidad$",fontsize=18)
plt.ylabel("$E/t$",fontsize=18)
plt.savefig("square/DOS_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")	
plt.show()
