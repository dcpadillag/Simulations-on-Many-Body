import numpy as np
from scipy.linalg import eig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

t=float(sys.argv[1]) #Hopping
#t=0.0+(0.76/38.)*int(sys.argv[1])
#print t
U=float(sys.argv[2]) #Interaction

#Lattice Parameters
L=int(sys.argv[3])
N=L*L
Ne=float(sys.argv[4])
#Ne=0.0+0.1*int(sys.argv[2])
#print Ne
Ntotal=int(Ne*N) #Total Number of particles

a=1.#Red parameter
a1x=a
a1y=0.
a2x=0.5*a
a2y=0.5*np.sqrt(3)*a

Deltak=2*np.pi/L
M=int(2*np.pi/(a*Deltak))
print M

#State densities
#density=np.load("density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+"_L_"+str(L)+".npy")
density=np.load("triangular/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".npy")
Nup=density[0]
Ndown=density[1]


kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)


def DOS(E,E0,dE):
    sigma=1.5*dE
    delta =((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-E0)/sigma)**2)).sum()
    return delta


#########################################################################

#Diagonalization
Kx, Ky = np.meshgrid(kx, ky)
Eup=-2.*t*np.cos(Kx*a1x+Ky*a1y)-2.*t*np.cos(Kx*a2x+Ky*a2y)-2.*t*np.cos(Kx*(a2x-a1x)+Ky*(a2y-a1y))+U*Ndown #Enegy level up (mean field)
Edown=-2.*t*np.cos(Kx*a1x+Ky*a1y)-2.*t*np.cos(Kx*a2x+Ky*a2y)-2.*t*np.cos(Kx*(a2x-a1x)+Ky*(a2y-a1y))+U*Nup #Enegy level down (mean field)

if U==0:
	cp = plt.contourf(Kx, Ky, Eup,cmap=plt.cm.bone)
	plt.colorbar(cp)
	plt.xlabel("$k_{x}$",fontsize=18)
	plt.ylabel("$k_{y}$",fontsize=18)
	plt.savefig("triangular/brioullin_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")	
	plt.close()
#	plt.show()


E=np.linspace(-7.*t,20.*t,M)
dE=(max(E)-min(E))/M
Dup=np.array([])
Ddown=np.array([])
for i in range(len(E)):
     Dup=np.append(Dup,DOS(E[i],Eup,dE))
     Ddown=np.append(Ddown,DOS(E[i],Edown,dE))

Dup=Dup/(N)
Ddown=Ddown/(N)

plt.plot(Dup,E/t,'red',label='$Espin$ $1/2$')
plt.plot(Ddown,E/t,'blue',label='$Espin$ $-1/2$')
plt.ylim([-8.,22.0])
plt.legend(loc='best')
plt.xlabel("$DOS$ $por$ $celda$ $unidad$",fontsize=18)
plt.ylabel("$E/t$",fontsize=18)
plt.savefig("triangular/DOS_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")	
plt.show()
