import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
import sys

def gamma(k):
	return -(1+np.exp(-1j*np.dot(k,b1+b2))+np.exp(-1j*np.dot(k,b1))+np.exp(-1j*np.dot(k,b2)))

t=float(sys.argv[1]) #Hopping
#t=0.0+((0.76/38.)*int(sys.argv[1]))
U=float(sys.argv[2]) #Interaction

#Lattice Parameters
L=int(sys.argv[3])#Lattice side length
N=L*L
#print N
Ne=float(sys.argv[4])
#Ne=0.0+(0.1*int(sys.argv[4]))
Ntotal=int(Ne*N) #Total Number of particles

a=1.#Red parameter
a1x=a
a1y=0.
a2x=0.5*a
a2y=0.5*np.sqrt(3)*a


Deltak=2*np.pi/L #Delta of k
M=int(2*np.pi/(a*Deltak)) 

#Momentum
kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)
Kx, Ky = np.meshgrid(kx, ky)


EF=np.array([]) #Array for the total energy
mF=np.array([]) #Array for magnetization
Em={} #Dictionary for minimal state
densityup={} #Dictionary for densities of each state

for Nup in range(Ntotal+1): #Run over the total number of Nup
	Ndown=Ntotal-Nup #Ndown particles

	nup=float(Nup)/float(N) #up spin density 
	ndown=float(Ndown)/float(N) #down spin density 

	mF=np.append(mF,(nup-ndown)/(nup+ndown)) #magnetization of the state

	Eup=-2.*t*np.cos(Kx*a1x+Ky*a1y)-2.*t*np.cos(Kx*a2x+Ky*a2y)-2.*t*np.cos(Kx*(a2x-a1x)+Ky*(a2y-a1y))+U*ndown #Enegy level up (mean field)
	Edown=-2.*t*np.cos(Kx*a1x+Ky*a1y)-2.*t*np.cos(Kx*a2x+Ky*a2y)-2.*t*np.cos(Kx*(a2x-a1x)+Ky*(a2y-a1y))+U*nup #Enegy level down (mean field)

	Eup=np.sort(Eup.reshape([1,M*M])[0]) #Sort the energy leves up
	Edown=np.sort(Edown.reshape([1,M*M])[0]) #Sort the energy leves down

	EF=np.append(EF,((Eup[:Nup].sum()+Edown[:Ndown].sum())/(N))-U*nup*ndown) # total enegy of state

	Em.update({EF[Nup]:mF[Nup]})
	densityup.update({mF[Nup]:nup})

m_min=Em[(np.sort(EF))[0]] #Magnetization of ground state
nup=densityup[m_min] #Up density of ground state

if abs(m_min)==1.:
#	print "Ferromagnetic phase"
	print t, Ne, "F"
#	density=np.array([0.5*nup, 0.5*(Ne-nup),0.5*nup, 0.5*(Ne-nup)])
#	np.save("density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+"_L_"+str(L),density)
	density=np.array([nup, (Ne-nup)])
	np.save("triangular/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne),density)
	plt.plot(mF,EF,'black')
	plt.xlabel("$m$",fontsize=18)
	plt.ylabel("$E$",fontsize=18)
	plt.savefig("triangular/energy_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")
	plt.close()
#	plt.show()
else:
#	print "Paramagnetic phase"
	print t, Ne, "P"
	density=np.array([nup, (Ne-nup)])
	np.save("triangular/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne),density)
	plt.plot(mF,EF,'black')
	plt.xlabel("$m$",fontsize=18)
	plt.ylabel("$E$",fontsize=18)
	plt.savefig("triangular/energy_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")	
	plt.close()
#	plt.show()







