import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
import sys
import itertools

def gamma(k):
	return -(1+np.exp(-1j*np.dot(k,b1+b2))+np.exp(-1j*np.dot(k,b1))+np.exp(-1j*np.dot(k,b2)))

t=float(sys.argv[1]) #Hopping
U=float(sys.argv[2]) #Interaction

#Lattice Parameters
L=int(sys.argv[3])#Lattice side length
N=L*L
#print N
Ne=float(sys.argv[4]) #Number of particles per site
Ntotal=int(Ne*N) #Total Number of particles

a=1.#Red parameter

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

	Eup=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)+U*ndown #Enegy level up (mean field)
	Edown=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)+U*nup #Enegy level down (mean field)
          
	Eup=np.sort(Eup.reshape([1,M*M])[0]) #Sort the energy leves up
	Edown=np.sort(Edown.reshape([1,M*M])[0]) #Sort the energy leves down

	EF=np.append(EF,((Eup[:Nup].sum()+Edown[:Ndown].sum())/(N))-U*nup*ndown) # total enegy of state

	Em.update({EF[Nup]:mF[Nup]})
	densityup.update({mF[Nup]:nup})

m_min=Em[(np.sort(EF))[0]] #Magnetization of ground state
nup=densityup[m_min] #Up density of ground state

if abs(m_min)==1.:
	print "Ferromagnetic phase"
#	density=np.array([0.5*nup, 0.5*(Ne-nup),0.5*nup, 0.5*(Ne-nup)])
#	np.save("density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+"_L_"+str(L),density)
	density=np.array([nup, (Ne-nup),nup, (Ne-nup)])
	np.save("square/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne),density)
	plt.plot(mF,EF,'black')
	plt.xlabel("$m$",fontsize=18)
	plt.ylabel("$E$",fontsize=18)
	plt.savefig("square/energy_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")
	plt.show()
else:
	b1=np.array([1.,1.])
	b2=np.array([1.,-1.])


	E=np.array([]) #Array for the total energy
	m=np.linspace(-1.,1.,Ntotal) #Array for magnetization
	Em={} #Dictionary for minimal state
	n=Ne
	for l in range(Ntotal):
		nup1=0.5*n*(1.+m[l]) #up spin density of site 1
		ndown1=0.5*n*(1.-m[l]) #down spin density of site 1

		nup2=0.5*n*(1.-m[l]) #up spin density of site 2
		ndown2=0.5*n*(1.+m[l]) #down spin density of site 2


		#Total number of particles for configuration of magnetization m
		Nup1=int((nup1)*N)
		Nup2=int((nup2)*N)
		Ndown1=int((ndown1)*N)
		Ndown2=int((ndown2)*N)
		Nup=Nup1+Nup2
		Ndown=Ndown1+Ndown2

		#Energy leves
		Eup1=np.zeros([M,M])
		Edown1=np.zeros([M,M])
		Eup2=np.zeros([M,M])
		Edown2=np.zeros([M,M])
		for i,j in itertools.product(range(M),range(M)):
			k=np.array([kx[i],ky[j]])

			Hup=np.array([[U*ndown1,t*gamma(k)],[t*gamma(k).conjugate(),U*ndown2]])
			vup, wup=eig(Hup)
			Eup1[i][j]=(vup[0]).real
			Eup2[i][j]=(vup[1]).real

			Hdown=np.array([[U*nup1,t*gamma(k)],[t*gamma(k).conjugate(),U*nup2]])
			vdown, wdown=eig(Hdown)
			Edown1[i][j]=(vdown[0]).real
			Edown2[i][j]=(vdown[1]).real

		#Sort the energy leves
		Eup1=Eup1.reshape((1,M*M))
		Eup2=Eup2.reshape((1,M*M))
		Eup=np.hstack((Eup1,Eup2))
		Eup=np.sort(Eup)
		Edown1=Edown1.reshape((1,M*M))
		Edown2=Edown2.reshape((1,M*M))
		Edown=np.hstack((Edown1,Edown2))
		Edown=np.sort(Edown)

		#Total energy of the system
		E=np.append(E,((Eup[0][:Nup].sum()+Edown[0][:Ndown].sum())/(N))-U*(nup1*ndown1+nup2*ndown2))
		Em.update({E[l]:m[l]})

	m_min=Em[(np.sort(E))[0]]#Magnetization of ground state

	if m[Ntotal/2]<=abs(m_min)<=m[1+Ntotal/2]:
		print "Paramagnetic phase"
		density=np.array([nup1,ndown1,nup2,ndown2])
		np.save("square/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne),density)
		plt.plot(mF,EF,'black')
		plt.xlabel("$m$",fontsize=18)
		plt.ylabel("$E$",fontsize=18)
		plt.savefig("square/energy_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")
		plt.show()

	else:
		density=np.array([nup1,ndown1,nup2,ndown2])
#		np.save("density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+"_L_"+str(L),density)
		np.save("square/density_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne),density)
		print "Antiferromagnetic phase"
		plt.plot(m,E,'black')
		plt.xlabel("$m$",fontsize=18)
		plt.ylabel("$E$",fontsize=18)
		plt.savefig("square/energy_t_"+str(t)+"_U_"+str(U)+"_Ne_"+str(Ne)+".png")
		plt.show()







