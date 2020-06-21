import numpy as np
from scipy.linalg import eig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys



t=0.077 #Hopping
U=1. #Interaction
beta=333. #1/kbT

#Lattice Parameters
L=10#Lattice side length
Ntotal=L*L #total number of sites
Ncell=2*(Ntotal/4)
Ne=Ntotal/Ncell
a=1.
b=2*np.pi/a
a1=a*np.array([1,0])
a2=a*np.array([0,1])

b1=np.array([1.,1.])
b2=np.array([1.,-1.])

Deltak=2*np.pi/L
M=int(2*np.pi/(a*Deltak))
print M

#Initial State
Ne=1.6 #total number of particles per cell
d=1.*Ne/2
Nup1=Ne/4+d
Ndown1=Ne/4-d
Nup2=Ne/4+d
Ndown2=Ne/4-d

#mu=0.419049263
#mu=0.0465594100952
mu=0.5*U*Ne
#Nup1=(Ne/2.)+d
#Ndown1=0.
#Nup2=(Ne/2.)-d
#Ndown2=0.

kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)

#def gamma(k):
#    return -(2.*np.cos(np.dot(k,a1))+2.*np.cos(np.dot(k,a2)))

def gamma(k):
	return -(1+np.exp(-1j*np.dot(k,b1+b2))+np.exp(-1j*np.dot(k,b1))+np.exp(-1j*np.dot(k,b2)))

def Eup(k,Nup1,Nup2,Ndown1,Ndown2,mu):
    EU=-U*(Nup1*Ndown1+Nup2*Ndown2)
    Eupp=0.5*U*(Ndown1+Ndown2)+EU+np.sqrt((U*(Ndown1-Ndown2))**2+4*gamma(k)**2)
    Eupm=0.5*U*(Ndown1+Ndown2)+EU-np.sqrt((U*(Ndown1-Ndown2))**2+4*gamma(k)**2)
    return Eupp, Eupm

def Edown(k,Nup1,Nup2,Ndown1,Ndown2,mu):
    EU=-U*(Nup1*Ndown1+Nup2*Ndown2)
    Edownp=0.5*U*(Nup1+Nup2)+EU+np.sqrt((U*(Nup1-Nup2))**2+4*gamma(k)**2)
    Edownm=0.5*U*(Nup1+Nup2)+EU-np.sqrt((U*(Nup1-Nup2))**2+4*gamma(k)**2)
    return Edownp, Edownm

def diag_up(k,Nup1,Nup2,Ndown1,Ndown2,mu):
#    Hup=np.array([[U*Ndown1-mu,t*gamma(k)],[t*gamma(k).conjugate(), U*Ndown2-mu]])
    Hup=np.array([[U*Ndown1,t*gamma(k)],[t*gamma(k).conjugate(), U*Ndown2]])
    Eup , Psiup= eig(Hup)
    return Eup, Psiup

def diag_down(k,Nup1,Nup2,Ndown1,Ndown2,mu):
#    Hdown=np.array([[U*Nup1-mu,t*gamma(k)],[t*gamma(k).conjugate(), U*Nup2-mu]])
    Hdown=np.array([[U*Nup1,t*gamma(k)],[t*gamma(k).conjugate(), U*Nup2]])
    Edown , Psidown= eig(Hdown)
    return Edown, Psidown

def State(kx,ky,Nup1,Nup2,Ndown1,Ndown2,mu):
	M=len(kx)
	Eup1=np.zeros([M,M])
	Eup2=np.zeros([M,M])
	Edown1=np.zeros([M,M])
	Edown2=np.zeros([M,M])
	
	Aup1=[np.zeros([M,M]),np.zeros([M,M])]
	Aup2=[np.zeros([M,M]),np.zeros([M,M])]
	Adown1=[np.zeros([M,M]),np.zeros([M,M])]
	Adown2=[np.zeros([M,M]),np.zeros([M,M])]
	EU=-U*(Nup1*Ndown1+Nup2*Ndown2)
	for i in range(M):
		for j in range(M):
			k=np.array([kx[i],ky[j]])
			wup, vup = diag_up(k,Nup1,Nup2,Ndown1,Ndown2,mu)
			wdown, vdown = diag_down(k,Nup1,Nup2,Ndown1,Ndown2,mu)
			vup=vup.transpose()
			vdown=vdown.transpose()
			Eup1[i][j]=wup[0].real
			Eup2[i][j]=wup[1].real
			Edown1[i][j]=wdown[0].real
			Edown2[i][j]=wdown[1].real
			aux=vup[0]*vup[0].conjugate()
			Aup1[0][i][j]=aux[0].real
			Aup1[1][i][j]=aux[1].real
			aux=vup[1]*vup[1].conjugate()
			Aup2[0][i][j]=aux[0].real
			Aup2[1][i][j]=aux[1].real
			aux=vdown[0]*vdown[0].conjugate()
			Adown1[0][i][j]=aux[0].real
			Adown1[1][i][j]=aux[1].real
			aux=vdown[1]*vdown[1].conjugate()
			Adown2[0][i][j]=aux[0].real
			Adown2[1][i][j]=aux[1].real
	return Eup1+EU-mu, Eup2+EU-mu, Edown1+EU-mu, Edown2+EU-mu, Aup1, Aup2, Adown1, Adown2
#	return Eup1+Ne*EU-mu, Eup2+Ne*EU-mu, Edown1+Ne*EU-mu, Edown2+Ne*EU-mu, Aup1, Aup2, Adown1, Adown2


def DOS(E,Eup1,Eup2,Edown1,Edown2,dE):
    sigma=1.5*dE
    suma=0.
    delta =((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-Eup1)/sigma)**2)).sum()
    delta+=((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-Eup2)/sigma)**2)).sum()
    delta+=((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-Edown1)/sigma)**2)).sum()
    delta+=((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-Edown2)/sigma)**2)).sum()
    suma=suma+delta
    return suma

def F(mu,D,E,Ne):
    integral=0.
    dE=(max(E)-min(E))/len(E)
    for i in range(len(E)):
        integral+=(D[i]/(np.exp(beta*(E[i]-mu))+1))*dE
    return integral-Ne
#    return integral-Ntotal

def biseccion(mu1,mu2,Nup1,Nup2,Ndown1,Ndown2):
    #Diagonalization
    Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2 = State(kx,ky,Nup1,Nup2,Ndown1,Ndown2,mu1)
    #Density of State
    E=np.linspace(-10.*t,20.*t,M)
    dE=(max(E)-min(E))/M
    D1=np.array([])
    for i in range(len(E)):
        D1=np.append(D1,DOS(E[i],Eup1,Eup2,Edown1,Edown2,dE))

    D1=D1/(Ncell)

    Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2 = State(kx,ky,Nup1,Nup2,Ndown1,Ndown2,mu2)
    D2=np.array([])
    for i in range(len(E)):
        D2=np.append(D2,DOS(E[i],Eup1,Eup2,Edown1,Edown2,dE))

    D2=D2/(Ncell)

    Ne=Nup1+Nup2+Ndown1+Ndown2
    F1=F(mu1,D1,E,Ne)
    F2=F(mu2,D2,E,Ne)
    if (F2*F1<0):
        mu_m=0.5*(mu1+mu2)

	Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2 = State(kx,ky,Nup1,Nup2,Ndown1,Ndown2,mu_m)
	Dm=np.array([])
	for i in range(len(E)):
		Dm=np.append(Dm,DOS(E[i],Eup1,Eup2,Edown1,Edown2,dE))
        Dm=Dm/(Ncell)
        Fb=F(mu_m,Dm,E,Ne)

        while(abs(Fb)>10.E-6):
            if (Fb>0):
                mu2=mu_m
                mu_m=0.5*(mu1+mu2)
            else:
                mu1=mu_m
                mu_m=0.5*(mu1+mu2)

	    Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2 = State(kx,ky,Nup1,Nup2,Ndown1,Ndown2,mu_m)
	    Dm=np.array([])
	    for i in range(len(E)):
		Dm=np.append(Dm,DOS(E[i],Eup1,Eup2,Edown1,Edown2,dE))
            Dm=Dm/(Ncell)
            Fb=F(mu_m,Dm,E,Ne)
    else:
        print ("NO PANA")

    return mu_m, Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2

def occupation(mu,Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2):
    suma1=0.
    suma2=0.
    suma3=0.
    suma4=0.

    suma1+=((Aup1[0])/(np.exp(beta*(Eup1-mu))+1)).sum()+((Aup1[1])/(np.exp(beta*(Eup2-mu))+1)).sum()
    suma2+=((Aup2[0])/(np.exp(beta*(Eup1-mu))+1)).sum()+((Aup2[1])/(np.exp(beta*(Eup2-mu))+1)).sum()
    suma3+=((Adown1[0])/(np.exp(beta*(Edown1-mu))+1)).sum()+((Adown1[1])/(np.exp(beta*(Edown2-mu))+1)).sum()
    suma4+=((Adown2[0])/(np.exp(beta*(Edown1-mu))+1)).sum()+((Adown2[1])/(np.exp(beta*(Edown2-mu))+1)).sum()
    return suma1/Ncell, suma2/Ncell, suma3/Ncell, suma4/Ncell
#    return suma1/Ntotal, suma2/Ntotal, suma3/Ntotal, suma4/Ntotal
#    return suma1, suma2, suma3, suma4
#########################################################################

for i in range(10):
	#Chemical Potential
	mu, Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2=biseccion(-10.*t,20.*t,Nup1,Nup2,Ndown1,Ndown2)
#	N=F(mu,D,E,Ne)+Ne
#	print mu, N

	#New Occupation
	Nup1, Nup2, Ndown1, Ndown2 = occupation(mu,Eup1, Eup2, Edown1, Edown2, Aup1, Aup2, Adown1, Adown2)
#	Ne=Nup1+Nup2+Ndown1+Ndown2
	NT=Nup1+Nup2+Ndown1+Ndown2
	print (Nup1, Nup2, Ndown1, Ndown2, Ne,NT)

#plt.show()
