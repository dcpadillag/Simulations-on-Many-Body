import numpy as np
import matplotlib.pyplot as plt
import sys

t=1./4.4 #Hopping
U=1. #Interaction
beta=333.33 #1/kbT

#Lattice Parameters
L=100#Lattice side length
Ntotal=L*L #total number of sites
n=1.3720699982192717
#Ncell=2*(Ntotal/4)
#Ne=Ntotal/Ncell
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
Ne=0.5 #total number of particles per cell
d=float(sys.argv[1])
Nup=Ne/2+d
Ndown=Ne/2-d

kx=np.linspace(-np.pi/a,np.pi/a,M)
ky=np.linspace(-np.pi/a,np.pi/a,M)

def E(k):
	return 2.*np.cos(np.dot(k,a1))+2.*np.cos(np.dot(k,a2))


def DOS(E,Eup,Edown,dE,mu):
    sigma=1.5*dE
    delta =((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-Eup+mu)/sigma)**2)).sum()
    delta+=((1./(np.sqrt(2.*np.pi*sigma**2)))*np.exp(-((E-Edown+mu)/sigma)**2)).sum()
    return delta


def F(mu,D,E,Ne):
    integral=0.
    dE=(max(E)-min(E))/len(E)
    for i in range(len(E)):
        integral+=(D[i]/(np.exp(beta*(E[i]-mu))+1))*dE
    return integral-Ne


def biseccion(mu1,mu2,Eup,Edown,Ne):
    #Density of State
    E=np.linspace(-20.*t,20.*t,M)
    dE=(max(E)-min(E))/M

    D1=np.array([])
    for i in range(len(E)):
         D1=np.append(D1,DOS(E[i],Eup,Edown,dE,mu1))
    D1=D1/Ntotal
    D2=np.array([])
    for i in range(len(E)):
         D2=np.append(D2,DOS(E[i],Eup,Edown,dE,mu2))
    D2=D2/Ntotal
#    plt.plot(D1,E/t)
#    plt.plot(D2,E/t)
#    plt.show()
    F1=F(mu1,D1,E,Ne)
    F2=F(mu2,D2,E,Ne)
    if (F2*F1<0):
        mu_m=0.5*(mu1+mu2)
        Dm=np.array([])
        for i in range(len(E)):
             Dm=np.append(Dm,DOS(E[i],Eup,Edown,dE,mu_m))
        Dm=Dm/Ntotal
        Fb=F(mu_m,Dm,E,Ne)
        while(abs(Fb)>10.E-6):
            if (Fb>0):
                mu2=mu_m
                mu_m=0.5*(mu1+mu2)
            else:
                mu1=mu_m
                mu_m=0.5*(mu1+mu2)
            Dm=np.array([])
            for i in range(len(E)):
                 Dm=np.append(Dm,DOS(E[i],Eup,Edown,dE,mu_m))
            Dm=Dm/Ntotal
            Fb=F(mu_m,Dm,E,Ne)
    else:
        print ("Wrong Chemical Potential Interval")
	mu_mu=0.0
    return mu_m

def biseccion0(mu1,mu2,D,E,Ne):
    F1=F(mu1,D,E,Ne)
    F2=F(mu2,D,E,Ne)
    if (F2*F1<0):
        mu_m=0.5*(mu1+mu2)
        Fb=F(mu_m,D,E,Ne)
        while(abs(Fb)>10.E-6):
            if (Fb>0):
                mu2=mu_m
                mu_m=0.5*(mu1+mu2)
            else:
                mu1=mu_m
                mu_m=0.5*(mu1+mu2)
            Fb=F(mu_m,D,E,Ne)
    else:
        print ("NO PANA")
	mu_m=0.
    return mu_m

def occupation(mu,Eup,Edown):
    suma1=0.
    suma2=0.

    suma1+=(0.5/(np.exp(beta*(Eup-mu))+1)).sum()
    suma2+=(0.5/(np.exp(beta*(Edown-mu))+1)).sum()
    return suma1/Ntotal, suma2/Ntotal

Kx, Ky = np.meshgrid(kx, ky)
for i in range(10):
	EU=-U*(Nup*Ndown)
	Eup=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)+U*Ndown+EU
	Edown=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)+U*Nup+EU

	mu=biseccion(-5.*t,5.*t,Eup,Edown,Ne)
	print mu

	E=np.linspace(-5.*t,5.*t,M)
	dE=(max(E)-min(E))/M

	D=np.array([])
        for i in range(len(E)):
   	     D=np.append(D,DOS(E[i],Eup,Edown,dE,mu))
	D=D/Ntotal
	N=F(mu,D,E,0)
	print N

	Nup, Ndown = occupation(mu,Eup,Edown)
#	Ne=Nup+Ndown
	print Nup, Ndown, Nup+Ndown


EU=-U*(Nup*Ndown)
Eup=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)+U*Ndown+EU
Edown=-2.*t*np.cos(Kx*a)-2.*t*np.cos(Ky*a)+U*Nup+EU

Energy=Eup.sum()+Edown.sum()
print Energy/Ntotal
plt.plot(D,E/t)
plt.show()
#cp = plt.contourf(Kx, Ky, Eup)
#plt.colorbar(cp)
#plt.show()

#Density of State
#E=np.linspace(-5.*t,5.*t,M)
#dE=(max(E)-min(E))/len(E)
#D=np.array([])
#for i in range(len(E)):
#	D=np.append(D,DOS(E[i],Eup,Edown,dE,0))
#D=D/Ntotal
#plt.plot(D,E/t)
#plt.show()

#N=0.
#for i in range(len(E)):
#	N+=D[i]*dE
#print N


