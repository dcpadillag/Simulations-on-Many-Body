import numpy as np
import matplotlib.pyplot as plt

a=open("results.txt",'r')

tF=[]
NeF=[]

tP=[]
NeP=[]

for i in a.readlines():
	lista=i.split()
	phase=str(lista[2])
	if phase == "F":
		tF.append(float(lista[0]))
		NeF.append(float(lista[1]))
	else:
		tP.append(float(lista[0]))
		NeP.append(float(lista[1]))

#plt.plot(NeF,tF,'s',color='blue',label="$Ferromagnetico$")
#plt.plot(NeP,tP,'o',color='red',label="$Paramagnetico$")
plt.plot(NeF,tF,'s',color='blue')
plt.plot(NeP,tP,'o',color='red')
#plt.legend(loc='best')
plt.xlabel("$n_{e}$",fontsize=18)
plt.ylabel("$t/U$",fontsize=18)
plt.xlim([0.05,1.05])
plt.ylim([-0.05,0.8])
plt.savefig("phase_diagram_triangular.pdf")
plt.show()
