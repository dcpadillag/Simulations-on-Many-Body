#!/bin/bash

echo "-> start"

#rm dmrg dmrg.o
#make

t=$1
U=$2
D=$3
geo=$4 #lattice geometry --> 0: Chain, 1:Honeycomb Ly=2, 2:Honeycomb Ly=4
Ly=$5
measure=$6 

#sed -i "6s/.*/    lattice = $geo/" datos
sed -i "7s/.*/    t1 = $t/" datos
sed -i "8s/.*/    U = $U/" datos
sed -i "9s/.*/    D = $D/" datos

sed -i "11s/.*/    quiet = yes/" datos

if [ $geo == 0 ]; then
#	arg=(64 128 200 256 300) #L for chain
	arg=(8 10)
else
	arg=(4 7) #Lx for honeycomb geometry
#	arg=(10 13 16 19 22 25)
fi

for Lx in ${arg[*]};do
	rm -r $(($geo))_t_$(bc -l <<< "scale=2;$t")_L_$(($Lx))_U_$(bc -l <<< "scale=2;$U")_D_$D
	mkdir $(($geo))_t_$(bc -l <<< "scale=2;$t")_L_$(($Lx))_U_$(bc -l <<< "scale=2;$U")_D_$D
	
	if [ $measure == 1 ]; then
		cd Measurements
		cp measure measure.o ..
		cd ..
		mv measure measure.o $(($geo))_t_$(bc -l <<< "scale=2;$t")_L_$(($Lx))_U_$(bc -l <<< "scale=2;$U")_D_$D
	else
		echo "No measurements"
 
	fi 

	cp -a datos dmrg dmrg.o $(($geo))_t_$(bc -l <<< "scale=2;$t")_L_$(($Lx))_U_$(bc -l <<< "scale=2;$U")_D_$D

	cd $(($geo))_t_$(bc -l <<< "scale=2;$t")_L_$(($Lx))_U_$(bc -l <<< "scale=2;$U")_D_$D

	if [ $geo == 0 ]; then
		L=$Lx #for chain
	else
		L=$(($Ly+$Ly*$Lx*2)) #for honeycomb geometry
	fi

	N=$(($L/2))

	sed -i "3s/.*/    N = $L/" datos
	sed -i "4s/.*/    Nup = $N/" datos
	sed -i "5s/.*/    Ndown = $N/" datos

	echo 'Geometry='$geo, 't='$t, 'U='$U, 'Delta='$D, 'Lx='$Lx, 'L='$L 


	(./dmrg datos > $N$N.txt

################################SPIN_GAP###################################
#	sed -i "4s/.*/    Nup = $(($N+1))/" datos
#	sed -i "5s/.*/    Ndown = $(($N-1))/" datos

#	./dmrg datos > $(($N+1))$(($N-1)).txt

#	rm dmrg dmrg.o) &
#	cd ..
#############################################################################
################################CHARGE_GAP###################################
	sed -i "4s/.*/    Nup = $(($N+1))/" datos
	sed -i "5s/.*/    Ndown = $(($N))/" datos
	./dmrg datos > $(($N+1))$N.txt

	sed -i "4s/.*/    Nup = $(($N-1))/" datos
	sed -i "5s/.*/    Ndown = $(($N))/" datos
	./dmrg datos > $(($N-1))$N.txt
	rm dmrg dmrg.o
##############################################################################	
########################Do Measurements#######################################	
	if [ $measure == 1 ]; then
		sed -i "4s/.*/    Nup = $N/" datos
		sed -i "5s/.*/    Ndown = $N/" datos
		./measure datos > measure_$N$N.txt
		rm measure measure.o
	fi
	) &
	cd ..
##############################################################################
done

echo "-> finished"
#exit
