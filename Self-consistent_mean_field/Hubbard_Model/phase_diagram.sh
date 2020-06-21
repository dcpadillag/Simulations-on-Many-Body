#!/bin/bash

echo "-> start"
rm -r triangular
rm results.txt
mkdir triangular

N=$1

for i in $(seq 0 38);do
	for j in $(seq 1 10);do
		echo $i,$j
		time python hubbard_mf_triangular.py $i 1.0 $N $j >> results.txt
	done
done

echo "-> finished"
#exit
