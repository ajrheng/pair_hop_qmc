#!/bin/bash
RANDOM=$$
for i in {0..4..2}
do
	
		mkdir "mu=$i"
		cp script.sh mu=$i
		cp 16x16code mu=$i 
		cd mu=$i
		echo 1 >> read.in
		echo 1 >> read.in
		echo 0 >> read.in
		echo 8 >> read.in
		echo $i >> read.in
		echo 16 >> read.in
		echo 1000000 >> read.in
		echo 10 >> read.in
		echo 100000 >> read.in
		echo 0 >> read.in
		echo $RANDOM >> rand.in
		echo $RANDOM >> rand.in
		echo $RANDOM >> rand.in
		echo $RANDOM >> rand.in
		qsub script.sh
		cd ..
		
	
done
exit 0
