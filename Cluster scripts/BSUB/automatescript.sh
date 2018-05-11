#!/bin/bash
RANDOM=$$
for i in {0..64..2}
do
	
		mkdir "mu=$i"
		cp script.sh mu=$i
		cp 16x16cy2001 mu=$i 
		cd mu=$i
		echo 1 >> read.in
		echo 0 >> read.in
		echo 4 >> read.in
		echo 4 >> read.in
		echo $i >> read.in
		echo 32 >> read.in
		echo 10000 >> read.in
		echo 10 >> read.in
		echo 10000 >> read.in
		echo 0 >> read.in
		echo $RANDOM > rand.in
		echo $RANDOM >> rand.in
		echo $RANDOM >> rand.in
		echo $RANDOM >> rand.in
		bsub < script.sh
		cd ..
		
	
done
exit 0
