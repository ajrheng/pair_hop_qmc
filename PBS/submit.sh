#!/bin/bash

#PBS -P Personal
#PBS -q normal 
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l place=free
#PBS -l walltime=120:00:00
#PBS -V
#PBS -m bea

## pre-processing script
cd $PBS_O_WORKDIR
./16x16testBed
