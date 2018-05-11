# Stochastic Series Expansion (SSE) Quantum Monte Carlo (QMC)

## Description
A QMC algorithm employing Sandvik's SSE method implemented in Fortran 90, adapted for a Bose-Hubbard model. 

See [Sandvik, 1999](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.59.R14157) on the algorithm, and [Sandvik, 1997](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.11678) for quantities measured in this program.

Program is split into 3 .f90 files: gfmain, gfhb and gfupd. gfmain contains the main program. Initialization of operating string, lattice, probability tables etc. are split between gfmain and gfhb. gfupd contains the main operator-loop updates.

## Installation
It is best to deploy the code on a high-performance computing cluster, although it is possible to run it on your local machine for a smaller number of Monte Carlo steps at smaller system sizes. Computing clusters should come with the `intel` compiler, which can be compiled as
```
intel -o [name of executable] gfmain_16x16.f90 gfupd.f90 gfhb.f90
```
for the 16x16 lattice. If you wish to run the simulation on your local machine, it is advisable to download the free GNU Fortran compiler `gfortran`, and compile the program as
```
gfortran -o [name of executable] gfmain_16x16.f90 gfupd.f90 gfhb.f90
```


