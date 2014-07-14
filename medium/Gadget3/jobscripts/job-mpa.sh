#!/usr/local/bin/tcsh
#$ -cwd
#$ -e llrunme0.err
#$ -o llrunme0.out
#$ -pe openib 4
#$ -l h_vmem=2048M
#$ -l arch=*-amd64
#$ -N GalaxyB

mpiexec P-Gadget3/P-Gadget3 galaxy.par
