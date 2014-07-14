#!/bin/tcsh
#$ -j n
#$ -cwd
#$ -pe impi4 8
#$ -m be
#$ -M username@mpa-garching.mpg.de
#$ -N test
#$ -l h_rt=...

mpiexec -np $NSLOTS  ./P-Gadget3  param.txt 

if [ -f ../cont ] 
then
 rm -f "../cont"
 qsub job-opa.sh
fi
    





