#!/bin/sh
#$ -j n
#$ -cwd
#$ -pe mvapich2 32
#$ -v THREADS_PER_MPI_TASK=4      #  affects placement of MPI tasks over nodes
#$ -m be
#$ -M volker@mpa-garching.mpg.de
#$ -N P3
#$ -l h_rt=24:00:00
#

# In this example, 32 Cores are used, with 8 MPI Tasks and 4 threads each

module unload mvapich2-1.2-sdr-gnu/4.1.2

module load intel/11.0
module load mvapich2-1.2-sdr-intel/11.0


#disables pinning of Tasks and their threads to CPUs
export MV2_ENABLE_AFFINITY=0

#tells OpenMP how many threads it should use in OpenMP directives
export OMP_NUM_THREADS=4


mpiexec -np 8  ./P-Gadget3  param.txt

if [ -f "../cont" ]
then
#    rm -f "../cont"
#    qsub job-opa-smp.sh
fi
