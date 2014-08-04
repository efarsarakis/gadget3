#!/bin/bash
#$ -S /bin/bash
#$ -j n
#$ -cwd
#$ -pe mvapich 96
#$ -q standard.q
#$ -m be
#$ -M volker.springel@h-its.org
#$ -N test
##$ -binding linear:32
#$ -l h_rt=24:00:00
#

cat $PE_HOSTFILE

source /etc/profile.d/modules.sh
module load sge
module load mvapich2/gcc/64/1.5.1-qlc
module load hydra


export OMP_NUM_THREADS=6
export NTASKS=$((NSLOTS/OMP_NUM_THREADS))
export TASKS_PER_HOST=$((NTASKS/NHOSTS))

mpistart -np $NTASKS -ppn $TASKS_PER_HOST  ./P-Gadget3 param.txt 


if [ -f "../cont" ]; then
  echo 
  
#  rm -f "../cont"
#  qsub job-magny.sh

fi





