#!/bin/sh
#$ -j n
#$ -cwd
#$ -pe mvapich 16
#$ -q standard.q
#$ -m be
#$ -M volker.springel@h-its.org
#$ -N run0
###$ -binding linear:16
#$ -l h_rt=24:00:00
#

cat $PE_HOSTFILE

. /etc/profile.d/modules.sh
module load sge
module load mvapich2/gcc/64/1.5.1-qlc
module load hydra

mpistart -np $NSLOTS  ./P-Gadget3 param.txt 


if [ -f "../cont" ] 
then
#    rm -f "../cont"
#    qsub job-magny.sh
fi





