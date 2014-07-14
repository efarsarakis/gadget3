#!/bin/sh
#$ -j n
#$ -cwd
#$ -pe fullnode-mpich2-* 32
#$ -v THREADS_PER_MPI_TASK=1
#$ -m be
#$ -M volker@mpa-garching.mpg.de
#$ -N run0
#$ -l h_rt=24:00:00
#

echo $NSLOTS
echo $TMPDIR/machines
cat $TMPDIR/machines


/opt/software/mvapich2-1.4/bin/mpiexec -np $NSLOTS  ./P-Gadget3 param.txt 


if [ -f "../cont" ] 
then
#    rm -f "../cont"
#    qsub job-mbm.sh
fi





