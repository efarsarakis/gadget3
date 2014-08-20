#!/bin/bash --login

# uncomment the following two lines if you need to debug this script
# set -v      # Print script lines as they are read.
# set -x      # Print commands and their arguments as they are executed.

# PBS job options (name, compute nodes, job time)
#PBS -N Gadget3
#PBS -l nodes=6:ppn=8
#PBS -q gpu_nodes
#PBS -l walltime=00:50:00

# Replace [budget code] below with your project code (e.g. t01)
#commented PBS -A d59

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)               

# Change to the directory that the job was submitted from
# (remember this should be on the /work filesystem)
cd $PBS_O_WORKDIR

# Set the number of threads to 1
#   This prevents any system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

module load cray-mpich
module add perftools
module load fftw/2.1.5.7.1
module load craype-accel-nvidia35

export CUDA_PROFILE=1
export CRAY_CUDA_MPS=1
export LD_LIBRARY_PATH=/lus/scratch/p02045/lib/mygsl/lib:LD_LIBRARY_PATH


make clean

make Config=Config.sh


# Launch the parallel job
#   Using 1536 MPI processes and 24 MPI processes per node
aprun -n 48 -N 8 ./GadgetMedium param.txt
