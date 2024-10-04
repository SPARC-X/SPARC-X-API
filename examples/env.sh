#!/bin/bash
# set path to any shared libraries such as mkl, blas, scalapack etc.
ml intel-oneapi-compilers intel-oneapi-mkl intel-oneapi-mpi
source ~/p-amedford6-0/venvs/joss-sparc-x-api/bin/activate

cd ~
export HOME=$(pwd)
cd -
export SCRATCH=$HOME/scratch
export DATA=$HOME/p-amedford6-0

export PYTHONPATH=$DATA/active/JOSS-SPARC-X-API/SPARC-X-API:$PYTHONPATH

# Check if Slurm is installed and running
if command -v sinfo &> /dev/null; then
  export ASE_SPARC_COMMAND="srun ${DATA}/active/JOSS-SPARC-X-API/dev_SPARC/lib/sparc --export=ALL -K"
else
  nprocs=$(nproc)
  export ASE_SPARC_COMMAND="mpirun -np $nprocs ${DATA}/SPARC/lib/sparc --export=ALL -K" 
fi
# Set pseudopotential path is you did not download with python -m sparc.download_data
export PLUMED_KERNEL=/storage/home/hcoda1/9/ltimmerman3/p-amedford6-0/active/JOSS-SPARC-X-API/plumed-2.9.2/src/lib/libplumedKernel.so
