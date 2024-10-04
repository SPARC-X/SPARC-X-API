#!/bin/bash
#SBATCH --job-name=socket_tests
#SBATCH --account=gts-amedford6
#SBATCH --mail-user=ltimmerman3@gatech.edu
#SBATCH -N 1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH --time=8:00:00
#SBATCH --output=Report_%A-%a.out
#SBATCH -qembers

cd $SLURM_SUBMIT_DIR
source ~/.bashrc
source env.sh

# cd socket/singlepoint
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd socket/relax/ab-initio/cell
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd socket/relax/ab-initio/coords
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd socket/relax/mlff/coords
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd socket/md/ab-initio/ase
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd socket/md/ab-initio/plumed
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd socket/md/mlff/plumed
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# FileIO

# cd FileIO/relax/relax_cell/ab-initio
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd FileIO/relax/relax_coords/ab-initio
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd FileIO/relax/relax_coords/mlff
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

# cd FileIO/md/ab-initio
# python run.py > run.out
# cd $SLURM_SUBMIT_DIR

cd FileIO/md/mlff
python run.py > run.out
cd $SLURM_SUBMIT_DIR