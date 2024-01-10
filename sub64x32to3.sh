#!/bin/bash

#SBATCH --account=mul-tra
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=24
#SBATCH --threads-per-core=1
#SBATCH --output=mpi_out_%j.txt
#SBATCH --error=mpi_err_%j.txt
#SBATCH --time=04:59:00
#SBATCH --partition=batch

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}



. run16 -i sample_64x32to3.ini
