#!/bin/bash

#SBATCH --account=mul-tra
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=mpi_out_%j.txt
#SBATCH --error=mpi_err_%j.txt
#SBATCH --time=00:59:00
#SBATCH --partition=batch



srun slepc_test
