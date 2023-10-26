#!/bin/bash

#SBATCH --account=hwu29
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=mpi_out_%j.txt
#SBATCH --error=mpi_err_%j.txt
#SBATCH --time=00:59:00
#SBATCH --partition=batch



srun a.out -fn_scale 1.0,1.0 -n 10 -fn_method 0 -verbose 1
