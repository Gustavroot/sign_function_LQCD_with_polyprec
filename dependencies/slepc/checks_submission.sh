#!/bin/bash

#SBATCH --account=mul-tra
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=mpi_out_%j.txt
#SBATCH --error=mpi_err_%j.txt
#SBATCH --time=00:19:00
#SBATCH --partition=batch


export BASE_SLEPC_DIR=`pwd`
make SLEPC_DIR=$BASE_SLEPC_DIR/../dir PETSC_DIR=/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex PETSC_ARCH="" check
