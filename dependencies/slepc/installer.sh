# EasyBuild installation
#export USERINSTALLATIONS=/p/home/jusers/ramirez1/juwels/installs/slepc/dir/
#module --force purge
#module load Stages/2024
#ml UserInstallations
#eb SLEPc-3.20.0-gpsfbf-2023a-complex.eb


# Manual Installation

export BASE_SLEPC_DIR=`pwd`
rm -Rf slepc-*
rm -Rf dir/
mkdir dir/
wget https://slepc.upv.es/download/distrib/slepc-3.20.0.tar.gz
tar -xvzf slepc-3.20.0.tar.gz
export SLEPC_DIR=$BASE_SLEPC_DIR/slepc-3.20.0
# no need to set PETSC_DIR, already set
# we're doing a prefix-based installation of SLEPc, so no PETSC_ARCH setting needed
cd $SLEPC_DIR
./configure --prefix=$BASE_SLEPC_DIR/dir
make -j 8
make SLEPC_DIR=$BASE_SLEPC_DIR/slepc-3.20.0 PETSC_DIR=/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex install
sbatch ../checks_submission.sh
cd ..
echo "Check $SLEPC_DIR for mpi_out and mpi_err to see the output of <make check>"
