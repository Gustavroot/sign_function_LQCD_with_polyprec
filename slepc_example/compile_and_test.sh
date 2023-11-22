# load modules on JUWELS
#. ../load_JUWELS_modules.sh
# compile
#mpicc slepc_test.c -DPETSC_USE_COMPLEX -lpetsc -lslepc -lm
# test
#sbatch subSLEPC_TEST.sh

make
sbatch subSLEPC_TEST.sh
