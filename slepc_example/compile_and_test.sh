# load modules on JUWELS
. ../load_JUWELS_modules.sh

# compile
mpicc slepc_test.c -lpetsc -lslepc -lm

# test
sbatch subSLEPC_TEST.sh

#./a.out -fn_scale 1.0,1.0 -n 10 -fn_method 0 -verbose 1
