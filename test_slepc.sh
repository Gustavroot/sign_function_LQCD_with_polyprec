# run this by doing : . test_slepc.sh

cd dependencies/
. load_JUWELS_modules.sh
cd ../

cd slepc_example/
. compile_and_test.sh
export BASE_SLEPC_DIR=`pwd`
echo "Check $BASE_SLEPC_DIR for mpi_out and mpi_err from invsqrt SLEPc test"
cd ../
