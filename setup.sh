# run this by doing : . setup.sh

# own local installation of SLEPc
cd dependencies/
. load_JUWELS_modules.sh
unset SLEPC_DIR
cd slepc/
. installer.sh
cd ../../

# own local installation of LAPACKE and BLAS
cd dependencies/
. load_JUWELS_modules.sh
rm -Rf lapack-3.9.0 v3.9.0.tar.gz
. install_LAPACKE_and_BLAS.sh

