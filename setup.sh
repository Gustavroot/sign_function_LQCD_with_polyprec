# IMPORTANT : run this by doing : . setup.sh A B

#             where A,B are 0 or 1
#             A : if 1, install SLEPc from scratch
#             B : if 1, install LAPACKE from scratch

if [ $1 -eq 1 ]
then
    echo "Installing SLEPc"
    # own local installation of SLEPc
    cd dependencies/
    . load_JUWELS_modules.sh
    unset SLEPC_DIR
    cd slepc/
    . installer.sh
    cd ../../
fi

if [ $2 -eq 1 ]
then
    echo "Installing LAPACKE"
    # own local installation of LAPACKE and BLAS
    cd dependencies/
    . load_JUWELS_modules.sh
    rm -Rf lapack-3.9.0 v3.9.0.tar.gz
    . install_LAPACKE_and_BLAS.sh
fi

if [ $1 -ne 1 -a $2 -ne 1 ]
then
    echo "Not installing SLEPc or LAPACKE, loading modules only"
    cd dependencies/
    . load_JUWELS_modules.sh
    cd ../
fi
