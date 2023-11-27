# run this by doing : . test_polyprec_invsqrt.sh

cd dependencies/
. load_JUWELS_modules.sh
cd ../

make -j 8
sbatch sub64x32to3.sh
