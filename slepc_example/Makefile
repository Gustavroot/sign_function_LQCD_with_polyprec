default: slepc_test

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

slepc_test: slepc_test.o
	mpicc -o slepc_test slepc_test.o ${SLEPC_EPS_LIB}
	${RM} slepc_test.o
