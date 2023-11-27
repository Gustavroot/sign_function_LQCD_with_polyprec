# --- COMPILER ----------------------------------------
CC = mpicc -std=gnu99 -Wall -pedantic
CPP = cpp
MAKEDEP = $(CPP) -MM

# --- DO NOT CHANGE -----------------------------------
SRCDIR = src
BUILDDIR = build
GSRCDIR = $(BUILDDIR)/gsrc
SRC = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.c,$(wildcard $(SRCDIR)/*.c)))
SRCGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.c))
GSRCFLT = $(patsubst %_generic.c,$(GSRCDIR)/%_float.c,$(SRCGEN))
GSRCDBL = $(patsubst %_generic.c,$(GSRCDIR)/%_double.c,$(SRCGEN))
GSRC = $(patsubst %,$(GSRCDIR)/%,$(SRC)) $(GSRCFLT) $(GSRCDBL)
HEA = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.h,$(wildcard $(SRCDIR)/*.h)))
HEAGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.h))
GHEAFLT = $(patsubst %_generic.h,$(GSRCDIR)/%_float.h,$(HEAGEN))
GHEADBL = $(patsubst %_generic.h,$(GSRCDIR)/%_double.h,$(HEAGEN))
GHEA = $(patsubst %,$(GSRCDIR)/%,$(HEA)) $(GHEAFLT) $(GHEADBL)
OBJ = $(patsubst $(GSRCDIR)/%.c,$(BUILDDIR)/%.o,$(GSRC))
OBJDB = $(patsubst %.o,%_db.o,$(OBJ))
DEP = $(patsubst %.c,%.dep,$(GSRC))

# --- FLAGS -------------------------------------------
# make this empty to disable POLYPREC and/or SLEPC
POLYPREC_ENABLER = #-DPOLYPREC
SLEPC_ENABLER = -DUSE_SLEPC

OPT_FLAGS = -fopenmp -msse4.2 -DORTH_MGS $(POLYPREC_ENABLER) $(SLEPC_ENABLER) #-DSSE -DOPENMP
CFLAGS = -DPARAMOUTPUT -DTRACK_RES -DFGMRES_RESTEST -DPROFILING
# -DSINGLE_ALLREDUCE_ARNOLDI
# -DCOARSE_RES -DSCHWARZ_RES -DTESTVECTOR_ANALYSIS
OPT_VERSION_FLAGS =$(OPT_FLAGS) -O3 -ffast-math
DEBUG_VERSION_FLAGS = $(OPT_FLAGS)

# --- FLAGS FOR HDF5 ---------------------------------
# H5HEADERS=-DHAVE_HDF5 /usr/include
# H5LIB=-lhdf5 -lz

# --- FLAGS FOR LIME ---------------------------------
# LIMEH=-DHAVE_LIME -I$(LIMEDIR)/include
# LIMELIB= -L$(LIMEDIR)/lib -llime

# use own installation of LAPACKE in case of compiling with polynomial preconditioner
ifeq ($(POLYPREC_ENABLER),-DPOLYPREC)
LAPACK_DIR = dependencies/lapack-3.9.0
LAPACKE_DIR = $(LAPACK_DIR)/LAPACKE
LAPACKE_INCLUDE = -I$(LAPACKE_DIR)/include
BLASLIB      = $(LAPACK_DIR)/librefblas.a
LAPACKLIB    = $(LAPACK_DIR)/liblapack.a
LAPACKELIB   = $(LAPACK_DIR)/liblapacke.a
LAPACK_LIBRARIES = $(LAPACKELIB) $(LAPACKLIB) $(BLASLIB)
else
LAPACKE_INCLUDE =
LAPACK_LIBRARIES =
endif

ifeq ($(SLEPC_ENABLER),-DUSE_SLEPC)
#include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
SLEPC_FLAGS_COMP = -O2 -ftree-vectorize -march=native -fno-math-errno -fPIC -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64  -I/p/home/jusers/ramirez1/juwels/projects/sign_function_LQCD_with_polyprec/dependencies/slepc/dir/include -I/p/home/jusers/ramirez1/juwels/projects/sign_function_LQCD_with_polyprec/dependencies/slepc/dir/include  -I/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex/include -I/p/software/juwels/stages/2024/software/ParMETIS/4.0.3-gpsmpi-2023a/include -I/p/software/juwels/stages/2024/software/METIS/5.1.0-GCC-12.3.0/include -I/p/software/juwels/stages/2024/software/HDF5/1.14.2-gpsmpi-2023a/include
SLEPC_FLAGS_LINK = -Wl,-rpath,/p/home/jusers/ramirez1/juwels/projects/sign_function_LQCD_with_polyprec/dependencies/slepc/dir/lib -L/p/home/jusers/ramirez1/juwels/projects/sign_function_LQCD_with_polyprec/dependencies/slepc/dir/lib -lslepc  -Wl,-rpath,/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex/lib -L/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex/lib -L/p/software/juwels/stages/2024/software/PETSc/3.20.0-gpsfbf-2023a-complex/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/ScaLAPACK/2.2.0-gpsmpi-2023a-fb/lib -L/p/software/juwels/stages/2024/software/ScaLAPACK/2.2.0-gpsmpi-2023a-fb/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/FlexiBLAS/3.3.1-GCC-12.3.0/lib -L/p/software/juwels/stages/2024/software/FlexiBLAS/3.3.1-GCC-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/ParMETIS/4.0.3-gpsmpi-2023a/lib -L/p/software/juwels/stages/2024/software/ParMETIS/4.0.3-gpsmpi-2023a/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/METIS/5.1.0-GCC-12.3.0/lib -L/p/software/juwels/stages/2024/software/METIS/5.1.0-GCC-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/HDF5/1.14.2-gpsmpi-2023a/lib -L/p/software/juwels/stages/2024/software/HDF5/1.14.2-gpsmpi-2023a/lib -Wl,-rpath,/opt/ddn/ime/lib -L/opt/ddn/ime/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/psmpi/5.9.2-1-GCC-12.3.0/lib -L/p/software/juwels/stages/2024/software/psmpi/5.9.2-1-GCC-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/ParMETIS/4.0.3-gpsmpi-2023a/lib64 -L/p/software/juwels/stages/2024/software/ParMETIS/4.0.3-gpsmpi-2023a/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/METIS/5.1.0-GCC-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/METIS/5.1.0-GCC-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/HDF5/1.14.2-gpsmpi-2023a/lib64 -L/p/software/juwels/stages/2024/software/HDF5/1.14.2-gpsmpi-2023a/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/Szip/2.1.1-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/Szip/2.1.1-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/libarchive/3.6.2-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/libarchive/3.6.2-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/cURL/8.0.1-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/cURL/8.0.1-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/bzip2/1.0.8-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/bzip2/1.0.8-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/ncurses/6.4-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/ncurses/6.4-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/ScaLAPACK/2.2.0-gpsmpi-2023a-fb/lib64 -L/p/software/juwels/stages/2024/software/ScaLAPACK/2.2.0-gpsmpi-2023a-fb/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/FFTW.MPI/3.3.10-gpsmpi-2023a/lib64 -L/p/software/juwels/stages/2024/software/FFTW.MPI/3.3.10-gpsmpi-2023a/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/FFTW/3.3.10-GCC-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/FFTW/3.3.10-GCC-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/FlexiBLAS/3.3.1-GCC-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/FlexiBLAS/3.3.1-GCC-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/OpenBLAS/0.3.23-GCC-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/OpenBLAS/0.3.23-GCC-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/BLIS/0.9.0-GCC-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/BLIS/0.9.0-GCC-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/psmpi/5.9.2-1-GCC-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/psmpi/5.9.2-1-GCC-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/PMIx/4.2.6-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/PMIx/4.2.6-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/hwloc/2.9.1-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/hwloc/2.9.1-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/libpciaccess/0.17-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/libpciaccess/0.17-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/libevent/2.1.12-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/libevent/2.1.12-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/OpenSSL/1.1/lib64 -L/p/software/juwels/stages/2024/software/OpenSSL/1.1/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/libxml2/2.11.4-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/libxml2/2.11.4-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/XZ/5.4.2-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/XZ/5.4.2-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/pscom/5-default-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/pscom/5-default-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/UCX/default-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/UCX/default-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/CUDA/12/stubs/lib64 -L/p/software/juwels/stages/2024/software/CUDA/12/stubs/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/numactl/2.0.16-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/numactl/2.0.16-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/binutils/2.40-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/binutils/2.40-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/zlib/1.2.13-GCCcore-12.3.0/lib64 -L/p/software/juwels/stages/2024/software/zlib/1.2.13-GCCcore-12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/GCCcore/12.3.0/lib/gcc/x86_64-pc-linux-gnu/12.3.0 -L/p/software/juwels/stages/2024/software/GCCcore/12.3.0/lib/gcc/x86_64-pc-linux-gnu/12.3.0 -Wl,-rpath,/p/software/juwels/stages/2024/software/GCCcore/12.3.0/lib64 -L/p/software/juwels/stages/2024/software/GCCcore/12.3.0/lib64 -Wl,-rpath,/p/software/juwels/stages/2024/software/Szip/2.1.1-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/Szip/2.1.1-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/libarchive/3.6.2-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/libarchive/3.6.2-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/cURL/8.0.1-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/cURL/8.0.1-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/bzip2/1.0.8-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/bzip2/1.0.8-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/ncurses/6.4-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/ncurses/6.4-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/FFTW.MPI/3.3.10-gpsmpi-2023a/lib -L/p/software/juwels/stages/2024/software/FFTW.MPI/3.3.10-gpsmpi-2023a/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/FFTW/3.3.10-GCC-12.3.0/lib -L/p/software/juwels/stages/2024/software/FFTW/3.3.10-GCC-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/imkl/2023.2.0/mkl/2023.2.0/lib/intel64 -L/p/software/juwels/stages/2024/software/imkl/2023.2.0/mkl/2023.2.0/lib/intel64 -Wl,-rpath,/p/software/juwels/stages/2024/software/imkl/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -L/p/software/juwels/stages/2024/software/imkl/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -Wl,-rpath,/p/software/juwels/stages/2024/software/OpenBLAS/0.3.23-GCC-12.3.0/lib -L/p/software/juwels/stages/2024/software/OpenBLAS/0.3.23-GCC-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/BLIS/0.9.0-GCC-12.3.0/lib -L/p/software/juwels/stages/2024/software/BLIS/0.9.0-GCC-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/PMIx/4.2.6-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/PMIx/4.2.6-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/hwloc/2.9.1-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/hwloc/2.9.1-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/libpciaccess/0.17-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/libpciaccess/0.17-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/libevent/2.1.12-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/libevent/2.1.12-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/OpenSSL/1.1/lib -L/p/software/juwels/stages/2024/software/OpenSSL/1.1/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/libxml2/2.11.4-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/libxml2/2.11.4-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/XZ/5.4.2-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/XZ/5.4.2-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/pscom/5-default-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/pscom/5-default-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/UCX/default-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/UCX/default-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/CUDA/12/lib -L/p/software/juwels/stages/2024/software/CUDA/12/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/numactl/2.0.16-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/numactl/2.0.16-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/binutils/2.40-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/binutils/2.40-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/zlib/1.2.13-GCCcore-12.3.0/lib -L/p/software/juwels/stages/2024/software/zlib/1.2.13-GCCcore-12.3.0/lib -Wl,-rpath,/p/software/juwels/stages/2024/software/GCCcore/12.3.0/lib -L/p/software/juwels/stages/2024/software/GCCcore/12.3.0/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lflexiblas -lgfortran -lsuperlu -lsuperlu_dist -lflexiblas -lgfortran -lpthread -lparmetis -lmetis -lhdf5_hl -lhdf5 -lchaco -ltriangle -lm -lim_client -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -lquadmath -ldl
endif

all:: wilson library #documentation
wilson: dd_alpha_amg #dd_alpha_amg_db
library: lib/libdd_alpha_amg.a include/dd_alpha_amg_parameters.h include/dd_alpha_amg.h
documentation: doc/user_doc.pdf

.PHONY: all wilson library
.SUFFIXES:
.SECONDARY:

dd_alpha_amg : $(OBJ)
ifeq ($(SLEPC_ENABLER),-DUSE_SLEPC)
	$(CC) $(OPT_VERSION_FLAGS) $(LIMEH) -o $@ $(OBJ) $(H5LIB) $(LIMELIB) $(LAPACK_LIBRARIES) $(SLEPC_FLAGS_LINK) -lm -lgfortran
else
	$(CC) $(OPT_VERSION_FLAGS) $(LIMEH) -o $@ $(OBJ) $(H5LIB) $(LIMELIB) $(LAPACK_LIBRARIES) -lm -lgfortran
endif

dd_alpha_amg_db : $(OBJDB)
	$(CC) -g $(DEBUG_VERSION_FLAGS) $(LIMEH) -o $@ $(OBJDB) $(H5LIB) $(LIMELIB) $(LAPACK_LIBRARIES) -lm -lgfortran

lib/libdd_alpha_amg.a: $(OBJ)
	ar rc $@ $(OBJ)
	ar d $@ main.o
	ranlib $@

doc/user_doc.pdf: doc/user_doc.tex doc/user_doc.bib
	( cd doc; pdflatex user_doc; bibtex user_doc; pdflatex user_doc; pdflatex user_doc; )

include/dd_alpha_amg.h: src/dd_alpha_amg.h
	cp src/dd_alpha_amg.h $@

include/dd_alpha_amg_parameters.h: src/dd_alpha_amg_parameters.h
	cp src/dd_alpha_amg_parameters.h $@

$(BUILDDIR)/%.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
ifeq ($(SLEPC_ENABLER),-DUSE_SLEPC)
	$(CC) $(CFLAGS) $(OPT_VERSION_FLAGS) $(LAPACKE_INCLUDE) $(H5HEADERS) $(LIMEH) $(SLEPC_FLAGS_COMP) -c $< -o $@
else
	$(CC) $(CFLAGS) $(OPT_VERSION_FLAGS) $(LAPACKE_INCLUDE) $(H5HEADERS) $(LIMEH) -c $< -o $@
endif

$(BUILDDIR)/%_db.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) -g $(CFLAGS) $(DEBUG_VERSION_FLAGS) $(LAPACKE_INCLUDE) $(H5HEADERS) $(LIMEH) -DDEBUG -c $< -o $@

$(GSRCDIR)/%.h: $(SRCDIR)/%.h $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

$(GSRCDIR)/%.c: $(SRCDIR)/%.c $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

%.dep: %.c $(GHEA)
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1.o $@ : ,g' > $@
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1_db.o $@ : ,g' >> $@
clean::
	rm -f $(BUILDDIR)/*.o
	rm -f $(GSRCDIR)/*
	rm -f dd_alpha_amg
	rm -f dd_alpha_amg_db

-include $(DEP)
