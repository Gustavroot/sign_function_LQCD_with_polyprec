/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepcfn.h>
#include <complex.h>

typedef double complex complex_double;
typedef double complex *vector_double;

typedef struct {
    vector_double x, b, r, w, wy, wx, wz, *V, *Z;
    complex_double **H, *y, *gamma, *c, *s, shift;
    double tol;
    int num_restart, restart_length, timing, print, kind,
        initial_guess_zero, layout, v_start, v_end, total_storage;
} gmres_double_struct;

gmres_double_struct gp;

static char help[] = "Test matrix inverse square root.\n\n";

/*
   Compute matrix inverse square root B = inv(sqrtm(A))
   Check result as norm(B*B*A-I)
*/

PetscErrorCode MatInvSqrt(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscScalar    tau,eta;
  PetscReal      nrm;
  PetscBool      set,flg;
  PetscInt       n;
  Mat            S,R,Acopy;
  Vec            v,f0;

  PetscFunctionBeginUser;
  PetscCall(MatGetSize(A,&n,NULL));
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&S));
  PetscCall(PetscObjectSetName((PetscObject)S,"S"));
  PetscCall(FNGetScale(fn,&tau,&eta));

  PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&Acopy));
  PetscCall(FNEvaluateFunctionMat(fn,A,S));
  // check that A has not been modified
  PetscCall(MatAXPY(Acopy,-1.0,A,SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(Acopy,NORM_1,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the input matrix has changed by %g\n",(double)nrm));
  PetscCall(MatDestroy(&Acopy));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n"));
    PetscCall(MatView(A,viewer));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed inv(sqrtm(A)) - - - - - - -\n"));
    PetscCall(MatView(S,viewer));
  }

  // check error ||S*S*A-I||_F
  PetscCall(MatMatMult(S,S,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R));
  if (eta!=1.0) PetscCall(MatScale(R,1.0/(eta*eta)));
  PetscCall(MatCreateVecs(A,&v,&f0));
  PetscCall(MatGetColumnVector(S,f0,0));
  PetscCall(MatCopy(R,S,SAME_NONZERO_PATTERN));
  PetscCall(MatDestroy(&R));
  if (tau!=1.0) PetscCall(MatScale(S,tau));
  PetscCall(MatMatMult(S,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R));
  PetscCall(MatShift(R,-1.0));
  PetscCall(MatNorm(R,NORM_FROBENIUS,&nrm));
  PetscCall(MatDestroy(&R));
  if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S*A-I||_F < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S*A-I||_F = %g\n",(double)nrm));
  // check FNEvaluateFunctionMatVec()
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&S));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&f0));
}

PetscErrorCode small_dense_invsqrt( int argcx,char **argvx, gmres_double_struct *p )
{
  FN             fn;
  Mat            A=NULL;
  PetscInt       i,j;
  PetscScalar    x,y,yp,*As;
  PetscViewer    viewer;
  PetscBool      verbose,inplace,matcuda;
  PetscRandom    myrand;
  PetscReal      v;
  char           strx[50],str[50];

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argcx,&argvx,(char*)0,help));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-matcuda",&matcuda));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix inverse square root, n=%" PetscInt_FMT ".\n",p->restart_length));

  // Create function object
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNINVSQRT));
  PetscCall(FNSetFromOptions(fn));

  // Set up viewer
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(FNView(fn,viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,p->restart_length,p->restart_length,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));

  // compute invsqrt for non-symmetic A
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&myrand));
  PetscCall(PetscRandomSetFromOptions(myrand));
  PetscCall(PetscRandomSetInterval(myrand,0.1,1.0));
  PetscCall(MatDenseGetArray(A,&As));

  for (i=0;i<p->restart_length;i++) As[i+i*p->restart_length]=2.5;
  for (j=1;j<3;j++) {
    for (i=0;i<p->restart_length-j;i++) { As[i+(i+j)*p->restart_length]=1.0; As[(i+j)+i*p->restart_length]=1.0; }
  }
  for (j=1;j<3;j++) {
    for (i=0;i<p->restart_length-j;i++) As[(i+j)+i*p->restart_length]=0.0;
  }
  for (j=1;j<3;j++) {
    for (i=0;i<p->restart_length-j;i++) {
      PetscCall(PetscRandomGetValueReal(myrand,&v));
      As[(i+j)+i*p->restart_length]=v+2*v*PETSC_i;
    }
  }

  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(PetscRandomDestroy(&myrand));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(MatInvSqrt(fn,A,viewer,verbose,inplace));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());
}

int main( int argc, char **argv ) {

  int i;

  gmres_double_struct *p = &gp;
  p->restart_length = 10;

  int argcx=7;
  char **argvx = (char**) malloc( 7*sizeof(char*) );
  argvx[0] = (char*) malloc( 7*50*sizeof(char) );
  for( i=1;i<7;i++ ) {
    argvx[i] = argvx[0] + i*50;
  }

  // options in static form
  char str0[50] = "a.out";
  char str1[50] = "-fn_scale";
  char str2[50] = "1.0,1.0";
  char str3[50] = "-fn_method";
  char str4[50] = "0";
  char str5[50] = "-verbose";
  char str6[50] = "1";

  // options in dynamic form
  strcpy( argvx[0],str0 );
  strcpy( argvx[1],str2 );
  strcpy( argvx[2],str2 );
  strcpy( argvx[3],str3 );
  strcpy( argvx[4],str4 );
  strcpy( argvx[5],str5 );
  strcpy( argvx[6],str6 );

  small_dense_invsqrt( argcx, argvx, p );

  free( argvx[0] );
  free( argvx );

  return 0;
}

/*TEST

   testset:
      args: -fn_scale 0.9,0.5 -n 10
      filter: grep -v "computing matrix functions"
      requires: !__float128
      output_file: output/test8_1.out
      test:
         suffix: 1
         args: -fn_method {{0 1 2 3}}
      test:
         suffix: 1_cuda
         args: -fn_method 2 -matcuda
         requires: cuda
      test:
         suffix: 1_magma
         args: -fn_method {{1 3}} -matcuda
         requires: cuda magma
      test:
         suffix: 2
         args: -inplace -fn_method {{0 1 2 3}}
      test:
         suffix: 2_cuda
         args: -inplace -fn_method 2 -matcuda
         requires: cuda
      test:
         suffix: 2_magma
         args: -inplace -fn_method {{1 3}} -matcuda
         requires: cuda magma

TEST*/
