/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include "main.h"



#ifdef USE_SLEPC

#include <slepcfn.h>
#include <petscsys.h>
#include <complex.h>
#include <math.h>
#include <slepceps.h>



static char help[] = "Test matrix inverse square root.\n\n";

/*
   Compute matrix inverse square root B = inv(sqrtm(A))
   Check result as norm(B*B*A-I)
*/

PetscErrorCode MatInvSqrt(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace,complex_double** His)
{
  PetscScalar    tau,eta,*Sx;
  //PetscComplex   *Sx;
  PetscReal      nrm;
  //PetscBool      set,flg;
  PetscInt       n;
  Mat            S,R,Acopy;
  Vec            v,f0;
  int            i,j;

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
    //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n"));
    //PetscCall(MatView(A,viewer));
    //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed inv(sqrtm(A)) - - - - - - -\n"));
    //PetscCall(MatView(S,viewer));
  }

  // finally, set His from S
  PetscCall(MatDenseGetArray(S,&Sx));
  // set His from Sx
  for( j=0;j<n;j++ ) {
    for( i=0;i<n;i++ ) {
      His[j][i] = Sx[i+j*n];
    }
  }
  PetscCall(MatDenseRestoreArray(S,&Sx));

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
  //if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S*A-I||_F < 100*eps\n"));
  //else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S*A-I||_F = %g\n",(double)nrm));
  // check FNEvaluateFunctionMatVec()
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&S));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&f0));

  return 0;
}

PetscErrorCode _small_dense_invsqrt( int argcx,char **argvx, complex_double **His, complex_double **H, int n )
{
  FN             fn;
  Mat            A=NULL;
  PetscScalar    *As;
  //PetscComplex   *As;
  PetscViewer    viewer;
  PetscInt       i,j;
  PetscBool      verbose,inplace;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argcx,&argvx,(char*)0,help));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix inverse square root, n=%" PetscInt_FMT ".\n",n));

  // Create function object
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNINVSQRT));
  PetscCall(FNSetFromOptions(fn));

  // Set up viewer
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  //PetscCall(FNView(fn,viewer));
  //if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));

  // compute invsqrt for non-symmetic A, in our case will be the Hessenberg matrix

  PetscCall(MatDenseGetArray(A,&As));

  // IMPORTANT : j is columns, i is rows. As is assumed column major, and
  //             the same for H

  // set As from H
  for( j=0;j<n;j++ ) {
    for( i=0;i<n;i++ ) {
      As[i+j*n] = H[j][i];
    }
  }

  PetscCall(MatDenseRestoreArray(A,&As));

  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(MatInvSqrt(fn,A,viewer,verbose,inplace,His));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());

  return 0;
}

void small_dense_invsqrt( complex_double **His, complex_double **H, int n ) {

  int i,j,k;
  //PetscRandom    myrand;
  //double         v;

  //complex_double **His;
  //gmres_double_struct *p = &(g.p);

  // TODO : switch to MALLOC and FREE macros

  // allocate His (will contain the invsqrt solution)
  //His = (complex_double**) malloc( p->restart_length*sizeof(complex_double*) );
  //His[0] = (complex_double*) malloc( p->restart_length*p->restart_length*sizeof(complex_double) );
  //for( i=1;i<p->restart_length;i++ ) {
  //  His[i] = His[0] + i*p->restart_length;
  //}

  // mimic input via command line
  int argcx=7;
  char **argvx = (char**) malloc( 7*sizeof(char*) );
  argvx[0] = (char*) malloc( 7*50*sizeof(char) );
  for( i=1;i<7;i++ ) {
    argvx[i] = argvx[0] + i*50;
  }

  // 'input' options in static form
  char str0[50] = "a.out";
  char str1[50] = "-fn_scale";
  char str2[50] = "1.0,1.0";
  char str3[50] = "-fn_method";
  char str4[50] = "0";
  char str5[50] = "-verbose";
  char str6[50] = "0";

  // 'input' options in dynamic form
  strcpy( argvx[0],str0 );
  strcpy( argvx[1],str1 );
  strcpy( argvx[2],str2 );
  strcpy( argvx[3],str3 );
  strcpy( argvx[4],str4 );
  strcpy( argvx[5],str5 );
  strcpy( argvx[6],str6 );

  _small_dense_invsqrt( argcx, argvx, His, H, n );

  // check that : H * His * His = I

  // buffers for the multiplication
  complex_double **B;
  complex_double **C;
  B = (complex_double**) malloc( n*sizeof(complex_double*) );
  C = (complex_double**) malloc( n*sizeof(complex_double*) );
  B[0] = (complex_double*) malloc( n*n*sizeof(complex_double) );
  C[0] = (complex_double*) malloc( n*n*sizeof(complex_double) );
  for( i=1;i<n;i++ ) {
    B[i] = B[0] + i*n;
  }
  for( i=1;i<n;i++ ) {
    C[i] = C[0] + i*n;
  }

  // first : B = His*His

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      B[j][i] = 0.0;
      for( k=0;k<n;k++ ) {
        B[j][i] += His[k][i]*His[j][k];
      }
    }
  }

  // second : C = H*B

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      C[j][i] = 0.0;
      for( k=0;k<n;k++ ) {
        C[j][i] += H[k][i]*B[j][k];
      }
    }
  }

  // compute Frobenius norm of identity
  double frob_norm_id = sqrt((double)n);

  // compute Frobenius norm of C-I

  double frob_norm_CminI = 0.0;
  complex_double buffx;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( i==j ) {
        buffx = C[j][i] - 1.0;
      }
      else {
        buffx = C[j][i];
      }
      buffx = buffx*conj(buffx);
      frob_norm_CminI += creal(buffx);
    }
  }
  frob_norm_CminI = sqrt(frob_norm_CminI);

  printf0( "relative error from invsqrt : %.8e\n",frob_norm_CminI/frob_norm_id );

  free( argvx[0] );
  free( argvx );
  //free( p->H[0] );
  //free( His[0] );
  free( B[0] );
  free( C[0] );
  //free( p->H );
  //free( His );
  free( B );
  free( C );

  //return His;
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

PetscErrorCode eig_op( Mat A, Vec x, Vec y )
{
  g.eig_ctr++;
  printf0( "Call to eig_op number %d\n",g.eig_ctr );

  void              *ctx;
  int                lo,i,j;
  const PetscScalar *px;
  PetscScalar       *py;
  level_struct *l;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  l = (level_struct*)ctx;

  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));

  int start,end;
  gmres_double_struct *p = &(g.p);
  start = p->v_start;
  end = p->v_end;

  // this mimics the identity matrix
  //vector_double_copy( py, px, start, end, l );

  //apply_operator_double( py, px, p, l, g.threading ); // p->wy

  apply_operator_double( p->wy, px, p, l, g.threading ); // wx = D*wy
  apply_operator_double( p->wx, p->wy, p, l, g.threading ); // w = D*wx
  sign_function_prec_pow2( p->wy, p->wx, p, l, g.threading ); // Z[j] = q(D)*V[j]
  sign_function_prec_pow2( py, p->wy, p, l, g.threading ); // wy = q(D)*Z[j]

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int _eig_via_slepc( int argc, char **argv, level_struct *l ){

  Vec            v0;              /* initial vector */
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscInt       N,nev;
  level_struct   *lx = l;
  PetscMPIInt    rank;
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  // create the shell for the operator

  gmres_double_struct *p = &(g.p);

  int nx = l->inner_vector_size;
  int Nx = nx*g.num_processes;
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,Nx,Nx,lx,&A));
  PetscCall(MatShellSetOperation(A,MATOP_MULT,(void(*)(void))eig_op));

  // create eigensolver context
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  // set operators. In this case, it is a standard eigenvalue problem
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));

  // set solver parameters at runtime
  PetscCall(EPSSetFromOptions(eps));

  // solve the eigensystem
  PetscCall(EPSSolve(eps));

  // optional: Get some information from the solver and display it
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  // display solution and clean up

  // show detailed info unless -terse option is given by user
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  // free data
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  //PetscCall(VecDestroy(&v0));
  PetscCall(SlepcFinalize());

  return 0;
}

void eig_via_slepc( level_struct *l ){
  int i;

  // mimic input via command line
  int argcx=16;
  int chr_lngth=100;
  char **argvx = (char**) malloc( argcx*sizeof(char*) );
  argvx[0] = (char*) malloc( argcx*chr_lngth*sizeof(char) );
  for( i=1;i<argcx;i++ ) {
    argvx[i] = argvx[0] + i*chr_lngth;
  }

  // 'input' options in static form
  char str0[100] = "a.out";
  char str1[100] = "-eps_smallest_magnitude";
  char str2[100] = "-eps_nev";
  char str3[100] = "5";
  char str4[100] = "-eps_two_sided";
  char str5[100] = "0";
  char str6[100] = "-eps_krylovschur_locking";
  char str7[100] = "0";
  char str8[100] = "-ds_parallel";
  char str9[100] = "distributed";
  char str10[100] = "-eps_ncv";
  char str11[100] = "50";
  char str12[100] = "-eps_tol";
  char str13[100] = "1.0e-9";
  char str14[100] = "-eps_max_it";
  char str15[100] = "100";

  // 'input' options in dynamic form
  strcpy( argvx[0],str0 );
  strcpy( argvx[1],str1 );
  strcpy( argvx[2],str2 );
  strcpy( argvx[3],str3 );
  strcpy( argvx[4],str4 );
  strcpy( argvx[5],str5 );
  strcpy( argvx[6],str6 );
  strcpy( argvx[7],str7 );
  strcpy( argvx[8],str8 );
  strcpy( argvx[9],str9 );
  strcpy( argvx[10],str10 );
  strcpy( argvx[11],str11 );
  strcpy( argvx[12],str12 );
  strcpy( argvx[13],str13 );
  strcpy( argvx[14],str14 );
  strcpy( argvx[15],str15 );

  _eig_via_slepc( argcx, argvx, l );

  free( argvx[0] );
  free( argvx );
}

#endif
