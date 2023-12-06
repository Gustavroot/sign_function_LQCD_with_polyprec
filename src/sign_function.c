/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include "main.h"



void export_matrix_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ){

  int i,j,start,end;
  vector_double out,in;

  PUBLIC_MALLOC( out, complex_double, l->vector_size );
  PUBLIC_MALLOC( in, complex_double, l->vector_size );

  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );

  FILE *fp;
  fp = fopen("../../../../../../../p/project/chwu29/ramirez1/matr_data.txt", "w");

  for ( i=0;i<p->v_end;i++ ) {

    START_MASTER(threading);
    printf0( "column : %d\n",i );
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);

    vector_double_define( in, 0, start, end, l );
    START_MASTER(threading);
    in[i] = 1.0;
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);

    apply_operator_double( out, in, p, l, threading );

    START_MASTER(threading);
    // saving data for the i-th column of D
    for ( j=0;j<p->v_end;j++ ) {
      if ( cabs_double(out[j]) > 1.0E-10 ) {
        fprintf(fp,"%d,%d,%.14f,%.14f\n",j,i,creal_double(out[j]),cimag_double(out[j]));
      }
    }
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);
  }

  fclose(fp);

  PUBLIC_FREE( out, complex_double, l->vector_size );
  PUBLIC_FREE( in, complex_double, l->vector_size );

  exit(0);
}


void check_rel_err( int j, gmres_double_struct *p ){

  //if ( j==10 ) {
  //  invsqrt_of_H( p->Hb1, p->H, j+1 );
  //}

  int check_fctr = 50;
  if ( g.use_polyprec==0 ) {
    //check_fctr = 20;
  } else {
    check_fctr /= p->polyprec_double->d_poly;
    check_fctr++;
  }

  if ( j>1 && j%check_fctr==0 ) {
    invsqrt_of_H( p->Hb1, p->H, j );
    invsqrt_of_H( p->Hb2, p->H, j+1 );

    // compute the indirect measure of the relative error from p->Hb1 and p->Hb2
    int ix;
    double normx1,normx2;
    complex_double diff_buff;
    normx1 = 0.0;
    for ( ix=0;ix<((j+1)-1);ix++ ) {
      diff_buff = p->Hb1[0][ix] - p->Hb2[0][ix];
      normx1 += diff_buff * conj_double(diff_buff);
    }
    diff_buff = 0.0 - p->Hb2[0][j];
    normx1 += diff_buff * conj_double(diff_buff);
    normx2 = 0.0;
    for ( ix=0;ix<(j+1);ix++ ) {
      diff_buff = p->Hb2[0][ix];
      normx2 += diff_buff * conj_double(diff_buff);
    }

    printf0( "indirect measure of relative error (m=%d) : %.8e\n",j+1,sqrt(normx1/normx2) );

    // save data for next iter
    // FIXME : should we just swap pointers here ?
    //memcpy( p->Hb1[0], p->Hb2[0], p->restart_length * p->restart_length );

    //if ( g.my_rank==0 && j==150 ) {
    //  FILE *fp1;
    //  FILE *fp2;
    //  fp1 = fopen("./matr_data1.txt", "w");
    //  fp2 = fopen("./matr_data2.txt", "w");
    //  int iy,jy;
    //  for ( iy=0;iy<(j+1);iy++ ) {
    //    for ( jy=0;jy<(j+1);jy++ ) {
    //      fprintf( fp1,"%d,%d,%.12f,%.12f\n",iy,jy,creal_double(p->H[jy][iy]),cimag_double(p->H[jy][iy]) );
    //      fprintf( fp2,"%d,%d,%.12f,%.12f\n",iy,jy,creal_double(p->Hb2[jy][iy]),cimag_double(p->Hb2[jy][iy]) );
    //    }
    //  }
    //  fclose(fp1);
    //  fclose(fp2);
    //  exit(0);
    //}
  }
}


int arnoldi_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ) {

  int m = p->restart_length;
  memset( p->H[0], 0, (m*(m+1))*sizeof(complex_double) );

  // start and end indices for vector functions depending on thread
  int start;
  int end;

  int j=-1, finish=0, iter=0, il, ol;
  complex_double gamma0 = 0;

  double norm_r0=1, gamma_jp1=1, t0=0, t1=0;
  START_LOCKED_MASTER(threading)
  if ( l->depth==0 && ( p->timing || p->print ) ) prof_init( l );
  if ( l->depth == 0 ) t0 = MPI_Wtime();
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  for( ol=0; ol<p->num_restart && finish==0; ol++ )  {

    if ( p->polyprec_double->apply==1 ) {
      sign_function_prec_pow2( p->r, p->b, p, l, threading );
    } else {
      vector_double_copy( p->r, p->b, start, end, l );
    }

    gamma0 = (complex_double) global_norm_double( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
    START_MASTER(threading)
    p->gamma[0] = gamma0;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if( ol == 0) {
      norm_r0 = creal(p->gamma[0]);
    }
    
    vector_double_real_scale( p->V[0], p->r, 1/p->gamma[0], start, end, l ); // v_0 = r / gamma_0
    
    for( il=0; il<p->restart_length && finish==0; il++) {
      j = il; iter++;

      // one step of Arnoldi
      if ( !arnoldi_step_double( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p->shift, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
        break;
      }
      
      //if ( cabs( p->H[j][j+1] ) > p->tol/10 ) {
      //  qr_update_double( p->H, p->s, p->c, p->gamma, j, l, threading );
      //  gamma_jp1 = cabs( p->gamma[j+1] );
      //  
      //  if( gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5 ) { // if satisfied ... stop
      //    finish = 1;
      //    START_MASTER(threading)
      //      if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_double, iter = %d, level=%d\n", iter, l->level );
      //    END_MASTER(threading)
      //  }
      //} else {
      //  printf0("depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", l->depth, iter, j+1, j, CSPLIT( p->H[j][j+1] ) );
      //  finish = 1;
      //  break;
      //}

      // after a certain number of Arnoldi steps, compute an indirect measure of the error
      check_rel_err( j, p );

    } // end of a single restart

    //compute_solution_double( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
    //                            p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );

  } // end of fgmres

  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.total_time = t1-t0; g.iter_count = iter; g.norm_res = gamma_jp1/norm_r0; }
  END_LOCKED_MASTER(threading)
  
  if ( p->print ) {
    printf0( "\nsome specific timings from Arnoldi :\n\n" );
  //  beta = gamma_jp1;
    START_MASTER(threading)
  //  g.norm_res = creal(beta)/norm_r0;
    printf0("+----------------------------------------------------------+\n");
  //  printf0("|       FGMRES iterations: %-6d coarse average: %-6.2lf   |\n", iter,
  //          ((double)g.coarse_iter_count)/((double)iter) );
  //  printf0("| exact relative residual: ||r||/||b|| = %e      |\n", creal(beta)/norm_r0 );
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", t1-t0 );
  //  if ( g.coarse_time > 0 ) 
  //    printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
  //            g.coarse_time, 100*(g.coarse_time/(t1-t0)) );
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  
  if ( l->level == 0 ) {
    START_LOCKED_MASTER(threading)
    g.coarse_iter_count += iter;
    END_LOCKED_MASTER(threading)
  }
    
  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", iter );
      printf0("solve time: %le seconds\n", t1-t0 );
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
    END_LOCKED_MASTER(threading)
    }
  }

  if ( l->depth == 0 && ( p->timing || p->print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
    if ( g.method != 6 ) prof_print( l );
    END_MASTER(threading)
  }

  return iter;
}


void timings_probe_matvec_vs_dotprod( level_struct* l, struct Thread* threading ){

  printf0(GREEN"\n--------------------------------------------------------\n");
  printf0("********************* CHECK #1 *************************\n");
  printf0("--------------------------------------------------------\n\n"RESET);

  //printf0( "sign function computation under construction ... timing matmuls and dot products first (average over 100 runs)\n\n" );

  // for now, let us measure some times for matmul and dot products

  int i, nr_calls=100;
  double t1=0, t2=0;

  for ( i=0; i<nr_calls; i++ ) {
    t1 -= MPI_Wtime();
    apply_operator_double( g.p.w, g.p.x, &(g.p), l, threading );
    t1 += MPI_Wtime();
  }

  for ( i=0; i<nr_calls; i++ ) {
    t2 -= MPI_Wtime();
    global_inner_product_double( g.p.w, g.p.x, g.p.v_start, g.p.v_end, l, threading );
    t2 += MPI_Wtime();
  }

  //printf0( "avg matmul time : %lf seconds \n", t1/nr_calls );
  //printf0( "avg dotprod time : %lf seconds \n", t2/nr_calls );

  printf0( "ratio of matvecs/dotprods on your machine : %f\n", t1/t2 );

  printf0( "\n" );
}


void check_arnoldi_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ){

  printf0(GREEN"\n--------------------------------------------------------\n");
  printf0("********************* CHECK #2 *************************\n");
  printf0("--------------------------------------------------------\n\n"RESET);

  //printf0( "\nCHECK #2\n\n" );
  printf0( "running Arnoldi, checking the Arnoldi relation\n\n" );

  // call a whole Arnoldi
  arnoldi_double( p, l, threading );

  // print the Hessenberg matrix
  //int i,j;
  //for ( i=0;i<p->restart_length+1;i++ ) {
  //  for ( j=0;j<p->restart_length;j++ ) {
  //    printf0( "%f+i%f\t", CSPLIT(p->H[j][i]) );
  // }
  //  printf0( "\n" );
  //}

  int i,j,start,end;
  double frob_norm = 0.0;
  vector_double v1 = p->x;
  vector_double v2 = p->w;
  vector_double v3 = p->r;

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  for ( i=0;i<p->restart_length;i++ ) {

    // first, compute the vector from the LHS of the Arnoldi relation

    apply_operator_double( p->wy, p->V[i], p, l, threading );
    apply_operator_double( p->wx, p->wy, p, l, threading );
    if ( g.global_shift != 0.0 ) {
      int startx, endx;
      compute_core_start_end_custom( p->v_start, p->v_end, &startx, &endx, l, threading, l->num_lattice_site_var );
      vector_double_saxpy( p->wx, p->wx, p->V[i], g.global_shift, startx, endx, l );
    }
    sign_function_prec_pow2( p->wy, p->wx, p, l, threading );
    sign_function_prec_pow2( v1, p->wy, p, l, threading );

//    sign_function_prec_pow2( p->wy, p->V[i], p, l, threading );
//    sign_function_prec_pow2( p->wx, p->wy, p, l, threading );
//    apply_operator_double( p->w, p->wx, p, l, threading );
//    apply_operator_double( v1, p->w, p, l, threading );

    // second, the RHS one
    vector_double_define( v2, 0, start, end, l );
    for ( j=0;j<p->restart_length+1;j++ ) {
      vector_double_saxpy( v2, v2, p->V[j], p->H[i][j], start, end, l );
    }

    vector_double_minus( v3, v1, v2, start, end, l );

    frob_norm += global_inner_product_double( v3, v3, p->v_start, p->v_end, l, threading );
  }

  printf0( "relative error in Arnoldi relation :%.14e\n\n", frob_norm/(12*l->inner_vector_size*g.num_processes) );

  // check orthonormality of Arnoldi basis
  check_orthonormality( p->V, p, l, threading );
}


void sign_function_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ) {

  int i,j,k,start,end;
  double t0,t1;

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  printf0(GREEN"\n--------------------------------------------------------\n");
  printf0("************** COMPUTING SIGN FUNCTION *****************\n");
  printf0("--------------------------------------------------------\n\n"RESET);

  //printf0( "\nCOMPUTING THE SIGN FUNCTION NOW\n\n" );

  // call the Arnoldi relation on the preconditioned system
  t0 = MPI_Wtime();
  arnoldi_double( p, l, threading );
  t1 = MPI_Wtime();
  printf0( "total time spent on Arnoldi : %f seconds\n\n", t1-t0 );

  // create the first unit vector and buffer vector
  complex_double* e1 = NULL;
  complex_double* b1 = NULL;
  MALLOC( e1, complex_double, p->restart_length );
  MALLOC( b1, complex_double, p->restart_length );

  // create buffer small matrix for storing H^{-1/2}
  complex_double** His = NULL;
  MALLOC( His, complex_double*, p->restart_length );
  His[0] = NULL;
  MALLOC( His[0], complex_double, p->restart_length*p->restart_length );
  for ( i=1;i<p->restart_length;i++ ) {
    His[i] = His[0] + i*p->restart_length;
  }

  //// and create another buffer to store H^{2}
  //complex_double** Hsq = NULL;
  //MALLOC( Hsq, complex_double*, p->restart_length );
  //Hsq[0] = NULL;
  //MALLOC( Hsq[0], complex_double, p->restart_length*p->restart_length );
  //for ( i=1;i<p->restart_length;i++ ) {
  //  Hsq[i] = Hsq[0] + i*p->restart_length;
  //}

  //// before computing (H^{2})^{-1/2}, compute H^{2}
  //for( i=0;i<p->restart_length;i++ ) {
  //  for( j=0;j<p->restart_length;j++ ) {
  //    Hsq[j][i] = 0.0;
  //    for( k=0;k<p->restart_length;k++ ) {
  //      Hsq[j][i] += p->H[k][i]*p->H[j][k];
  //    }
  //  }
  //}

  // compute H^{-1/2}
  t0 = MPI_Wtime();
  START_MASTER(threading)
  invsqrt_of_H( His, p->H, p->restart_length );
  //invsqrt_of_H( His, Hsq, p->restart_length );
  END_MASTER(threading)
  t1 = MPI_Wtime();
  printf0( "\ntime spent on invsqrt_of_H : %f\n", t1-t0 );

  // do (H^{2})^{-1/2}*e1
  memcpy( b1, His[0], p->restart_length*sizeof(complex_double) );

  // do Vm*b1 and store this in p->w
  t0 = MPI_Wtime();
  vector_double_define( p->w, 0, start, end, l );
  for ( i=0;i<p->restart_length;i++ ) {
    vector_double_saxpy( p->w, p->w, p->V[i], b1[i], start, end, l );
  }
  t1 = MPI_Wtime();
  printf0( "\ntime spent on Vm*b1 : %f\n\n", t1-t0 );

  apply_operator_double( p->x, p->w, p, l, threading ); // x = D*w

  // finally, scale with p->gamma[0]
  vector_double_scale( p->x, p->x, p->gamma[0], start, end, l );

  FREE( e1, complex_double, p->restart_length );
  FREE( b1, complex_double, p->restart_length );
  FREE( His[0], complex_double, p->restart_length*p->restart_length );
  FREE( His, complex_double*, p->restart_length );
}


void sign_function_prec_pow1( vector_double out, vector_double in, gmres_double_struct* p, level_struct* l, struct Thread* threading ) {

  // UNDER CONSTRUCTION

  int start,end;
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  vector_double_copy( out, in, start, end, l );
}


void sign_function_prec_pow2( vector_double out, vector_double in, gmres_double_struct* p, level_struct* l, struct Thread* threading ) {

  if ( g.use_polyprec==0 ) {
    int start,end;
    compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
    vector_double_copy( out, in, start, end, l );
  } else {
    apply_polyprec_double( out, NULL, in, 0, l, threading );
  }
}


void invsqrt_of_H( complex_double** His, complex_double** H, int n ) {

  small_dense_invsqrt( His, H, n );

  //// TODO : implement this correctly
  //int i,j;
  //for ( i=0;i<n;i++ ) {
  //  for ( j=0;j<n;j++ ) {
  //    His[j][i] = (double) (((double)rand()/(double)RAND_MAX));
  //  }
  //}
}


void check_orthonormality( vector_double* V, gmres_double_struct* p, level_struct* l, struct Thread* threading ) {

  int i,j;

  // variable to store the relative error
  double frob_norm = 0.0;

  complex_double buff;

  // buffer for dot products
  complex_double* dotprods1 = NULL;
  complex_double* dotprods2 = NULL;
  MALLOC( dotprods1, complex_double, p->restart_length+1 );
  MALLOC( dotprods2, complex_double, p->restart_length+1 );

  for ( i=0;i<p->restart_length+1;i++ ) {

    process_multi_inner_product_double( p->restart_length+1, dotprods1, V, V[i], 0, l->inner_vector_size, l, threading );

    MPI_Allreduce( dotprods1, dotprods2, p->restart_length+1, MPI_COMPLEX_double, MPI_SUM, g.comm_cart );

    for ( j=0;j<p->restart_length+1;j++ ) {
      if ( j==i ) {
        buff = dotprods2[j] - 1.0;
      }
      else {
        buff = dotprods2[j];
      }
      frob_norm += buff*conj(buff);
    }
  }

  printf0( "relative error in orthonormality : %.14e\n\n", frob_norm/(p->restart_length+1) );

  FREE( dotprods1, complex_double, p->restart_length+1 );
  FREE( dotprods2, complex_double, p->restart_length+1 );
}
