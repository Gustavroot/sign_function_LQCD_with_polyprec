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

void rhs_define( vector_double rhs, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading here
  if(threading->thread != 0)
    return;

  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];

  if ( g.rhs == 0 ) {
    vector_double_define( rhs, 1, start, end, l );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = ones\n");
    END_MASTER(threading)
  } else if ( g.rhs == 1 )  {
    vector_double_define( rhs, 0, start, end, l );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      rhs[0] = 1.0;
      END_LOCKED_MASTER(threading)
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
    END_MASTER(threading)
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = random\n");
    END_MASTER(threading)
  } else if ( g.rhs == 3 ) {
    vector_double_define( rhs, 0, start, end, l );
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
}


int wilson_driver( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];
  
  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
  vector_double sol = g.mixed_precision==2?g.p_MP.dp.x:g.p.x;

#ifdef WILSON_BENCHMARK
  START_MASTER(threading)
  prof_init( l );
  END_MASTER(threading)
  double t = -MPI_Wtime();
  double t_min = 1000;
  for ( int i=0; i<100; i++ ) {
    double tmp_t = -MPI_Wtime();
#endif

  vector_double_copy( rhs, source, start, end, l );

  int save_solution_to_file = g.save_solution_to_file; //0;
  int check_sign2_eq_Id = g.check_sign2_eq_Id; //1;
  //char save_solution_filename[] = g.; //"sign_invsqrt_solution.dat";
  char save_solution_filename[500]; strcpy( save_solution_filename,g.save_solution_filename );

  if ( g.check_with_large_vecs==1 ) {
    START_MASTER(threading)
    lime_read_vector( (double*)g.p.invsqrt_sol,save_solution_filename );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
  }

  // simple timings test
  timings_probe_matvec_vs_dotprod( l, threading );

#ifdef POLYPREC
  // construct the polynomial preconditioner and test it
  if ( g.use_polyprec==1 ) {
    int buffz = g.check_do;
    START_MASTER(threading)
    g.check_do = 0;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    set_up_polynomial_and_test_double( &(g.p), l, threading );

    START_MASTER(threading)
    g.check_do = buffz;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
  }
#endif

  // calling a simple Arnoldi
  int buffx = g.p.restart_length;
  int buffy = g.check_do;
  START_MASTER(threading)
  g.p.restart_length = 10;
  g.check_do = 0;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  check_arnoldi_double( &(g.p), l, threading );
  START_MASTER(threading)
  g.p.restart_length = buffx;
  g.check_do = buffy;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // print matrix out, to do checks in MATLAB
  //export_matrix_double( &(g.p), l, threading );

  if ( g.read_rhs_from_file==1 ) {
    START_MASTER(threading)
    lime_read_vector( (double*)g.p.b,g.read_rhs_filename );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
  }

  // computing the sign function
  sign_function_double( &(g.p), l, threading );

  if ( save_solution_to_file==1 ) {
    // save g.p.x
    START_MASTER(threading)
    lime_write_vector( (double*)g.p.x,save_solution_filename );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    //START_MASTER(threading)
    //lime_read_vector( (double*)g.p.r,save_solution_filename );
    //END_MASTER(threading)
    //SYNC_MASTER_TO_ALL(threading)

    //vector_double_minus( g.p.r, g.p.r, g.p.x, start, end, l );
    //double normbff = global_norm_double( g.p.r, g.p.v_start, g.p.v_end, l, threading );
    //START_MASTER(threading)
    //printf0( "diff vecs = %.12f\n",normbff );
    //END_MASTER(threading)
  }

  if ( check_sign2_eq_Id==1 ) {
    // and repeat, to verify that sign^2=1

    // for this, first back up the RHS from the previous "solve"
    vector_double_copy( g.p.wz, g.p.b, start, end, l );
    // then, set the RHS to the previous solution
    vector_double_copy( g.p.b, g.p.x, start, end, l );
    // then, "solve"
    sign_function_double( &(g.p), l, threading );
    // after that sign function application, x should be "similar" to wz
    double norm_bef = global_norm_double( g.p.wz, g.p.v_start, g.p.v_end, l, threading );
    vector_double_minus( g.p.wz, g.p.wz, g.p.x, start, end, l );
    double norm_aft = global_norm_double( g.p.wz, g.p.v_start, g.p.v_end, l, threading );
    START_MASTER(threading)
    printf0( "relative error in sign^2=1 : %.8e\n",norm_aft/norm_bef );
    END_MASTER(threading)
  }

  /*

  if ( g.method == -1 ) {
    cgn_double( &(g.p), l, threading );
  } else if ( g.mixed_precision == 2 ) {
    iter = fgmres_MP( &(g.p_MP), l, threading );
  } else {
    iter = fgmres_double( &(g.p), l, threading );
  }

  */

  vector_double_copy( solution, sol, start, end, l );
#ifdef WILSON_BENCHMARK
    tmp_t += MPI_Wtime();
    if ( tmp_t < t_min )
      t_min = tmp_t;
  }
  t +=MPI_Wtime();
  START_MASTER(threading)
  printf0("average over 100 solves: %lf seconds\n", t/100 );
  printf0("minimum out of 100 solves: %lf seconds\n", t_min );
  prof_print( l );
  END_MASTER(threading)
#endif
  
  return iter;
}


void solve( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;

  if ( g.vt.evaluation ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    scan_var( &(g.vt), l );
    END_LOCKED_MASTER(threading)
  } else {
    wilson_driver( solution, source, l, threading );
  }
}


void solve_driver( level_struct *l, struct Thread *threading ) {
  
  vector_double solution = NULL, source = NULL;
  
  START_MASTER(threading)
  MALLOC( solution, complex_double, l->inner_vector_size );
  MALLOC( source, complex_double, l->inner_vector_size );
  // use threading->workspace to distribute pointer to newly allocated memory to all threads
  ((vector_double *)threading->workspace)[0] = solution;
  ((vector_double *)threading->workspace)[1] = source;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  solution = ((vector_double *)threading->workspace)[0];
  source   = ((vector_double *)threading->workspace)[1];
  
  rhs_define( source, l, threading );
  
  solve( solution, source, l, threading );

  START_LOCKED_MASTER(threading)
  FREE( solution, complex_double, l->inner_vector_size );
  FREE( source, complex_double, l->inner_vector_size );
  END_LOCKED_MASTER(threading)
}

