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

#ifndef SIGN_FUNCTION_HEADER
  #define SIGN_FUNCTION_HEADER

  struct Thread;

  // this function does a subset of the operations of fgmres_double(...)
  int arnoldi_double( gmres_double_struct *p, level_struct *l, struct Thread *threading );

  void check_arnoldi_double( gmres_double_struct *p, level_struct *l, struct Thread *threading );
  void timings_probe_matvec_vs_dotprod( level_struct* l, struct Thread* threading );

#endif
