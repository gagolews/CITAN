/* ************************************************************************* *
 *   This file is part of the CITAN library.                                 *
 *                                                                           *
 *   Copyright 2011 Marek Gagolewski                                         *
 *                                                                           *
 *                                                                           *
 *   CITAN is free software: you can redistribute it and/or modify           *
 *   it under the terms of the GNU Lesser General Public License             *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   CITAN is distributed in the hope that it will be useful,                *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the             *
 *   GNU Lesser General Public License for more details.                     *
 *                                                                           *
 *   You should have received a copy of the GNU Lesser General Public        *
 *   License along with CITAN. If not, see <http://www.gnu.org/licenses/>.   *
 * ************************************************************************* */




#ifndef __pareto2_h
#define __pareto2_h


/* #define CITAN_DEBUG */


#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>




void pareto2_phirsch(double* x, int* m, double* n, double* k, double* s);
void pareto2_dhirsch(double* x, int* m, double* n, double* k, double* s);



#endif
