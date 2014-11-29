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


#include "impact_functions.h"
#include "pareto2.h"

static const R_CMethodDef cMethods[] = {
	{"index_h", (DL_FUNC)&index_h, 3},
	{"index_h_log", (DL_FUNC)&index_h_log, 3},
	{"index_g", (DL_FUNC)&index_g, 3},
	{"Sstat2", (DL_FUNC)&Sstat2, 3},
	{"index_rp_finite", (DL_FUNC)&index_rp_finite, 4},
	{"index_lp_finite", (DL_FUNC)&index_lp_finite, 5},
	{"index_rp_infinite", (DL_FUNC)&index_rp_infinite, 3},
	{"index_lp_infinite", (DL_FUNC)&index_lp_infinite, 3},
	{"pareto2_phirsch", (DL_FUNC)&pareto2_phirsch, 5},
	{"pareto2_dhirsch", (DL_FUNC)&pareto2_dhirsch, 5},
	{NULL, NULL, 0}
};


void R_init_CITAN(DllInfo *dll)
{
	R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);

#ifdef CITAN_DEBUG
	fprintf(stderr, "DEBUG: Dynamic library CITAN loaded.\n");
#endif
}

