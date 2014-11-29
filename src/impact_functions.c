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



/** Function to compute the h-index, O(n) time.
 *  @param x vector of non-negative reals, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param out pointer to the result (return value)
 */
void index_h(double* x, int* n, double* out)
{
	int i = 0;
	int k = *n;
	while (i<k)
	{
		if (x[i] < (double)i+1) break;
		++i;
	}
	*out = (double)i;
}


/** *internal* function to compute the h-index, O(log n) time.
 *  @param x vector of non-negative reals, sorted non-increasingly
 *  @param n number of observations, n >= 1
 *  @return result
 */
int __index_h_log(double* x, int n)
{
	int h1 = 0;
	int h2 = n-1;

	if (x[0] < 1.0) { return 0; }

	while (1)
	{
		int m = (h2+h1+1)/2;
		if (x[m] == (double)m+1 || h1 == h2) {return m+1; }
		if (x[m] < (double)m+1) h2 = m-1;
		else h1 = m;
	}
}


/** Function to compute the h-index, O(log n) time.
 *  @param x vector of non-negative reals, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param out pointer to the result (return value)
 */
void index_h_log(double* x, int* n, double* out)
{
	*out = (double) __index_h_log(x, *n);
}


/** Function to compute the g-index, O(n) time.
 *  @param x vector of non-negative reals, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param out pointer to the result (return value)
 */
void index_g(double* x, int* n, double* out)
{
	int i = 0;
	int k = *n;
	double sum = 0.0;
	while (i<k)
	{
		sum += x[i];
		if (sum < (double)((i+1)*(i+1))) break;
		++i;
	}
	*out = (double)i;
}


/** Function to compute the r_p-index, O(n) time.
 *  @param x vector of non-negative reals, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param p pointer to the index order, 1 <= p < oo
 *  @param out pointer to the result (return value)
 */
void index_rp_finite(double* x, int* n, double *p, double* out)
{
	int N = *n;
	double P = *p;
	double r2p = pow((double)N, P);
	int i;

	if (x[0] <= 0.0) { *out = 0.0; return; }

	for (i=0; i<N; ++i)
	{
		double xip = pow(x[i],P);
		double ip  = pow((double)i,P);
		if (r2p-ip > xip)
			r2p = ip + xip;
	}

	*out = pow(r2p, 1.0/P);
}



/** Function to compute the r_oo-index, O(log n) time.
 *  @param x vector of non-negative numbers, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param out pointer to the result (return value)
 */
void index_rp_infinite(double* x, int* n, double* out)
{
	int N = *n;
	int h =__index_h_log(x, N);

	if (h < N)
	{
		if ((double)h < x[h])
			*out = x[h];
		else
			*out = (double)h;
	}
	else
		*out = (double)h;
}



/** *internal* function for index_lp_finite
 *  check if L^p ellipse interpolating (ui,vi) and (uj,vj) contains (uk,vk)
 *  @param a2p return value #1
 *  @param b2p return value #2
 */
void __index_lp_finite_getAB(double p, double ui, double vi, double uj, double vj, double* a2p, double* b2p)
{
	double uip = pow(ui,p);
	double ujp = pow(uj,p);
	double vip = pow(vi,p);
	double vjp = pow(vj,p);
	double c = uip*vjp-ujp*vip;


#ifdef CITAN_DEBUG
	if (c < 1e-15 && c > 1e-15)
		fprintf(stderr, "CITAN_DEBUG: __index_lp_finite_getAB: c==0\n");
#endif

	*a2p =  c/(vjp-vip);
	*b2p = -c/(ujp-uip);
}



/** *internal* function for index_lp_finite
 *  check if L^p ellipse interpolating (ui,vi) and (uj,vj) contains (uk,vk)
 *  @return 1 if true or 0 otherwise
 */
int __index_lp_finite_testContains(double uk, double vk, double p, double ui, double vi, double uj, double vj)
{
	double a2p, b2p;
	__index_lp_finite_getAB(p,ui,vi,uj,vj,&a2p,&b2p);

#ifdef CITAN_DEBUG
	if (a2p < 0.0)
		fprintf(stderr, "CITAN_DEBUG: __index_lp_finite_testContains: a2p<0\n");
	else if (a2p < 1e-15)
		fprintf(stderr, "CITAN_DEBUG: __index_lp_finite_testContains: a2p==0\n");
#endif

	if (b2p*(1.0 - pow(uk,p)/a2p) >= pow(vk,p))
		return 1;
	else
		return 0;
}


/** Function to compute the l_p-index, O(n) time.
 *  The procedure bases on Graham's scan for determining the convex hull
 *  of a planar set, see (Gagolewski,Debski,Nowakiewicz,2009b).
 *  @param x vector of non-negative numbers, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param p pointer to the index order, 1 <= p < oo
 *  @param s array of length *n+1, used as stack
 *  @param out pointer to the result (return value)
 */
void index_lp_finite(double* x, int* n, double *p, int* s, double* out)
{
	double P = *p;
	double N = *n;
	int m = 0; /* stack s is empty */
	int i = 1;
	int j; /*, k */
	double a2p, b2p;

#define LP_EPS 1e-9

	if (x[0] <= 0.0) { out[0] = 0.0; out[1] = 0.0; return; }

	while (i<N && x[i] >= x[0]-LP_EPS) ++i;

	s[m++] = 0; /* push 0 */
	s[m++] = i; /* push i */

	for (j=i+1; j<=N; ++j)
	{
		double vj = (j<N)?x[j]:0.0;
		if (vj >= x[s[m-1]]-LP_EPS) continue;

		while (m >= 2 && __index_lp_finite_testContains((double)s[m-2],x[s[m-2]], P, (double)s[m-1],x[s[m-1]], (double)j,vj))
		{
			--m; /* pop */
		}
		s[m++] = j; /* push j */
	}

	/*  now, selectMaxPair()  */
/*  	k = 0; */
	__index_lp_finite_getAB(P, (double)s[0],x[s[0]], (double)s[1],((s[1]<N)?x[s[1]]:0.0), &a2p, &b2p);
	for (j=1; j<m-1; ++j)
	{
		double a22p, b22p;
		__index_lp_finite_getAB(P, (double)s[j],x[s[j]], (double)s[j+1],((s[j+1]<N)?x[s[j+1]]:0.0), &a22p, &b22p);

#ifdef CITAN_DEBUG
		if (a22p < 1e-15 || b22p < 1e-15)
			fprintf(stderr, "CITAN_DEBUG: index_lp_finite: a22p < 1e-15 || b22p < 1e-15\n");
#endif

		if (a2p*b2p < a22p*b22p)
		{
			a2p = a22p;
			b2p = b22p;
/*			k = j; */
		}
	}

#ifdef CITAN_DEBUG
	if (a2p < 1e-15 || b2p < 1e-15)
		fprintf(stderr, "CITAN_DEBUG: index_lp_finite: a2p < 1e-15 || b2p < 1e-15\n");
#endif

	out[0] = pow(a2p,1.0/P);
	out[1] = pow(b2p,1.0/P);
}


/** Function to compute the l_p-index, O(n) time.
 *  @param x vector of non-negative numbers, sorted non-increasingly
 *  @param n pointer to the number of observations, n >= 1
 *  @param p pointer to the index order
 *  @param out two-dimensional array which stores the result, (a,b)
 */
void index_lp_infinite(double* x, int* n, double* out)
{
	int N = *n;
	int k = 0;
	double max = x[0];
	int i;

	if (x[0] <= 0.0) { out[0] = 0.0; out[1] = 0.0; return; }

	for (i=1; i<N; i++)
	{
		double maxcand = x[i]*(double)(i+1);
		if (max < maxcand)
		{
			k = i;
			max = maxcand;
		}
	}

	out[0] = (double)(k+1);
	out[1] = x[k];
}







/** Function to compute the S-statistic for kappa=id, O(log n) time.
 *  @param x vector of numbers, 0<=x[i]<=1, sorted non-increasingly
 *  @param n pointer to the number of observations
 *  @param out one-dimensional array which stores the result
 */
void Sstat2(double* x, int* n, double* out)
{
	double N = (double)(*n);

	int h1 = 0;
	int h2 = (*n)-1;
	int m;
	double xmulN;
	double mp1;

	if (x[0] < 1.0/N) { *out = x[0]; return; }

	while (1)
	{
		m = (h2+h1+1)/2;
		mp1 = (double)(m+1);
		xmulN = N*x[m];
		if (xmulN == mp1 || h1 == h2) {break;}
		if (xmulN < mp1) h2 = m-1;
		else h1 = m;
	}

#ifdef CITAN_DEBUG
	if (!(m+1 <= *n && m+1>=0)) fprintf(stderr, "CITAN_DEBUG: Sstat2: !(m+1 <= *n && m+1>=0)\n");
	if (m>=0 && x[m]<mp1/N) fprintf(stderr, "CITAN_DEBUG: Sstat2: m>=0 && x[m]<mp1/N\n");
#endif

	if (m+1 < *n)
	{
		if (mp1 > N*x[m+1]) *out = mp1/N;
		else                *out = x[m+1];
	}
	else
	{
		*out = mp1/N;
	}
}



/*
void Sstat2(double* x, int* n, double* out) -- OFTEN SLOWER THAN THE ABOVE
{
	int i = 0;
	int k = *n;
	double d = 1.0/(double)k;

	while (i<k)
	{
		if (x[i] < d)
		{
			if (x[i] >= (double)i/(double)k)
				*out = x[i];
			else
				*out = (double)i/(double)k;
			return;
		}
		++i;
		d = (double)(i+1)/(double)k;
	}

	// i == k
	*out = x[k-1];
}
*/








