/*******************************************************************
GEODES/GEODES.C        Geodesic Computation          August 19, 1996

This routine approximates a minimum-length geodesic between two
given endpoints in multidimensional parameter space. It returns
the shortest curve found as a discrete set of parameter points
ordered by arc length.

********************************************************************

The C language function prototype is:
void geodes(int n, int m,
void yeval(double *,double *,double *,double *),
void ambmet(double *,double **),
int pmax, double accel, double eps, int itmax,
double *xend, double *xgeod);

Arguments are:
n ...... Dimension of parameter space x^i.
m ...... Dimension of ambient space y^r.
yeval .. Address of manifold function y^r = y^r(x^i).
ambmet.. Address of ambient space metric function
pmax ... Number of points on geodesic, not counting endpoints.
accel .. Acceleration factor.
eps .... Arc length ratio convergence criteria.
itmax .. Maximum number of iterations without convergence.
xend ... Given geodesic endpoints.
xgeod .. Resulting shortest length geodesic found.

Performance controls are provided by arguments pmax, accel, eps,
and itmax.

********************************************************************

This routine solves a system of n nonlinear ordinary differential
equations

d^2x^i/ds^2 + Gamma^i_{jk} dx^j/ds d^k/ds = 0,

where s is arc length, i = 1,2,...,n, and with implied summation on
j and k. Gamma is the Christoffel symbol of the second kind.
Boundary conditions are two geodesic endpoints. This routine uses a
relaxation (iterative) solving method.

There are two issues the user should understand. There can
be more than one minimum-length geodesic. This is common if
endpoints are located symmetrically on a symmetric geometry. This
routine approximates only one minimum-length geodesic.

This routine solves for any geodesic but returns the shortest curve
found. The distinction is important if the routine converges on a
geodesic other than a minimum-length geodesic. Success is indicated
by monitoring iteration statistics provided by compiler option
-D_DEBUG. Best results occur when convergence results in shorter
curves. If convergence is slow, the shortest curve found is the best
approximation of a minimum-length geodesic since a theorem proves an
actual shortest curve(s) is always a geodesic.

Reference:
Manfredo do Carmo, Riemannian Geometry, Birkhauser, Boston, 1992

********************************************************************

This software can be modified or extended by:
* Allowing indefinite metric with null geodesic(s).
* Optimizing grid points dynamically. Using recursion.
* Optimizing performance control values dynamically.
* Finding families of geodesics.
* Add localized salients (bumps) to the manifold.
* Plotting geodesic in both parameter and ambient space.

Refer comments and questions to:
William L. Anderson
Elements Research
7501 Windyrush Road
Charlotte, NC 28226
704-543-9180
elements@ix.netcom.com
http://www.netcom.com/~elements/

********************************************************************

Copyright (c) 1994-1996 Elements Research, William L. Anderson.

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose without fee is hereby granted,
provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear
in supporting documentation.

This software is provided "as is" without express or implied
warranty.

*******************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <malloc.h>
#include "geodesic.h"



AmbientMetric::AmbientMetric(int m_arg) {
	int r,s;
	m = m_arg;
	ambmet_matrix = (double *)malloc(sizeof(double)*(m*m));
	/* initialize as euclidean metric */
	for (r = 0; r < m; r++)
		for (s = 0; s < m; s++)
			if (r == s)
				ambmet_matrix[SUBS2(m,r,s)] = 1.0;
			else
				ambmet_matrix[SUBS2(m,r,s)] = 0.0;
}

AmbientMetric::~AmbientMetric() {
	free(ambmet_matrix);
	ambmet_matrix = NULL;
	m = 0;
}

void EuclideanMetric::metric(double *y, double **ambmet_matrix_arg) {
	*ambmet_matrix_arg = ambmet_matrix;
}

void HyperbolicMetric::metric(double *y, double **ambmet_matrix_arg) {
	int r;
	double hyperbolic;
	hyperbolic = 1.0/(y[0]*y[0]);
	for (r = 0; r < m; r++)
		ambmet_matrix[SUBS2(m,r,r)] = hyperbolic;
	*ambmet_matrix_arg = ambmet_matrix;
}

void SurfaceGeodesic::mult(double *a, double *b, int m, int n,int h, double *c) {
	/* multiply two matrices */
	int i,j,k;
	double sum;
	for (i = 0; i < m; i++)
		for (j = 0; j < h; j++) {
			sum = 0.0;
			for (k = 0; k < n; k++)
				sum += a[SUBS2(n,i,k)]*b[SUBS2(h,k,j)];
			c[SUBS2(h,i,j)] = sum;
		}
}

void SurfaceGeodesic::multt(double *a, double *b, int m, int n,int h, double *c) {
	/* multiply matrix transpose with another matrix */
	int i,j,k;
	double sum;
	for (i = 0; i < n; i++)
		for (j = 0; j < h; j++) {
			sum = 0.0;
 			for (k = 0; k < m; k++)
				sum += a[SUBS2(n,k,i)]*b[SUBS2(h,k,j)];
			c[SUBS2(h,i,j)] = sum;
		}
}

void SurfaceGeodesic::metric_comp(int n, int m,
						double *jacobian, double *y, AmbientMetric *ambientMetric,
						double *jacobianTG, double *metric) 
{
	double *ambient_metric_matrix;
	int i,j,r;
	double sum;
	ambientMetric->metric(y,&ambient_metric_matrix);
	multt(jacobian,ambient_metric_matrix,m,n,m,jacobianTG);
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++) {
			sum = 0.0;
			for (r = 0; r < m; r++)
				sum += jacobianTG[SUBS2(m,i,r)]*jacobian[SUBS2(n,r,j)];
			metric[SUBS2(n,i,j)] = sum;
			if (j != i)
				metric[SUBS2(n,j,i)] = sum;
		}
}

void SurfaceGeodesic::metric_inv_comp(int n, double *metric,
							double *work, double *diag, double *metric_inv) 
{
	int i,j,k;
	double sum;
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)
			work[SUBS2(n,i,j)] = metric[SUBS2(n,i,j)];
	/* Cholesky decomposition */
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++) {
			sum = work[SUBS2(n,i,j)];
			for (k = 0; k < i; k++)
				sum -= work[SUBS2(n,i,k)]*work[SUBS2(n,j,k)];
			if (i == j) {
				if (sum <= 0.0)
					printf("Metric tensor is not positive definite\n");
				diag[i] = sqrt(sum);
			}
			else
				work[SUBS2(n,j,i)] = sum/diag[i];
		}
		/* Cholesky L inverse */
		for (i = 0; i < n; i++) {
			work[SUBS2(n,i,i)] = 1.0/diag[i];
			for (j = i+1; j < n; j++) {
				sum = 0.0;
				for (k = i; k < j; k++)
					sum -= work[SUBS2(n,j,k)]*work[SUBS2(n,k,i)];
				work[SUBS2(n,j,i)] = sum/diag[j];
				work[SUBS2(n,i,j)] = 0.0;
			}
		}
		/* metric inverse */
		for (i = 0; i < n; i++)
			for (j = i; j < n; j++) {
				sum = 0.0;
				for (k = 0; k < n; k++)
					sum += work[SUBS2(n,k,i)]*work[SUBS2(n,k,j)];
				metric_inv[SUBS2(n,i,j)] = sum;
				if (j != i)
					metric_inv[SUBS2(n,j,i)] = sum;
			}
}

void SurfaceGeodesic::christoffel_comp(int n, int m,
							 double *jacobianTG, double *hessian, double *christoffel_first) 
{
	int i,j,k,r;
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)
			for (k = 0; k < n; k++) {
				christoffel_first[SUBS3(n,n,k,i,j)] = 0.0;
				for (r = 0; r < m; r++)
					christoffel_first[SUBS3(n,n,k,i,j)] +=
					jacobianTG[SUBS2(m,k,r)]*hessian[SUBS3(n,n,r,i,j)];
				if (j != i)
					christoffel_first[SUBS3(n,n,k,j,i)] =
					christoffel_first[SUBS3(n,n,k,i,j)];
			}
}

void SurfaceGeodesic::tridiagonal_init(int pmax) 
{
	int p;
	pivot = (double *) malloc(sizeof(double)*pmax);
	for (p = 1; p < pmax-1; p++)
		pivot[p] = (double)(p)/(double)(p+1);
}

void SurfaceGeodesic::tridiagonal_free() 
{
	free(pivot);
}

void SurfaceGeodesic::tridiagonal(int pmax, double *b) 
{
	/* tridiagonal (1; 1 -2 1; ...; 1) solver */
	int p;
	/* forward elimination */
	for (p = 1; p < pmax-1; p++) {
		b[p] -= b[p-1];
		b[p] *= -pivot[p];
	}
	/* back substitution */
	for (p = pmax-2; p > 0; p--)
		b[p] += pivot[p] * b[p+1];
}

void SurfaceGeodesic::geodes(int n, int m,
			GeodesicManifold *manifold,
			AmbientMetric *ambientMetric,
			int pmax, double accel, double eps, int itmax,
			double *xend, double *xgeod) 
{
	int it;
	int p;
	int i,j,k;
	double *x,*y,*jacobian,*hessian;
	double *jacobianTG;
	double *metric,*metric_work,*metric_diag,*metric_inv;
	double *christoffel_first,*christoffel_second;
	double *dx,*dxds,*csdxdsdxds;
	double *ds,*s;
	double *rhs;
	double prev_arc_len,prev_arc_len_ratio;
	double s_shortest;
	double *xshortest;
	x = (double *) malloc(sizeof(double)*pmax*n);
	y = (double *) malloc(sizeof(double)*pmax*m);
	jacobian = (double *) malloc(sizeof(double)*(m*n));
	hessian = (double *) malloc(sizeof(double)*(m*n*n));
	jacobianTG = (double *) malloc(sizeof(double)*(n*m));
	metric = (double *) malloc(sizeof(double)*(n*n));
	metric_work = (double *) malloc(sizeof(double)*(n*n));
	metric_diag = (double *) malloc(sizeof(double)*n);
	metric_inv = (double *) malloc(sizeof(double)*(n*n));
	christoffel_first = (double *) malloc(sizeof(double)*(n*n*n));
	christoffel_second = (double *) malloc(sizeof(double)*(n*n*n));
	dx = (double *) malloc(sizeof(double)*n);
	ds = (double *) malloc(sizeof(double)*pmax);
	dxds = (double *) malloc(sizeof(double)*n);
	csdxdsdxds = (double *) malloc(sizeof(double)*n);
	s = (double *) malloc(sizeof(double)*(pmax+2));
	rhs = (double *) malloc(sizeof(double)*n*(pmax+2));
	xshortest = (double *) malloc(sizeof(double)*pmax*n);
	tridiagonal_init(pmax+2);
	/* initial estimate is straight line in parameter space */
	for (p = 0; p < pmax; p++)
		for (i = 0; i < n; i++)
			x[SUBS2(n,p,i)] = xend[SUBS2(n,0,i)] +
			(xend[SUBS2(n,1,i)] - xend[SUBS2(n,0,i)])*
			(double)(p+1)/(double)(pmax+1);
	prev_arc_len = DBL_MAX;
	s_shortest = DBL_MAX;
	for (it = 0; it < itmax; it++) {
		for (p = 0; p < pmax; p++) {
			manifold->eval(&x[SUBS2(n,p,0)],&y[SUBS2(m,p,0)],jacobian,hessian);
			metric_comp(n,m,jacobian,&y[SUBS2(m,p,0)],ambientMetric,jacobianTG,metric);
			metric_inv_comp(n,metric,metric_work,metric_diag,metric_inv);
			christoffel_comp(n,m,jacobianTG,hessian,christoffel_first);
			mult(metric_inv,christoffel_first,n,n,n*n,christoffel_second);
			for (i = 0; i < n; i++) {
				if (p == 0)
					dx[i] = x[SUBS2(n,1,i)] - xend[SUBS2(n,0,i)];
				else if (p == pmax-1)
					dx[i] = xend[SUBS2(n,1,i)] - x[SUBS2(n,p-1,i)];
				else
					dx[i] = x[SUBS2(n,p+1,i)] - x[SUBS2(n,p-1,i)];
				dx[i] /= 2.0;
			}
			ds[p] = 0.0;
			for (i = 0; i < n; i++) {
				ds[p] += metric[SUBS2(n,i,i)]*dx[i]*dx[i];
				for (j = i+1; j < n; j++)
					ds[p] += 2.0*metric[SUBS2(n,i,j)]*dx[i]*dx[j];
			}
			ds[p] = sqrt(ds[p]);
			for (i = 0; i < n; i++)
				dxds[i] = dx[i]/ds[p];
			for (k = 0; k < n; k++) {
				csdxdsdxds[k] = 0.0;
				for (i = 0; i < n; i++) {
					csdxdsdxds[k] +=
						christoffel_second[SUBS3(n,n,k,i,i)]*dxds[i]*dxds[i];
					for (j = i+1; j < n; j++)
						csdxdsdxds[k] += 2.0*
						christoffel_second[SUBS3(n,n,k,i,j)]*dxds[i]*dxds[j];
				}
			}
			for (i = 0; i < n; i++)
				rhs[SUBS2(pmax+2,i,p+1)] = -csdxdsdxds[i]*ds[p]*ds[p];
		}
		/* arc length */
		s[0] = 0.0;
		s[1] = ds[0];
		for (p = 1; p < pmax; p++)
			s[p+1] = s[p] + (ds[p-1] + ds[p])/2.0;
		s[pmax+1] = s[pmax] + ds[pmax-1];
#ifdef _DEBUG
		/* iteration statistics */
		printf("it=%d arc length=%12.5g\n",it,s[pmax+1]);
#endif
		/* save curve with shortest arc length */
		if (s[pmax+1] < s_shortest) {
			s_shortest = s[pmax+1];
			for (p = 0; p < pmax; p++)
				for (i = 0; i < n; i++)
					xshortest[SUBS2(n,p,i)] = x[SUBS2(n,p,i)];
		}
		/* stopping criteria */
		prev_arc_len_ratio = (s[pmax+1] - prev_arc_len)/prev_arc_len;
		if (fabs(prev_arc_len_ratio) < eps)
			break;
		prev_arc_len = s[pmax+1];
		/* solve */
		for (i = 0; i < n; i++) {
			rhs[SUBS2(pmax+2,i,0)] = xend[SUBS2(n,0,i)];
			rhs[SUBS2(pmax+2,i,pmax+1)] = xend[SUBS2(n,1,i)];
			tridiagonal(pmax+2,&rhs[SUBS2(pmax+2,i,0)]);
		}
		/* update x for next iteration */            
		for (p = 0; p < pmax; p++)
			for (i = 0; i < n; i++)
				x[SUBS2(n,p,i)] +=
				accel*(rhs[SUBS2(pmax+2,i,p+1)] - x[SUBS2(n,p,i)]);
	}
	/* copy to xgeod */
	for (i = 0; i < n; i++) {
		xgeod[SUBS2(n,0,i)] = xend[SUBS2(n,0,i)];
		xgeod[SUBS2(n,pmax+1,i)] = xend[SUBS2(n,1,i)];
	}
	for (p = 0; p < pmax; p++)
		for (i = 0; i < n; i++)
			xgeod[SUBS2(n,p+1,i)] = xshortest[SUBS2(n,p,i)];
	tridiagonal_free();
	free(x);
	free(y);
	free(jacobian);
	free(hessian);
	free(jacobianTG);
	free(metric);
	free(metric_work);
	free(metric_diag);
	free(metric_inv);
	free(christoffel_first);
	free(christoffel_second);
	free(dx);
	free(ds);
	free(dxds);
	free(csdxdsdxds);
	free(s);
	free(rhs);
	free(xshortest);
}

