#ifndef _GEODESIC_H_
#define _GEODESIC_H_

/* macros reduce multidimensional array to C language array */
#define SUBS2(n,i,j) ((i)*(n)+(j))
#define SUBS3(n,m,i,j,k) (((i)*(n)+(j))*(m)+(k))

class AmbientMetric {
protected:
	int m;
	double *ambmet_matrix;
public:
	AmbientMetric(int m_arg);
	virtual ~AmbientMetric();
	
	virtual void metric(double *y, double **ambmet_matrix_arg) = 0;
};

class EuclideanMetric : public AmbientMetric {
public:
	EuclideanMetric(int m_arg) : AmbientMetric(m_arg) {};
	
	void metric(double *y, double **ambmet_matrix_arg);
};

class HyperbolicMetric : public AmbientMetric {
public:
	HyperbolicMetric(int m_arg) : AmbientMetric(m_arg) {};
	
	void metric(double *y, double **ambmet_matrix_arg);
};

class GeodesicManifold {
private:
	int n,m;
public:
	GeodesicManifold(int n,int m) {
		this->m = m;
		this->n = n;
	}

	int getDomainDim() {
		return n;
	}

	int getRangeDim() {
		return m;
	}

	virtual void eval(double *x,double *y,double *jacobian,double *hessian) = 0;
};

class SurfaceGeodesic {
private:
	double *pivot;

	void mult(double *a, double *b, int m, int n,int h, double *c);
	void multt(double *a, double *b, int m, int n,int h, double *c);

	void metric_comp(int n, int m,double *jacobian, double *y, AmbientMetric *ambientMetric,
				double *jacobianTG, double *metric);
	void metric_inv_comp(int n, double *metric,double *work, double *diag, double *metric_inv);


	void christoffel_comp(int n, int m, double *jacobianTG, double *hessian, double *christoffel_first);

	void tridiagonal_init(int pmax);
	void tridiagonal_free(); 	
	void tridiagonal(int pmax, double *b);

public:
	void geodes(int n, int m,GeodesicManifold *manifold,
				AmbientMetric *ambientMetric,
				int pmax, double accel, double eps, int itmax,
				double *xend, double *xgeod);
};

#endif // _GEODESIC_H_