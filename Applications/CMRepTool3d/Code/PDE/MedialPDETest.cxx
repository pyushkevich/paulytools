#include "MedialPDESolver.h"
#include <ctime>

// Define an arbitrary surface to test the differential geometry
// computation, and later to test the medial PDEs
class MyProblem1 : virtual public IMedialPDEProblem {
public:
  void ComputeJet2(double u, double v, double *X, double *Xu, double *Xv, 
    double *Xuu, double *Xuv, double *Xvv);
  double ComputeLaplacian(double u, double v)
    { return -1.0; }
};

void
MyProblem1
::ComputeJet2(double u, double v, double *X, double *Xu, double *Xv, 
  double *Xuu, double *Xuv, double *Xvv)
{
  X[0] = 3*u;
  X[1] = v;
  X[2] = 0.3*(-3 - 3*u + u*u + 4*v + 2*u*v - v*v);
  Xu[0] = 3;
  Xu[1] = 0;
  Xu[2] = 0.3*(-3 + 2*u + 2*v);
  Xv[0] = 0;
  Xv[1] = 1;
  Xv[2] = 0.3*(4 + 2*u - 2*v);
  Xuu[0] = 0;
  Xuu[1] = 0;
  Xuu[2] = 0.6;
  Xuv[0] = 0;
  Xuv[1] = 0;
  Xuv[2] = 0.6;
  Xvv[0] = 0;
  Xvv[1] = 0;
  Xvv[2] = -0.6;
}

int main(int argc, char *argv[])
{
  /* Step 1. Test the Differential Geometry code */
  double X[3], Xu[3], Xv[3], Xuu[3], Xuv[3], Xvv[3];
  double u = 0.3, v = 0.4;

  // Compute the Jet
  MyProblem1 p1;
  p1.ComputeJet2(u,v,X,Xu,Xv,Xuu,Xuv,Xvv);

  // Compute the differential geometry
  GeometryDescriptor gd(X,Xu,Xv,Xuu,Xuv,Xvv);
  gd.PrintSelf(cout);

  /* Step 2. Evaluate the equation at each site */
  MedialPDESolver mps(11,11);

  time_t t = clock();
  mps.Solve(&p1);
  cout << "*** Elapsed Time: " << (t-clock()) << " ms." << " ***" << endl;
}




