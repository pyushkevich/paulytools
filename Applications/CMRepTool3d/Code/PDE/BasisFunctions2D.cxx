#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include "BasisFunctions2D.h"
#include "BasisFunctions2D.txx"


double CosineBasisFunction::Evaluate(double u, size_t k, size_t d)
{
  // Shift u so we don't have zero derivatives
  double a = 0.8, b = 0.1;

  // Compute the basis at the point
  double T1 = M_PI * k * a, T2 = M_PI * k * (a * u + b);

  switch(d) 
    {
  case 0: return cos(T2);
  case 1: return -T1 * sin(T2);
  case 2: return T1 * T1 * cos(T2);
  default: return 0.0;
    }
}


void 
FourierSurface
::ApplyAffineTransform(const vnl_matrix<double> &A, 
  const vnl_vector<double> &b, const vnl_vector<double> &c)
{
  // Sample points from the medial surface and apply rotation
  for(size_t i = 0; i < ncu; i++) for(size_t j = 0; j < ncv; j++)
    {
    SMLVec3d y, x( C(0, i, j), C(1, i, j), C(2, i, j));

    if(i == 0 && j == 0)
      {
      y = A * (x - c) + b;
      }
    else
      {
      y = A * x;
      }

    C(0, i, j) = y[0]; C(1, i, j) = y[1]; C(2, i, j) = y[2];
    }

}




template class GenericBasisRepresentation2D< 3, 3, 
         CosineBasisFunction, CosineBasisFunction>;
