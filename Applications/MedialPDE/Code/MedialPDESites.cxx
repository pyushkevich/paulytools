#include "MedialPDESites.h"

double FDInternalSite::ComputeEquation(const Mat &Y)
{
  double Fu, Fv, Fuu, Fuv, Fvv;
  
  // Compute the partial derivatives of F with respect to u and v
  xMask->ComputeTwoJet(Y, Fu, Fv, Fuu, Fuv, Fvv);

  // Return the sum weighted by geometry
  return Cuu * Fuu + Cuv * Fuv + Cvv * Fvv + Cu * Fu + Cv * Fv - rho;
};

/** Compute the derivatives for Newton's method */
void FDInternalSite
::ComputeDerivative(const Mat &, double *A, unsigned int iRow)
{
  // Because the equation is linear, there is no need to use values of 
  // Phi. We just copy the weightings
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivative[i];
}

void FDInternalSite::SetGeometry(GeometryDescriptor *g, double rho)
{
  // Store the rho
  this->rho = rho;

  // Compute the factors applied to partials of Phi at the site
  Cuu = g->xContravariantTensor[0][0];
  Cvv = g->xContravariantTensor[1][1];
  Cuv = g->xContravariantTensor[1][0] * 2.0;
  Cu = - (
    g->xContravariantTensor[0][0] * g->xChristoffelSecond[0][0][0] + 
    g->xContravariantTensor[0][1] * g->xChristoffelSecond[0][1][0] * 2.0 +
    g->xContravariantTensor[1][1] * g->xChristoffelSecond[1][1][0]);
  Cv = - (
    g->xContravariantTensor[0][0] * g->xChristoffelSecond[0][0][1] + 
    g->xContravariantTensor[0][1] * g->xChristoffelSecond[0][1][1] * 2.0 +
    g->xContravariantTensor[1][1] * g->xChristoffelSecond[1][1][1]);

  // Compute the derivative with respect to each cite
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    // Get the contribution of node i to each partial derivative
    double Wu, Wv, Wuu, Wuv, Wvv;
    xMask->GetNodeComponentInTwoJet(i, Wu, Wv, Wuu, Wuv, Wvv);
    
    // Add them with weighting to produce the derivatives
    xDerivative[i] = Cu * Wu + Cv * Wv + Cuu * Wuu + Cuv * Wuv + Cvv * Wvv;
    }
}

/** 
 * Compute the variational derivative 
 */
void
FDInternalSite::
ComputeVariationalDerivative(const Mat &Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // The A[i] elements are computed in the same way as for Newton's method
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivative[i];

  // The rest of this routine computes the right hand side

  // First get shorthands for the relevant vectors
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Xu = xAtom->Xu;
  SMLVec3d &Nv = dAtom->Xv; SMLVec3d &Xv = xAtom->Xv;
  SMLVec3d &Nuu = dAtom->Xuu; SMLVec3d &Xuu = xAtom->Xuu;
  SMLVec3d &Nuv = dAtom->Xuv; SMLVec3d &Xuv = xAtom->Xuv;
  SMLVec3d &Nvv = dAtom->Xvv; SMLVec3d &Xvv = xAtom->Xvv;

  // Get shorthand for the differential geometric operators
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;
  double (*K)[2][2] = xAtom->G.xChristoffelFirst;
  double (*K2)[2][2] = xAtom->G.xChristoffelSecond;
  double &g = xAtom->G.g;
  
  // Compute the derivatives of the contravariant tensor
  double NuXu = dot_product(Nu, Xu);
  double NvXu = dot_product(Nv, Xu);
  double NuXv = dot_product(Nu, Xv);
  double NvXv = dot_product(Nv, Xv);

  // The derivative of 'g'
  double dgdN = dAtom->G.g = 
    2 * ( NuXu * G1[1][1] + NvXv * G1[0][0] - (NuXv + NvXu) * G1[0][1] );
  double g2 = g * g;

  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 
  z[0][0] = (2 * NvXv * g - G1[1][1] * dgdN) / g2;
  z[1][1] = (2 * NuXu * g - G1[0][0] * dgdN) / g2;
  z[0][1] = z[1][0] = - ((NuXv + NvXu) * g - G1[1][0] * dgdN) / g2;

  // Compute the derivative of the Christoffel symbols of the second kind
  double NuXuu = dot_product(Nu, Xuu), NvXuu = dot_product(Nv, Xuu);
  double NuXuv = dot_product(Nu, Xuv), NvXuv = dot_product(Nv, Xuv);
  double NuXvv = dot_product(Nu, Xvv), NvXvv = dot_product(Nv, Xvv);
  double NuuXu = dot_product(Nuu, Xu), NuuXv = dot_product(Nuu, Xv);
  double NuvXu = dot_product(Nuv, Xu), NuvXv = dot_product(Nuv, Xv);
  double NvvXu = dot_product(Nvv, Xu), NvvXv = dot_product(Nvv, Xv);
  
  // The derivative of Christofel symbols of second kind  
  double (*Q)[2][2] = dAtom->G.xChristoffelSecond;
  Q[0][0][0] = 
    z[0][0] * K[0][0][0] + G2[0][0] * ( NuXuu + NuuXu ) +
    z[0][1] * K[0][0][1] + G2[0][1] * ( NvXuu + NuuXv );
  Q[0][0][1] = 
    z[1][0] * K[0][0][0] + G2[1][0] * ( NuXuu + NuuXu ) +
    z[1][1] * K[0][0][1] + G2[1][1] * ( NvXuu + NuuXv );
  
  Q[0][1][0] = Q[1][0][0] = 
    z[0][0] * K[1][0][0] + G2[0][0] * ( NuXuv + NuvXu ) +
    z[0][1] * K[1][0][1] + G2[0][1] * ( NvXuv + NuvXv );
  Q[0][1][1] = Q[1][0][1] = 
    z[1][0] * K[1][0][0] + G2[1][0] * ( NuXuv + NuvXu ) +
    z[1][1] * K[1][0][1] + G2[1][1] * ( NvXuv + NuvXv );

  Q[1][1][0] = 
    z[0][0] * K[1][1][0] + G2[0][0] * ( NuXvv + NvvXu ) +
    z[0][1] * K[1][1][1] + G2[0][1] * ( NvXvv + NvvXv );
  Q[1][1][1] = 
    z[1][0] * K[1][1][0] + G2[1][0] * ( NuXvv + NvvXu ) +
    z[1][1] * K[1][1][1] + G2[1][1] * ( NvXvv + NvvXv );

  // Compute the partials of phi
  double Fu, Fv, Fuu, Fuv, Fvv;
  xMask->ComputeTwoJet(Y, Fu, Fv, Fuu, Fuv, Fvv);

  // Compute the derivative (this should be the derivative)
  *b = dAtom->xLapR - ( 
    z[0][0] * (Fuu - K2[0][0][0] * Fu - K2[0][0][1] * Fv) - 
    G2[0][0] * (Q[0][0][0] * Fu + Q[0][0][1] * Fv) +
    z[1][1] * (Fvv - K2[1][1][0] * Fu - K2[1][1][1] * Fv) - 
    G2[1][1] * (Q[1][1][0] * Fu + Q[1][1][1] * Fv) +
    2.0 * ( 
      z[0][1] * (Fuv - K2[0][1][0] * Fu - K2[0][1][1] * Fv) - 
      G2[0][1] * (Q[0][1][0] * Fu + Q[0][1][1] * Fv) ));
}


/** 
 * Compute the variational derivative with respect to N
 */
/*
void
FDInternalSite::
ComputeVariationalDerivativeX(const Mat &Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // The A[i] elements are computed in the same way as for Newton's method
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivative[i];

  // The rest of this routine computes the right hand side

  // First get shorthands for the relevant vectors
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Xu = xAtom->Xu;
  SMLVec3d &Nv = dAtom->Xv; SMLVec3d &Xv = xAtom->Xv;
  SMLVec3d &Nuu = dAtom->Xuu; SMLVec3d &Xuu = xAtom->Xuu;
  SMLVec3d &Nuv = dAtom->Xuv; SMLVec3d &Xuv = xAtom->Xuv;
  SMLVec3d &Nvv = dAtom->Xvv; SMLVec3d &Xvv = xAtom->Xvv;

  // Get shorthand for the differential geometric operators
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;
  double (*K)[2][2] = xAtom->G.xChristoffelFirst;
  double (*K2)[2][2] = xAtom->G.xChristoffelSecond;
  double &g = xAtom->G.g;
  
  // Compute the derivatives of the contravariant tensor
  double NuXu = dot_product(Nu, Xu);
  double NvXu = dot_product(Nv, Xu);
  double NuXv = dot_product(Nu, Xv);
  double NvXv = dot_product(Nv, Xv);

  // The derivative of 'g'
  double dgdN = dAtom->G.g = 
    2 * ( NuXu * G1[1][1] + NvXv * G1[0][0] - (NuXv + NvXu) * G1[0][1] );
  double g2 = g * g;

  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 
  z[0][0] = (2 * NvXv * g - G1[1][1] * dgdN) / g2;
  z[1][1] = (2 * NuXu * g - G1[0][0] * dgdN) / g2;
  z[0][1] = z[1][0] = - ((NuXv + NvXu) * g - G1[1][0] * dgdN) / g2;

  // Compute the derivative of the Christoffel symbols of the second kind
  double NuXuu = dot_product(Nu, Xuu), NvXuu = dot_product(Nv, Xuu);
  double NuXuv = dot_product(Nu, Xuv), NvXuv = dot_product(Nv, Xuv);
  double NuXvv = dot_product(Nu, Xvv), NvXvv = dot_product(Nv, Xvv);
  double NuuXu = dot_product(Nuu, Xu), NuuXv = dot_product(Nuu, Xv);
  double NuvXu = dot_product(Nuv, Xu), NuvXv = dot_product(Nuv, Xv);
  double NvvXu = dot_product(Nvv, Xu), NvvXv = dot_product(Nvv, Xv);
  
  // The derivative of Christofel symbols of second kind  
  double (*Q)[2][2] = dAtom->G.xChristoffelSecond;
  Q[0][0][0] = 
    z[0][0] * K[0][0][0] + G2[0][0] * ( NuXuu + NuuXu ) +
    z[0][1] * K[0][0][1] + G2[0][1] * ( NvXuu + NuuXv );
  Q[0][0][1] = 
    z[1][0] * K[0][0][0] + G2[1][0] * ( NuXuu + NuuXu ) +
    z[1][1] * K[0][0][1] + G2[1][1] * ( NvXuu + NuuXv );
  
  Q[0][1][0] = Q[1][0][0] = 
    z[0][0] * K[1][0][0] + G2[0][0] * ( NuXuv + NuvXu ) +
    z[0][1] * K[1][0][1] + G2[0][1] * ( NvXuv + NuvXv );
  Q[0][1][1] = Q[1][0][1] = 
    z[1][0] * K[1][0][0] + G2[1][0] * ( NuXuv + NuvXu ) +
    z[1][1] * K[1][0][1] + G2[1][1] * ( NvXuv + NuvXv );

  Q[1][1][0] = 
    z[0][0] * K[1][1][0] + G2[0][0] * ( NuXvv + NvvXu ) +
    z[0][1] * K[1][1][1] + G2[0][1] * ( NvXvv + NvvXv );
  Q[1][1][1] = 
    z[1][0] * K[1][1][0] + G2[1][0] * ( NuXvv + NvvXu ) +
    z[1][1] * K[1][1][1] + G2[1][1] * ( NvXvv + NvvXv );

  // Compute the partials of phi
  double Fu, Fv, Fuu, Fuv, Fvv;
  xMask->ComputeTwoJet(Y, Fu, Fv, Fuu, Fuv, Fvv);

  // Compute the derivative (this should be the derivative)
  *b = -( 
    z[0][0] * (Fuu - K2[0][0][0] * Fu - K2[0][0][1] * Fv) - 
    G2[0][0] * (Q[0][0][0] * Fu + Q[0][0][1] * Fv) +
    z[1][1] * (Fvv - K2[1][1][0] * Fu - K2[1][1][1] * Fv) - 
    G2[1][1] * (Q[1][1][0] * Fu + Q[1][1][1] * Fv) +
    2.0 * ( 
      z[0][1] * (Fuv - K2[0][1][0] * Fu - K2[0][1][1] * Fv) - 
      G2[0][1] * (Q[0][1][0] * Fu + Q[0][1][1] * Fv) ));
}
*/

/*
void
FDInternalSite::
ComputeVariationalDerivativeRho(const Mat &Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // The A[i] elements are computed in the same way as for Newton's method
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivative[i];

  // The right hand side is just the derivative of rho
  *b = dAtom->xLapR;
}
*/

void 
FDInternalSite::PrintReport()
{
	cout << "Internal Site at " << xMask->GetLocationU() 
		<< "," << xMask->GetLocationV() << endl;
	cout << "LB(F) = " 
		<< Cuu << " * Fuu + " 
		<< Cuv << " * Fuv + "
		<< Cvv << " * Fvv + "
		<< Cu  << " * Fu + "
		<< Cv  << " * Fv " << endl;
}

double FDBorderSite::ComputeEquation(const Mat &Y)
{
  double F, Fu, Fv;
  
  // Compute the partial derivatives of F with respect to u and v
  F = xMask->ComputeOneJet(Y, Fu, Fv);

  // Compute the generalized gradient magnitude using precomputed values
  return CuCu * Fu * Fu + CuCv * Fu * Fv + CvCv * Fv * Fv - 4.0 * F;
}

/** Compute the derivatives for Newton's method */
void FDBorderSite
::ComputeDerivative(const Mat &Y, double *A, unsigned int iRow)
{
  double F, Fu, Fv;
  
  // Compute the partial derivatives of F with respect to u and v
  F = xMask->ComputeOneJet(Y, Fu, Fv);

  // Because the equation is linear, there is no need to use values of 
  // Phi. We just copy the weightings
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivativeFu[i] * Fu + xDerivativeFv[i] * Fv + xDerivativeF[i];
}

void
FDBorderSite
::ComputeVariationalDerivative(const Mat &Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  double F, Fu, Fv;

  // First, let's compute grad phi
  F = xMask->ComputeOneJet(Y, Fu, Fv);

  // Next, compute the weights for Hu and Hv and H
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;

  // Set the entries in the derivative matrix. Should be the same as
  // in the usual case...
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    // Get the contribution of node i to the jet of H
    double Wu, Wv, W;
    W = xMask->GetNodeComponentInOneJet(i, Wu, Wv);
    
    // Compute the equation
    A[i] = -4.0 * W + 
      2.0 * (G2[0][0] * Wu + G2[0][1] * Wv) * Fu + 
      2.0 * (G2[1][0] * Wu + G2[1][1] * Wv) * Fv;
    }
    // A[i] = xDerivativeFu[i] * Fu + xDerivativeFv[i] * Fv + xDerivativeF[i];

  // Compute the derivatives of the contravariant tensor
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Nv = dAtom->Xv;
  SMLVec3d &Xu = xAtom->Xu; SMLVec3d &Xv = xAtom->Xv;

  double NuXu = dot_product(Nu, Xu);
  double NvXu = dot_product(Nv, Xu);
  double NuXv = dot_product(Nu, Xv);
  double NvXv = dot_product(Nv, Xv);
  
  // The derivative of 'g'
  double &g = xAtom->G.g;
  double dgdN = dAtom->G.g = 
    2 * ( NuXu * G1[1][1] + NvXv * G1[0][0] - (NuXv + NvXu) * G1[0][1] );
  double g2 = g * g;

  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 
  z[0][0] = (2 * NvXv * g - G1[1][1] * dgdN) / g2;
  z[1][1] = (2 * NuXu * g - G1[0][0] * dgdN) / g2;
  z[0][1] = z[1][0] = - ((NuXv + NvXu) * g - G1[1][0] * dgdN) / g2;

  // Compute the right hand side
  *b = - (z[0][0] * Fu * Fu + 2.0 * z[0][1] * Fu * Fv + z[1][1] * Fv * Fv);
}

/*
void
FDBorderSite
::ComputeVariationalDerivativeX(const Mat &Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  double F, Fu, Fv;

  // First, let's compute grad phi
  F = xMask->ComputeOneJet(Y, Fu, Fv);

  // Next, compute the weights for Hu and Hv and H
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;

  // Set the entries in the derivative matrix. Should be the same as
  // in the usual case...
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivativeFu[i] * Fu + xDerivativeFv[i] * Fv + xDerivativeF[i];

  // Compute the derivatives of the contravariant tensor
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Nv = dAtom->Xv;
  SMLVec3d &Xu = xAtom->Xu; SMLVec3d &Xv = xAtom->Xv;

  double NuXu = dot_product(Nu, Xu);
  double NvXu = dot_product(Nv, Xu);
  double NuXv = dot_product(Nu, Xv);
  double NvXv = dot_product(Nv, Xv);
  
  // The derivative of 'g'
  double &g = xAtom->G.g;
  double dgdN = dAtom->G.g = 
    2 * ( NuXu * G1[1][1] + NvXv * G1[0][0] - (NuXv + NvXu) * G1[0][1] );
  double g2 = g * g;

  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 
  z[0][0] = (2 * NvXv * g - G1[1][1] * dgdN) / g2;
  z[1][1] = (2 * NuXu * g - G1[0][0] * dgdN) / g2;
  z[0][1] = z[1][0] = - ((NuXv + NvXu) * g - G1[1][0] * dgdN) / g2;

  // Compute the right hand side
  *b = - (z[0][0] * Fu * Fu + 2.0 * z[0][1] * Fu * Fv + z[1][1] * Fv * Fv);
}

void
FDBorderSite
::ComputeVariationalDerivativeRho(const Mat &Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  double F, Fu, Fv;

  // First, let's compute grad phi
  F = xMask->ComputeOneJet(Y, Fu, Fv);

  // We need the contravariant tensor
  double (*G2)[2] = xAtom->G.xContravariantTensor;

  // Set the entries in the derivative matrix. Should be the same as
  // in the usual case...
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    // Get the contribution of node i to the jet of H
    double Wu, Wv, W;
    W = xMask->GetNodeComponentInOneJet(i, Wu, Wv);
    
    // Compute the equation
    A[i] = -2.0 * W + 
      (G2[0][0] * Wu + G2[0][1] * Wv) * Fu + 
      (G2[1][0] * Wu + G2[1][1] * Wv) * Fv;
    }

  // The right hand side is just 0
  *b = 0.0;
}
*/

/** Initialize the border site */
void FDBorderSite::SetGeometry(GeometryDescriptor *g, double)
{
  // Compute the weights of the terms in the equation
  CuCu = g->xContravariantTensor[0][0];
  CuCv = g->xContravariantTensor[0][1] * 2.0;
  CvCv = g->xContravariantTensor[1][1];

  // Compute the derivatives of Fu and Fv with respect to each node
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    double W, Wu, Wv;
    
    // Get the contribution of node i to the partials 
    W = xMask->GetNodeComponentInOneJet(i, Wu, Wv);
    
    // The derivative of Fu with respect to neighbor i
    xDerivativeFu[i] = 2.0 * CuCu * Wu + CuCv * Wv;
      
    // The derivative of Fv with respect to neighbor i
    xDerivativeFv[i] = 2.0 * CvCv * Wv + CuCv * Wu;
      
    // The derivative of F with respect to neighbor i
    xDerivativeF[i] = -4.0 * W;
    }
}

void 
FDBorderSite::PrintReport()
{
	cout << "Border Site at " << xMask->GetLocationU() 
		<< "," << xMask->GetLocationV() << endl;
	cout << "LB(F) = " 
		<< CuCu << " * Fu * Fu + " 
		<< CuCv << " * Fu * Fv + "
		<< CvCv << " * Fv * Fv + " << endl;
}

