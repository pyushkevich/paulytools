#include "MedialAtom.h"

void
MedialAtom
::SetAllDerivativeTermsToZero()
{
  X.fill(0.0); Xu.fill(0.0); Xv.fill(0.0); Xuu.fill(0.0); Xuv.fill(0.0); Xvv.fill(0.0);
  F = 0.0; Fu = 0.0; Fv = 0.0; R = 0.0; Ru = 0.0; Rv = 0.0;

  G.g = G.gInv = 0.0;
  for(size_t a1 = 0; a1 < 2; a1++)
    for(size_t a2 = 0; a2 < 2; a2++)
      {
      G.xCovariantTensor[a1][a2] = 0.0;
      G.xContravariantTensor[a1][a2] = 0.0;
      for(size_t a3 = 0; a3 < 2; a3++)
        {
        G.xChristoffelFirst[a1][a2][a3] = 0.0;
        G.xChristoffelSecond[a1][a2][a3] = 0.0;
        }
      }

  aelt = 0.0;

  N.fill(0.0);
  xGradR.fill(0.0);
  xLapR = 0.0;
  xGradRMagSqr = 0.0; xNormalFactor = 0.0;

  xBnd[0].X.fill(0.0); xBnd[0].N.fill(0.0);
  xBnd[1].X.fill(0.0); xBnd[1].N.fill(0.0);
}

void
MedialAtom
::ComputeDifferentialGeometry()
{
  G.SetTwoJet( X.data_block(), Xu.data_block(), Xv.data_block(),
    Xuu.data_block(), Xuv.data_block(), Xvv.data_block());
}

void MedialAtom ::ComputeNormalVector()
{
  aelt = sqrt(G.g);
  N = vnl_cross_3d(Xu, Xv) / aelt;
}

bool MedialAtom::ComputeBoundaryAtoms(bool flagEdgeAtom)
{
  // Compute the Riemannian gradient of Phi
  double (*G2)[2] = G.xContravariantTensor;
  double Cu = (Fu * G2[0][0] + Fv * G2[0][1]);
  double Cv = (Fu * G2[0][1] + Fv * G2[1][1]);
  SMLVec3d xGradPhi = Xu * Cu + Xv * Cv;

  // Compute the radius of the atom
  R = sqrt(F);

  // Compute the Riemannian gradient of R
  xGradR = 0.5 * xGradPhi / R;

  // Compute the gradient magnitude of R
  xGradRMagSqr = 0.25 * (Fu * Cu + Fv * Cv) / F;

  // Split depending on whether this is an end atom
  if(flagEdgeAtom)
    {
    // We are at a crest, valid atom
    flagCrest = true;

    // There is a zero badness value
    // xGradRMagSqr = 1.0;
    xNormalFactor = 0.0;

    // Make grad-R have unit length
    // xGradR = xGradPhi / xGradPhi.magnitude(); 

    // If the magnitude of gradR is zero, we are in real trouble
    if(xGradRMagSqr == 0.0)
      {
      flagValid = false;
      xBnd[0].X = xBnd[1].X = X;
      xBnd[0].N = -N;
      xBnd[1].N = N;
      }
    else
      {
      flagValid = true;

      // Use a different computation for the boundary atoms. This
      // should be a correct computation if |gradR|==1 but it will
      // different for brute force models where gradR may not be 1
      // on the boundary
      xBnd[0].N = xBnd[1].N = - xGradR / sqrt(xGradRMagSqr);
      xBnd[0].X = xBnd[1].X = X + R * xBnd[0].N;
      }
    }
  else
    {
    // Otherwise, we have a normal atom
    flagCrest = false;

    // Compute the Riemannian gradient of R
    xGradR = 0.5 * xGradPhi / R;

    // Compute the gradient magnitude of R
    xGradRMagSqr = 0.25 * (Fu * Cu + Fv * Cv) / F;
    // xGradRMagSqr = xGradR.squared_magnitude();

    // Check if this is greater than one - should never be very close either
    if(xGradRMagSqr > 1.0)
      {
      // This is an invalid atom
      flagValid = false;

      // Set the boundary atoms to zero
      xNormalFactor = 0;
      xBnd[0].X = xBnd[1].X = X;
      xBnd[0].N = -N;
      xBnd[1].N = N;
      }
    else
      {
      // Finally, we have a clean internal atom
      flagValid = true;

      // Compute the normal component of the sail vectors
      xNormalFactor = sqrt(1.0 - xGradRMagSqr);
      SMLVec3d CN = N * xNormalFactor;

      // Compute the position and normals of boundary points
      xBnd[0].N = - xGradR - CN;
      xBnd[0].X = X + xBnd[0].N * R;
      xBnd[1].N = - xGradR + CN;
      xBnd[1].X = X + xBnd[1].N * R;
      }
    }

  return flagValid;
}


void
MedialAtom
::ComputeCommonDerivativeTerms(MedialAtom::DerivativeTerms &dt) const
{
  // Get the elements of the first fundamental form and its derivative
  const double &g11 = G.xContravariantTensor[0][0];
  const double &g12 = G.xContravariantTensor[0][1];
  const double &g22 = G.xContravariantTensor[1][1];

  // Get the derivatives of R
  dt.x1_2R = 0.5 / R;
  dt.Ru = Fu * dt.x1_2R, dt.Rv = Fv * dt.x1_2R;
  dt.Ru_R = dt.Ru / R;
  dt.Rv_R = dt.Rv / R;
  dt.Ru_2F = 0.5 * dt.Ru / F;
  dt.Rv_2F = 0.5 * dt.Rv / F;
  
  // Terms used to compute the derivative of the normal vector
  dt.Xu_aelt = Xu / aelt; dt.Xv_aelt = Xv / aelt;
  dt.N_2g = 0.5 * N / G.g;

  // We will compute several intermediate terms
  dt.g1iRi = g11 * dt.Ru + g12 * dt.Rv;
  dt.g2iRi = g12 * dt.Ru + g22 * dt.Rv;

  // Compute the plus-minus term
  dt.N_2nt = 0.5 * N / xNormalFactor;

  // Compute the inverse magnitude of grad R
  dt.inv_mag_gradR = (flagValid) ?  1.0 / xGradR.magnitude() : 0.0;
}

// Compute directional derivative of the contravariant tensor given the
// directional derivatives of X, Xu and Xv contained in dAtom
void
MedialAtom
::ComputeMetricTensorDerivatives(MedialAtom &dAtom) const
{
  // Compute the derivatives of the covariant tensor
  dAtom.G.xCovariantTensor[0][0] = 2 * dot(dAtom.Xu, Xu);
  dAtom.G.xCovariantTensor[1][1] = 2 * dot(dAtom.Xv, Xv);
  dAtom.G.xCovariantTensor[0][1] = 
    dAtom.G.xCovariantTensor[1][0] = dot(dAtom.Xv, Xu) + dot(dAtom.Xu, Xv);

  // Compute the derivative of g
  dAtom.G.g = 
    dAtom.G.xCovariantTensor[0][0] * G.xCovariantTensor[1][1] +
    dAtom.G.xCovariantTensor[1][1] * G.xCovariantTensor[0][0] -
    2 * dAtom.G.xCovariantTensor[0][1] * G.xCovariantTensor[0][1];

  // Compute the derivatives of the contravariant tensor
  dAtom.G.xContravariantTensor[0][0] = G.gInv *
    (dAtom.G.xCovariantTensor[1][1] - G.xContravariantTensor[0][0] * dAtom.G.g);

  dAtom.G.xContravariantTensor[1][1] = G.gInv *
    (dAtom.G.xCovariantTensor[0][0] - G.xContravariantTensor[1][1] * dAtom.G.g);

  dAtom.G.xContravariantTensor[0][1] = dAtom.G.xContravariantTensor[1][0] = - G.gInv *
    (dAtom.G.xCovariantTensor[0][1] + G.xContravariantTensor[0][1] * dAtom.G.g);

  // Compute the area element
  dAtom.aelt = 0.5 * dAtom.G.g / aelt;
}

// Prerequisites:
//  * The first jet for the atom and datom
//  * The contravariant tensor for the atom and datom
//
void 
MedialAtom::ComputeBoundaryAtomDerivatives(
  MedialAtom &dAtom, const DerivativeTerms &dt) const
{
  // Get the relevant elements of the atoms
  SMLVec3d &Y = dAtom.X, &Yu = dAtom.Xu, &Yv = dAtom.Xv;
  
  // Get the elements of the first fundamental form and its derivative
  const double &g11 = G.xContravariantTensor[0][0];
  const double &g12 = G.xContravariantTensor[0][1];
  const double &g22 = G.xContravariantTensor[1][1];
  const double &z11 = dAtom.G.xContravariantTensor[0][0];
  const double &z12 = dAtom.G.xContravariantTensor[0][1];
  const double &z22 = dAtom.G.xContravariantTensor[1][1];

  // Get the g's
  const double &g = G.g; const double &z = dAtom.G.g;

  // Get the partials of Phi and its variational derivative
  double &H = dAtom.F, &Hu = dAtom.Fu, &Hv = dAtom.Fv;

  // Get the derivatives of R
  double P = H * dt.x1_2R;
  double Pu = Hu * dt.x1_2R - H * dt.Ru_2F;
  double Pv = Hv * dt.x1_2R - H * dt.Rv_2F;

  dAtom.R = P; 
  
  // This is the derivative of the normal vector
  dAtom.N =  vnl_cross_3d(dt.Xu_aelt, Yv);
  dAtom.N += vnl_cross_3d(Yu, dt.Xv_aelt);
  vmuladd(dAtom.N, dt.N_2g, -z);

  // We will compute several intermediate terms
  double z1iRi = z11 * dt.Ru + z12 * dt.Rv;
  double z2iRi = z12 * dt.Ru + z22 * dt.Rv;

  // Compute the derivative of grad R
  vmulset(dAtom.xGradR, Xu, z1iRi + g11 * Pu + g12 * Pv);
  vmuladd(dAtom.xGradR, Xv, z2iRi + g12 * Pu + g22 * Pv);
  vmuladd(dAtom.xGradR, Yu, dt.g1iRi);
  vmuladd(dAtom.xGradR, Yv, dt.g2iRi);

  // Compute the derivative of grad R . grad R
  dAtom.xGradRMagSqr = z1iRi * dt.Ru + z2iRi * dt.Rv 
    + 2.0 * (dt.g1iRi * Pu + dt.g2iRi * Pv);

  // Address the edge case first
  if(flagCrest) 
    {
    dAtom.xNormalFactor = 0.0;
    if(flagValid)
      {
      // Compute (gradR / |gradR|)'
      vmulset(dAtom.xBnd[0].N, dAtom.xGradR, -dt.inv_mag_gradR);
      vmuladd(dAtom.xBnd[0].N, xBnd[0].N, 
        - dot_product(dAtom.xBnd[0].N, xBnd[0].N));
      dAtom.xBnd[1].N = dAtom.xBnd[0].N;

      // Compute the boundary atoms
      dAtom.xBnd[1].X = dAtom.xBnd[0].X = 
        Y + R * dAtom.xBnd[0].N + P * xBnd[0].N;

      // The normal term vanishes
      // dAtom.xBnd[0].N = dAtom.xBnd[1].N = - dAtom.xGradR;
      // dAtom.xGradRMagSqr = 0.0;
      }
    else
      {
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = dAtom.X;
      dAtom.xBnd[0].N = -dAtom.N;
      dAtom.xBnd[1].N = dAtom.N;
      }
    }
  else
    {
    if(flagValid)
      {
      // Compute the plus-minus term
      SMLVec3d dNormalTerm;
      vmulset(dNormalTerm, dAtom.N, xNormalFactor);
      vmuladd(dNormalTerm, dt.N_2nt, -dAtom.xGradRMagSqr);

      // Compute the boundary atom normals
      dAtom.xBnd[0].N = dAtom.xBnd[1].N = - dAtom.xGradR;
      dAtom.xBnd[0].N -= dNormalTerm;
      dAtom.xBnd[1].N += dNormalTerm;

      // Compute the boundary atoms
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = Y;
      vmuladd(dAtom.xBnd[0].X, dAtom.xBnd[0].N, R);
      vmuladd(dAtom.xBnd[0].X, xBnd[0].N, P);
      vmuladd(dAtom.xBnd[1].X, dAtom.xBnd[1].N, R);
      vmuladd(dAtom.xBnd[1].X, xBnd[1].N, P);
      }
    else
      {
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = dAtom.X;
      dAtom.xBnd[0].N = -dAtom.N;
      dAtom.xBnd[1].N = dAtom.N;
      }
    }
}




void AddScaleMedialAtoms(
  const MedialAtom &A, const MedialAtom &B, double p, MedialAtom &C)
{
  // The u and v coordinates stay the same
  C.u = A.u;
  C.v = A.v;
  C.uIndex = A.uIndex;
  C.vIndex = A.vIndex;
  C.flagCrest = A.flagCrest;
  C.flagValid = A.flagValid;

  // The other objects are simply added and scaled
  C.X = A.X + p * B.X;
  C.Xu = A.Xu + p * B.Xu;
  C.Xv = A.Xv + p * B.Xv;
  C.Xuu = A.Xuu + p * B.Xuu;
  C.Xuv = A.Xuv + p * B.Xuv;
  C.Xvv = A.Xvv + p * B.Xvv;

  C.R = A.R + p * B.R;
  C.Ru = A.Ru + p * B.Ru;
  C.Rv = A.Rv + p * B.Rv;
  // C.Ruu = A.Ruu + p * B.Ruu;
  // C.Ruv = A.Ruv + p * B.Ruv;
  // C.Rvv = A.Rvv + p * B.Rvv;

  C.F = A.F + p * B.F;
  C.Fu = A.Fu + p * B.Fu;
  C.Fv = A.Fv + p * B.Fv;

  C.N = A.N + p * B.N;
  C.xGradR = A.xGradR + p * B.xGradR;
  // C.xGradPhi = A.xGradPhi + p * B.xGradPhi;
  C.xGradRMagSqr = A.xGradRMagSqr + p * B.xGradRMagSqr;
  C.xNormalFactor = A.xNormalFactor + p * B.xNormalFactor;

  // The differential geometry is also added and scaled 
  C.G.g = A.G.g + p * B.G.g;
  C.G.gInv = A.G.gInv + p * B.G.gInv;

  for(size_t i=0; i<2; i++) for(size_t j=0; j<2; j++)
    {
    C.G.xContravariantTensor[i][j] = 
      A.G.xContravariantTensor[i][j] + p * B.G.xContravariantTensor[i][j];
    C.G.xCovariantTensor[i][j] = 
      A.G.xCovariantTensor[i][j] + p * B.G.xCovariantTensor[i][j];
    C.G.xChristoffelSecond[i][j][0] = 
      A.G.xChristoffelSecond[i][j][0] + p * B.G.xChristoffelSecond[i][j][0];
    C.G.xChristoffelFirst[i][j][0] = 
      A.G.xChristoffelFirst[i][j][0] + p * B.G.xChristoffelFirst[i][j][0];
    C.G.xChristoffelSecond[i][j][1] = 
      A.G.xChristoffelSecond[i][j][1] + p * B.G.xChristoffelSecond[i][j][1];
    C.G.xChristoffelFirst[i][j][1] = 
      A.G.xChristoffelFirst[i][j][1] + p * B.G.xChristoffelFirst[i][j][1];
    }

  // The boundary atoms are scaled as well
  for(size_t k=0; k<2; k++)
    {
    C.xBnd[k].X = A.xBnd[k].X + p * B.xBnd[k].X;
    C.xBnd[k].N = A.xBnd[k].N + p * B.xBnd[k].N;
    }
}

void MedialAtomCentralDifference(
  const MedialAtom &A, const MedialAtom &B, double eps, MedialAtom &C)
{
  double p = 0.5 / eps;
  
  // The u and v coordinates stay the same
  C.u = A.u;
  C.v = A.v;
  C.uIndex = A.uIndex;
  C.vIndex = A.vIndex;
  C.flagCrest = A.flagCrest;
  C.flagValid = A.flagValid;

  // The other objects are simply added and scaled
  C.X = p * (A.X - B.X);
  C.Xu = p * (A.Xu - B.Xu);
  C.Xv = p * (A.Xv - B.Xv);
  C.Xuu = p * (A.Xuu - B.Xuu);
  C.Xuv = p * (A.Xuv - B.Xuv);
  C.Xvv = p * (A.Xvv - B.Xvv);

  C.R = p * (A.R - B.R);
  C.Ru = p * (A.Ru - B.Ru);
  C.Rv = p * (A.Rv - B.Rv);
  // C.Ruu = A.Ruu + p * B.Ruu;
  // C.Ruv = A.Ruv + p * B.Ruv;
  // C.Rvv = A.Rvv + p * B.Rvv;

  C.F = p * (A.F - B.F);
  C.Fu = p * (A.Fu - B.Fu);
  C.Fv = p * (A.Fv - B.Fv);

  C.N = p * (A.N - B.N);
  C.xGradR = p * (A.xGradR - B.xGradR);
  // C.xGradPhi = p * (A.xGradPhi - B.xGradPhi);
  C.xGradRMagSqr = p * (A.xGradRMagSqr - B.xGradRMagSqr);
  C.xNormalFactor = p * (A.xNormalFactor - B.xNormalFactor);

  // The differential geometry is also added and scaled 
  C.G.g = p * (A.G.g - B.G.g);
  C.G.gInv = p * (A.G.gInv - B.G.gInv);

  for(size_t i=0; i<2; i++) for(size_t j=0; j<2; j++)
    {
    C.G.xContravariantTensor[i][j] = 
      p * (A.G.xContravariantTensor[i][j] - B.G.xContravariantTensor[i][j]);
    C.G.xCovariantTensor[i][j] = 
      p * (A.G.xCovariantTensor[i][j] - B.G.xCovariantTensor[i][j]);
    C.G.xChristoffelSecond[i][j][0] = 
      p * (A.G.xChristoffelSecond[i][j][0] - B.G.xChristoffelSecond[i][j][0]);
    C.G.xChristoffelFirst[i][j][0] = 
      p * (A.G.xChristoffelFirst[i][j][0] - B.G.xChristoffelFirst[i][j][0]);
    C.G.xChristoffelSecond[i][j][1] = 
      p * (A.G.xChristoffelSecond[i][j][1] - B.G.xChristoffelSecond[i][j][1]);
    C.G.xChristoffelFirst[i][j][1] = 
      p * (A.G.xChristoffelFirst[i][j][1] - B.G.xChristoffelFirst[i][j][1]);
    }

  // The boundary atoms are scaled as well
  for(size_t k=0; k<2; k++)
    {
    C.xBnd[k].X = p * (A.xBnd[k].X - B.xBnd[k].X);
    C.xBnd[k].N = p * (A.xBnd[k].N - B.xBnd[k].N);
    }
}
