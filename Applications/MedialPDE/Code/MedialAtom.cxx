#include "MedialAtom.h"

void
MedialAtom
::ComputeDifferentialGeometry()
{
  G.SetJet( X.data_block(), Xu.data_block(), Xv.data_block(),
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
  xGradPhi = Xu * Cu + Xv * Cv;

  return ComputeBoundaryAtoms(xGradPhi, flagEdgeAtom);
}

bool MedialAtom::ComputeBoundaryAtoms(const SMLVec3d &xGradPhi, bool flagEdgeAtom)
{
  // Store the grad phi
  this->xGradPhi = xGradPhi;

  // Compute the radius of the atom
  R = sqrt(F);

  // Compute the Riemannian gradients of Phi and R
  xGradR = 0.5 * xGradPhi / R;

  // Split depending on whether this is an end atom
  if(flagEdgeAtom)
    {
    // We are at a crest, valid atom
    flagValid = true; flagCrest = true;

    // There is a zero badness value
    xGradRMagSqr = 1.0;
    xNormalFactor = 0.0;

    // Simpler geometry, save a square root!
    xBnd[0].N = xBnd[1].N = -xGradR;
    xBnd[0].X = xBnd[1].X = X - xGradR * R;

    }
  else
    {
    // Otherwise, we have a normal atom
    flagCrest = false;

    // Compute the gradient magnitude of R
    // xGradRMagSqr = 0.25 * (Fu * Cu + Fv * Cv) / F;
    xGradRMagSqr = xGradR.squared_magnitude();

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

  // Compute the partials of R
  // R = sqrt(F);
  // Ru = 0.5 * Fu / R;
  // Rv = 0.5 * Fv / R;

  return flagValid;
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
  // C.Ru = A.Ru + p * B.Ru;
  // C.Rv = A.Rv + p * B.Rv;
  // C.Ruu = A.Ruu + p * B.Ruu;
  // C.Ruv = A.Ruv + p * B.Ruv;
  // C.Rvv = A.Rvv + p * B.Rvv;

  C.F = A.F + p * B.F;
  C.Fu = A.Fu + p * B.Fu;
  C.Fv = A.Fv + p * B.Fv;

  C.N = A.N + p * B.N;
  C.xGradR = A.xGradR + p * B.xGradR;
  C.xGradPhi = A.xGradPhi + p * B.xGradPhi;
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
  // C.Ru = A.Ru + p * B.Ru;
  // C.Rv = A.Rv + p * B.Rv;
  // C.Ruu = A.Ruu + p * B.Ruu;
  // C.Ruv = A.Ruv + p * B.Ruv;
  // C.Rvv = A.Rvv + p * B.Rvv;

  C.F = p * (A.F - B.F);
  C.Fu = p * (A.Fu - B.Fu);
  C.Fv = p * (A.Fv - B.Fv);

  C.N = p * (A.N - B.N);
  C.xGradR = p * (A.xGradR - B.xGradR);
  C.xGradPhi = p * (A.xGradPhi - B.xGradPhi);
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
