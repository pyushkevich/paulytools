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
  N = vnl_cross_3d(Xu, Xv) * sqrt(G.gInv);
}

bool MedialAtom::ComputeBoundaryAtoms(bool flagEdgeAtom)
{
  // Get a shorthand for the contravariant tensor
  double (*G2)[2] = G.xContravariantTensor;

  // Compute the radius of the atom
  R = sqrt(F);

  // Compute the Riemannian gradients of Phi and R
  double Cu = (Fu * G2[0][0] + Fv * G2[0][1]);
  double Cv = (Fu * G2[0][1] + Fv * G2[1][1]);
  xGradPhi = Xu * Cu + Xv * Cv;
  xGradR = 0.5 * xGradPhi / R;

  // Split depending on whether this is an end atom
  if(flagEdgeAtom)
    {
    // We are at a crest, valid atom
    flagValid = true; flagCrest = true;

    // There is a zero badness value
    xGradRMagSqr = 1.0;

    // Simpler geometry, save a square root!
    xBnd[0].N = xBnd[1].N = -xGradR;
    xBnd[0].X = xBnd[1].X = X - xGradR * R;

    }
  else
    {
    // Otherwise, we have a normal atom
    flagCrest = false;

    // Compute the gradient magnitude of R
    xGradRMagSqr = 0.25 * (Fu * Cu + Fv * Cv) / F;

    // Check if this is greater than one - should never be very close either
    if(xGradRMagSqr > 1.0)
      {
      // This is an invalid atom
      flagValid = false;

      // Set the boundary atoms to zero
      xBnd[0].X = xBnd[1].X = X;
      xBnd[0].N = -N;
      xBnd[1].N = N;
      }
    else
      {
      // Finally, we have a clean internal atom
      flagValid = true;

      // Compute the normal component of the sail vectors
      SMLVec3d CN = N * sqrt(1.0 - xGradRMagSqr);

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
/*****************************************************************************
 * CELL VOLUME AND DERIVATIVE COMPUTATION
 ****************************************************************************/

inline double ScalarTripleProduct(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return
    A[0] * (B[1] * C[2] - B[2] * C[1]) + 
    A[1] * (B[2] * C[0] - B[0] * C[2]) + 
    A[2] * (B[0] * C[1] - B[1] * C[0]);
}

inline double ScalarTripleProductDerivative(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C,
  const SMLVec3d &DA, const SMLVec3d &DB, const SMLVec3d &DC)
{
  return
    DA[0] * (B[1] * C[2] - B[2] * C[1]) + 
    DA[1] * (B[2] * C[0] - B[0] * C[2]) + 
    DA[2] * (B[0] * C[1] - B[1] * C[0]) +
    A[0] * (DB[1] * C[2] - DB[2] * C[1] + B[1] * DC[2] - B[2] * DC[1]) + 
    A[1] * (DB[2] * C[0] - DB[0] * C[2] + B[2] * DC[0] - B[0] * DC[2]) + 
    A[2] * (DB[0] * C[1] - DB[1] * C[0] + B[0] * DC[1] - B[1] * DC[0]);
}


// This returns 6 x the volume of the tetrahedron
inline double TetrahedronVolume(
  const SMLVec3d &X1, const SMLVec3d &X2, const SMLVec3d &X3, const SMLVec3d X4)
{
  return dot_product(vnl_cross_3d(X2 - X1, X3 - X1), X4 - X1);
}

// This returns 6 x the derivative
inline double TetrahedronVolumePartialDerivative(
  const SMLVec3d &X1, const SMLVec3d &X2, const SMLVec3d &X3, const SMLVec3d X4,
  const SMLVec3d &D1, const SMLVec3d &D2, const SMLVec3d &D3, const SMLVec3d D4)
{
  return 
    dot_product(
      vnl_cross_3d( X1 - X3, D2 ) + 
      vnl_cross_3d( X2 - X1, D3 ) + 
      vnl_cross_3d( X3 - X2, D1 ) , X4 - X1 ) + 
    dot_product( 
      vnl_cross_3d( X2 - X1, X3 - X1 ) , D4 - D1 );
}

/**
 * 
double CellVolume(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111)
{
  // Compute the centerpoint inside the cell
  SMLVec3d C = (X000 + X001 + X010 + X011 + X100 + X101 + X110 + X111) * 0.125;

  // Compute each of the twelve tetrahedra inside the cuboid (group by the
  // face of the cell)
  double v00 = TetrahedronVolume(X000, X001, X100, C);
  double v01 = TetrahedronVolume(X101, X100, X001, C);

  double v10 = TetrahedronVolume(X100, X101, X110, C);
  double v11 = - TetrahedronVolume(X111, X110, X101, C);
  
  double v20 = - TetrahedronVolume(X110, X111, X010, C);
  double v21 = TetrahedronVolume(X011, X010, X111, C);
  
  double v30 = TetrahedronVolume(X010, X011, X000, C);
  double v31 = TetrahedronVolume(X001, X000, X011, C);
  
  double v40 = TetrahedronVolume(X011, X111, X001, C);
  double v41 = - TetrahedronVolume(X101, X001, X111, C);
  
  double v50 = - TetrahedronVolume(X010, X000, X110, C);
  double v51 = TetrahedronVolume(X100, X110, X000, C);

   cout << "VOLS: " 
    << (v00 > 0 ? "+1" : "-1") << ", " 
    << (v01 > 0 ? "+1" : "-1") << ", " 
    << (v10 > 0 ? "+1" : "-1") << ", " 
    << (v11 > 0 ? "+1" : "-1") << ", " 
    << (v20 > 0 ? "+1" : "-1") << ", " 
    << (v21 > 0 ? "+1" : "-1") << ", " 
    << (v30 > 0 ? "+1" : "-1") << ", " 
    << (v31 > 0 ? "+1" : "-1") << ", " 
    << (v40 > 0 ? "+1" : "-1") << ", " 
    << (v41 > 0 ? "+1" : "-1") << ", " 
    << (v50 > 0 ? "+1" : "-1") << ", " 
    << (v51 > 0 ? "+1" : "-1") << endl;
  
  // Return a sixth of these sums
  static const double SIXTH = 1.0 / 6.0;
  return - SIXTH *
    (v00 + v01 + v10 + v11 + v20 + v21 + v30 + v31 + v40 + v41 + v50 + v51);
}
*/

/**
 * This in another, less precise, formula for the volume element. We just
 * model the solid as a parallellepiped with sides computed as averages
 */
double CellVolume(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111)
{
  SMLVec3d U = X111 + X110 + X101 + X100 - X011 - X010 - X001 - X000;
  SMLVec3d V = X111 + X110 - X101 - X100 + X011 + X010 - X001 - X000;
  SMLVec3d W = X111 - X110 + X101 - X100 + X011 - X010 + X001 - X000;
  return - 0.015625 * ScalarTripleProduct(U, V, W);
} 

double CellVolumePartialDerivative(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111,
  const SMLVec3d &D000, const SMLVec3d &D001, 
  const SMLVec3d &D010, const SMLVec3d &D011, 
  const SMLVec3d &D100, const SMLVec3d &D101, 
  const SMLVec3d &D110, const SMLVec3d &D111)
{
  SMLVec3d U  = X111 + X110 + X101 + X100 - X011 - X010 - X001 - X000;
  SMLVec3d V  = X111 + X110 - X101 - X100 + X011 + X010 - X001 - X000;
  SMLVec3d W  = X111 - X110 + X101 - X100 + X011 - X010 + X001 - X000;
  SMLVec3d dU = D111 + D110 + D101 + D100 - D011 - D010 - D001 - D000;
  SMLVec3d dV = D111 + D110 - D101 - D100 + D011 + D010 - D001 - D000;
  SMLVec3d dW = D111 - D110 + D101 - D100 + D011 - D010 + D001 - D000;
  return - 0.015625 * ScalarTripleProductDerivative(U, V, W, dU, dV, dW);
}

/* 
double CellVolumePartialDerivative(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111,
  const SMLVec3d &D000, const SMLVec3d &D001, 
  const SMLVec3d &D010, const SMLVec3d &D011, 
  const SMLVec3d &D100, const SMLVec3d &D101, 
  const SMLVec3d &D110, const SMLVec3d &D111)
{
  // Compute the centerpoint inside the cell
  SMLVec3d C = (X000 + X001 + X010 + X011 + X100 + X101 + X110 + X111) * 0.125;
  SMLVec3d DC = (D000 + D001 + D010 + D011 + D100 + D101 + D110 + D111) * 0.125;

  // Compute each of the twelve tetrahedra inside the cuboid (group by the
  // face of the cell)
  double v00 = TetrahedronVolumePartialDerivative(
    X000, X001, X100, C, D000, D001, D100, DC);
  double v01 = TetrahedronVolumePartialDerivative(
    X101, X100, X001, C, D101, D100, D001, DC);

  double v10 = TetrahedronVolumePartialDerivative(
    X100, X101, X110, C, D100, D101, D110, DC);
  double v11 = - TetrahedronVolumePartialDerivative(
    X111, X110, X101, C, D111, D110, D101, DC);
  
  double v20 = - TetrahedronVolumePartialDerivative(
    X110, X111, X010, C, D110, D111, D010, DC);
  double v21 = TetrahedronVolumePartialDerivative(
    X011, X010, X111, C, D011, D010, D111, DC);
  
  double v30 = TetrahedronVolumePartialDerivative(
    X010, X011, X000, C, D010, D011, D000, DC);
  double v31 = TetrahedronVolumePartialDerivative(
    X001, X000, X011, C, D001, D000, D011, DC);
  
  double v40 = TetrahedronVolumePartialDerivative(
    X011, X111, X001, C, D011, D111, D001, DC);
  double v41 = - TetrahedronVolumePartialDerivative(
    X101, X001, X111, C, D101, D001, D111, DC);
  
  double v50 = - TetrahedronVolumePartialDerivative(
    X010, X000, X110, C, D010, D000, D110, DC);
  double v51 = TetrahedronVolumePartialDerivative(
    X100, X110, X000, C, D100, D110, D000, DC);

  // Return a sixth of these sums
  static const double SIXTH = 1.0 / 6.0;
  return - SIXTH *
    (v00 + v01 + v10 + v11 + v20 + v21 + v30 + v31 + v40 + v41 + v50 + v51);
}
*/

/** Helper function to compute the coordinate of an internal medial point
 * that is referenced by an iterator */
SMLVec3d GetInternalPoint(MedialInternalPointIterator *it, MedialAtom *xAtoms)
{
  // If a medial atom, return the atom
  size_t iAtom = it->GetAtomIndex();
  size_t iDepth = it->GetDepth();

  // If depth is zero (we are on the medial axis) return the medial point
  if(iDepth == 0) return xAtoms[iAtom].X;

  // If depth is equal to max-depth, we can return the boundary site instead
  SMLVec3d Y = xAtoms[iAtom].xBnd[it->GetBoundarySide()].X;
  if(iDepth == it->GetMaxDepth()) return Y;

  // In the remaining case, the point is internal and should be interpolated
  SMLVec3d X = xAtoms[iAtom].X;
  return X + (Y - X) * it->GetRelativeDistanceToMedialAxis();
}


/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
/* 
double ComputeMedialBoundaryAreaWeights( MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, double *xWeights)
{
  // Clear the content of the weights
  memset(xWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  double xTotalArea = 0.0f;
  
  // Create a quad-based iterator
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();

  // For each quad, compute the area associated with it
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Access the four medial points
    const BoundaryAtom& A00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0);
    const BoundaryAtom& A01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1);
    const BoundaryAtom& A11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1);
    const BoundaryAtom& A10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0);

    // Compute the four vectors along the edges of this rectangle
    SMLVec3d U0 = A10.X - A00.X; SMLVec3d V0 = A01.X - A00.X;
    SMLVec3d U1 = A11.X - A01.X; SMLVec3d V1 = A11.X - A10.X;

    // Compute the contribution wo the weight elements of each node
    double W00 = -0.25 * ScalarTripleProduct(U0, V0, A00.N);
    double W10 = -0.25 * ScalarTripleProduct(U0, V1, A10.N);
    double W01 = -0.25 * ScalarTripleProduct(U1, V0, A01.N);
    double W11 = -0.25 * ScalarTripleProduct(U1, V1, A11.N);

    // Update the total area
    xTotalArea += W00 + W10 + W01 + W11;
    
    // Assign a third of each weight to each corner
    xWeights[itQuad->GetBoundaryIndex(0, 0)] += W00;
    xWeights[itQuad->GetBoundaryIndex(1, 1)] += W11;
    xWeights[itQuad->GetBoundaryIndex(0, 1)] += W01;
    xWeights[itQuad->GetBoundaryIndex(1, 0)] += W10;
    }

  // Delete the iterator object
  delete itQuad;

  // Return the total area
  return xTotalArea;
}
*/

// Old Method: always returns positive areas - dangerous!
double ComputeMedialBoundaryAreaWeights( MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, double *xWeights)
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  double xTotalArea = 0.0f;

  // Clear the content of the weights
  memset(xWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  
  // Create a quad-based iterator
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();

  // For each quad, compute the area associated with it
  size_t iQuad = 0;
  while(!itQuad->IsAtEnd())
    {
    // Access the four medial points
    SMLVec3d X00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0).X;
    SMLVec3d X01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1).X;
    SMLVec3d X11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1).X;
    SMLVec3d X10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0).X;

    // Compute the areas of the two triangles involved
    double A1 = THIRD * TriangleArea( X00, X01, X11 );
    double A2 = THIRD * TriangleArea( X00, X11, X10 );
    xTotalArea += (A1 + A2);

    // Assign a third of each weight to each corner
    xWeights[itQuad->GetBoundaryIndex(0, 0)] += A1 + A2;
    xWeights[itQuad->GetBoundaryIndex(1, 1)] += A1 + A2;
    xWeights[itQuad->GetBoundaryIndex(0, 1)] += A1;
    xWeights[itQuad->GetBoundaryIndex(1, 0)] += A2;

    // Go to the next boundary quad
    ++(*itQuad);
    }

  // Delete the iterator object
  delete itQuad;

  // Return the total area
  return 3.0f * xTotalArea;
}

/**
 * Compute the derivative of the boundary area elements given the derivative
 * of the atoms themselves with respect to some parameter 
 */
/*
double ComputeMedialBoundaryAreaPartialDerivative(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights)
{
  // Clear the content of the weights
  memset(dWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  double dTotalArea = 0.0f;
  
  // Create a quad-based iterator
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();

  // For each quad, compute the area associated with it
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Access the four medial points
    const BoundaryAtom& A00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0);
    const BoundaryAtom& A01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1);
    const BoundaryAtom& A11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1);
    const BoundaryAtom& A10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0);

    // Also access their derivatives
    const BoundaryAtom& DA00 = GetBoundaryPoint(itQuad, dAtoms, 0, 0);
    const BoundaryAtom& DA01 = GetBoundaryPoint(itQuad, dAtoms, 0, 1);
    const BoundaryAtom& DA11 = GetBoundaryPoint(itQuad, dAtoms, 1, 1);
    const BoundaryAtom& DA10 = GetBoundaryPoint(itQuad, dAtoms, 1, 0);

    // Compute the four vectors along the edges of this rectangle
    SMLVec3d U0 = A10.X - A00.X; SMLVec3d V0 = A01.X - A00.X;
    SMLVec3d U1 = A11.X - A01.X; SMLVec3d V1 = A11.X - A10.X;
    SMLVec3d DU0 = DA10.X - DA00.X; SMLVec3d DV0 = DA01.X - DA00.X;
    SMLVec3d DU1 = DA11.X - DA01.X; SMLVec3d DV1 = DA11.X - DA10.X;

    // Compute the contribution wo the weight elements of each node
    double DW00 = -0.25 * ScalarTripleProductDerivative(U0, V0, A00.N, DU0, DV0, DA00.N);
    double DW10 = -0.25 * ScalarTripleProductDerivative(U0, V1, A10.N, DU0, DV1, DA10.N);
    double DW01 = -0.25 * ScalarTripleProductDerivative(U1, V0, A01.N, DU1, DV0, DA01.N);
    double DW11 = -0.25 * ScalarTripleProductDerivative(U1, V1, A11.N, DU1, DV1, DA11.N);

    // Update the total area
    dTotalArea += DW00 + DW10 + DW01 + DW11;
    
    // Assign a third of each weight to each corner
    dWeights[itQuad->GetBoundaryIndex(0, 0)] += DW00;
    dWeights[itQuad->GetBoundaryIndex(1, 1)] += DW11;
    dWeights[itQuad->GetBoundaryIndex(0, 1)] += DW01;
    dWeights[itQuad->GetBoundaryIndex(1, 0)] += DW10;
    }

  // Delete the iterator object
  delete itQuad;

  // Return the total area
  return dTotalArea;
}
*/

double ComputeMedialBoundaryAreaPartialDerivative(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights)
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  double dTotalArea = 0.0f;

  // Clear the content of the weights
  memset(dWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  
  // Create a quad-based iterator
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();

  // For each quad, compute the area associated with it
  size_t iQuad = 0;
  while(!itQuad->IsAtEnd())
    {
    // Access the four medial points
    SMLVec3d X00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0).X;
    SMLVec3d X01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1).X;
    SMLVec3d X11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1).X;
    SMLVec3d X10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0).X;

    // Access the derivatives of the points
    SMLVec3d D00 = GetBoundaryPoint(itQuad, dAtoms, 0, 0).X;
    SMLVec3d D01 = GetBoundaryPoint(itQuad, dAtoms, 0, 1).X;
    SMLVec3d D11 = GetBoundaryPoint(itQuad, dAtoms, 1, 1).X;
    SMLVec3d D10 = GetBoundaryPoint(itQuad, dAtoms, 1, 0).X;

    // Get the areas of the two triangles
    double A1 = xWeights[itQuad->GetBoundaryIndex(0, 1)];
    double A2 = xWeights[itQuad->GetBoundaryIndex(1, 0)];

    // Compute the areas of the two triangles involved
    double dA1 = THIRD * 
      TriangleAreaPartialDerivative( X00, X01, X11, D00, D01, D11, 
        TriangleArea(X00,X01,X11));
    
    double dA2 = THIRD * 
      TriangleAreaPartialDerivative( X00, X11, X10, D00, D11, D10, 
        TriangleArea(X00, X11, X10));

    // Update the total area
    dTotalArea += dA1 + dA2;

    // Assign a third of each weight to each corner
    dWeights[itQuad->GetBoundaryIndex(0, 0)] += dA1 + dA2;
    dWeights[itQuad->GetBoundaryIndex(1, 1)] += dA1 + dA2;
    dWeights[itQuad->GetBoundaryIndex(0, 1)] += dA1;
    dWeights[itQuad->GetBoundaryIndex(1, 0)] += dA2;

    // Go to the next boundary quad
    ++(*itQuad);
    }

  // Delete the iterator object
  delete itQuad;

  // Return the total area
  return 3.0f * dTotalArea;
}

/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints)
{
  // Create an internal point iterator
  MedialInternalPointIterator *it = xGrid->NewInternalPointIterator(nCuts);
  while(!it->IsAtEnd())
    {
    // Compute the internal point
    xPoints[it->GetIndex()] = GetInternalPoint(it, xAtoms);
    ++(*it);
    }
  delete it;
}

class CellVolumeWeightComputer {
public:
  
  CellVolumeWeightComputer(SMLVec3d *xPoints)
    { this->xPoints = xPoints; }
  
  double ComputeWeight (
    size_t i000, size_t i001, size_t i010, size_t i011, 
    size_t i100, size_t i101, size_t i110, size_t i111)
    {
    return CellVolume(
      xPoints[i000], xPoints[i001], xPoints[i010], xPoints[i011], 
      xPoints[i100], xPoints[i101], xPoints[i110], xPoints[i111]); 
    }

private:
  SMLVec3d *xPoints;
};

class CellVolumePartialDerivativeWeightComputer {
public:
  
  CellVolumePartialDerivativeWeightComputer(SMLVec3d *xPoints, SMLVec3d *dPoints)
    { this->xPoints = xPoints; this->dPoints = dPoints; }
  
  double ComputeWeight (
    size_t i000, size_t i001, size_t i010, size_t i011, 
    size_t i100, size_t i101, size_t i110, size_t i111)
    {
    return CellVolumePartialDerivative (
      xPoints[i000], xPoints[i001], xPoints[i010], xPoints[i011], 
      xPoints[i100], xPoints[i101], xPoints[i110], xPoints[i111],
      dPoints[i000], dPoints[i001], dPoints[i010], dPoints[i011], 
      dPoints[i100], dPoints[i101], dPoints[i110], dPoints[i111]); 
    }

private:
  SMLVec3d *xPoints, *dPoints;
};

/**
 * This generic function reuses code to compute both the cell volume-based
 * weights and their derivatives.
 */
template<class TWeightComputer>
double ComputeMedialInteralWeights(
  MedialAtomGrid *xGrid, TWeightComputer *xComputer, size_t nCuts, 
  double *xWeights, double *xInternalProfileWeights)
{
  // Set all the weights to zero to begin with
  size_t nPoints = xGrid->GetNumberOfInternalPoints(nCuts);
  memset(xWeights, 0, sizeof(double) * nPoints);
  memset(xInternalProfileWeights, 0, sizeof(double) 
    * xGrid->GetNumberOfProfileIntervals(nCuts));
  
  // Create an internal cell iterator
  MedialInternalCellIterator *itCell = xGrid->NewInternalCellIterator(nCuts);

  // For each cell, compute the internal volume
  double xTotalVolume = 0.0;
  while(!itCell->IsAtEnd())
    {
    // Get the indices of the cell corners
    size_t i000 = itCell->GetInternalPointIndex(0, 0, 0);
    size_t i001 = itCell->GetInternalPointIndex(0, 0, 1);
    size_t i010 = itCell->GetInternalPointIndex(0, 1, 0);
    size_t i100 = itCell->GetInternalPointIndex(1, 0, 0);
    size_t i011 = itCell->GetInternalPointIndex(0, 1, 1);
    size_t i101 = itCell->GetInternalPointIndex(1, 0, 1);
    size_t i110 = itCell->GetInternalPointIndex(1, 1, 0);
    size_t i111 = itCell->GetInternalPointIndex(1, 1, 1);

    // Compute the cell's volume
    // cout.width(16);
    // cout << "Cell : " << itCell->GetAtomIndex(0,0) << " ";
    double V =  
      xComputer->ComputeWeight( i000, i001, i010, i011, i100, i101, i110, i111); 

    // Add to total volume of the structure
    xTotalVolume += V;

    // Scale by an eigth (since each cube has eight corners)
    V *= 0.125;

    // Assign the fractions of the volume to the corresponding internal points
    xWeights[i000] += V; xWeights[i110] += V;
    xWeights[i010] += V; xWeights[i100] += V;
    xWeights[i001] += V; xWeights[i111] += V;
    xWeights[i011] += V; xWeights[i101] += V;

    // Also assign fractions of the weights to the profile intervals
    if(xInternalProfileWeights)
      {
      size_t j00 = itCell->GetProfileIntervalIndex(0, 0);
      size_t j01 = itCell->GetProfileIntervalIndex(0, 1);
      size_t j10 = itCell->GetProfileIntervalIndex(1, 0);
      size_t j11 = itCell->GetProfileIntervalIndex(1, 1);
      V *= 2.0;
      xInternalProfileWeights[j00] += V;
      xInternalProfileWeights[j01] += V;
      xInternalProfileWeights[j10] += V;
      xInternalProfileWeights[j11] += V;
      }

    // Go to the next boundary quad
    ++(*itCell);
    }

  // Delete the iterator object
  delete itCell;

  // Return the total area
  return xTotalVolume; 
}


/**
 * This function computes the volume associated with each cell in the medial
 * interior and assigns a portion of this cell to each internal point that is
 * the vertex of the cell. The sum of these resulting internal point weights
 * is equal to the volume of the m-rep
 */
double ComputeMedialInternalVolumeWeights(
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, size_t nCuts, 
  double *xWeights, double *xInternalProfileWeights)
{
  // Create a specialization object
  CellVolumeWeightComputer xComputer(xPoints);

  // Compute the weights
  return ComputeMedialInteralWeights(xGrid, &xComputer, nCuts, 
    xWeights, xInternalProfileWeights);
}

double ComputeMedialInternalVolumePartialDerivativeWeights(
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, SMLVec3d *dPoints, size_t nCuts, 
  double *dWeights, double *dInternalProfileWeights)
{
  // Create a specialization object
  CellVolumePartialDerivativeWeightComputer xComputer(xPoints, dPoints);

  // Compute the weights
  return ComputeMedialInteralWeights(xGrid, &xComputer, nCuts, 
    dWeights, dInternalProfileWeights);
}

/** This method integrates any given function over the interior of mrep */
double IntegrateFunctionOverInterior (
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, 
  double *xWeights, size_t nCuts, EuclideanFunction *fMatch) 
{
  // Match accumulator
  double xMatch = 0.0;
  
  // Create an internal point iterator
  MedialInternalPointIterator *it = xGrid->NewInternalPointIterator(nCuts);
  while(!it->IsAtEnd())
    {
    // Evaluate the image at this location
    size_t iPoint = it->GetIndex();
    xMatch += xWeights[iPoint] * fMatch->Evaluate( xPoints[iPoint] );

    // On to the next point
    ++(*it); 
    }

  // Clean up
  delete it;

  // Return match scaled by total weight
  return xMatch;
}

/** Integrate any function over the interior of a medial atom */
double IntegrateFunctionOverBoundary (
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch)
{
  // Match accumulator
  double xMatch = 0.0;

  // Iterate through all boundary points
  MedialBoundaryPointIterator *itBoundary = xGrid->NewBoundaryPointIterator();
  while(!itBoundary->IsAtEnd())
    {
    double xLocal = fMatch->Evaluate(GetBoundaryPoint(itBoundary, xAtoms).X);
    xMatch += xLocal * xWeights[itBoundary->GetIndex()];
    ++(*itBoundary);
    }
  
  delete itBoundary;
  return xMatch;
}

double RecursiveGetTriangleMatch(
  SMLVec3d &A, SMLVec3d &B, SMLVec3d &C, 
  double xMinArea, EuclideanFunction *fMatch, double &xArea)
{
  // Compute the area of this triangle
  xArea = TriangleArea(A, B, C);

  // If the area is smaller or equal to the min area, return the current match
  if(xArea <= xMinArea)
    return 0.33333333333334 * 
      (fMatch->Evaluate(A) + fMatch->Evaluate(B) + fMatch->Evaluate(C));

  // Otherwise, subdivide the triangle into three points, and recurse for each
  // of the triangles
  double m[4], a[4];
  SMLVec3d AB = 0.5 * (A + B);
  SMLVec3d AC = 0.5 * (A + C);
  SMLVec3d BC = 0.5 * (B + C);

  m[0] = RecursiveGetTriangleMatch(A, AB, AC, xMinArea, fMatch, a[0]);
  m[1] = RecursiveGetTriangleMatch(AB, B, BC, xMinArea, fMatch, a[1]);
  m[2] = RecursiveGetTriangleMatch(AC, BC, C, xMinArea, fMatch, a[2]);
  m[3] = RecursiveGetTriangleMatch(AC, AB, BC, xMinArea, fMatch, a[3]);

  xArea = a[0] + a[1] + a[2] + a[3];
  return (a[0]*m[0] + a[1]*m[1] + a[2]*m[2] + a[3]*m[3]) / xArea;
}


/** 
 * Integrate a function over the boundary, making sure that irregular quads are 
 * super-sampled. This results in better quality integration 
 */
double AdaptivelyIntegrateFunctionOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double xMinQuadArea, EuclideanFunction *fMatch)
{
  // Match
  double xMatch = 0.0, xArea = 0.0;

  // Integrate the match over each quad
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Get the vector at each vertex
    SMLVec3d x00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0).X;
    SMLVec3d x01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1).X;
    SMLVec3d x11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1).X;
    SMLVec3d x10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0).X;

    // Call the recursive procedure for each sub-quad
    double A1, A2;
    double M1 = RecursiveGetTriangleMatch( x00, x01, x11, xMinQuadArea, fMatch, A1);
    double M2 = RecursiveGetTriangleMatch( x00, x11, x10, xMinQuadArea, fMatch, A2);

    // Add the weighted area to the match
    xMatch += A1 * M1 + A2 * M2;
    xArea += A1 + A2;
    }

  delete itQuad;

  // Scale the match by the area
  return xMatch;
}

/** Compute gradient of a function over the boundary */
double ComputeFunctionJetOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch, SMLVec3d *xOutGradient) 
{
  // Compute the image match
  double xMatch = IntegrateFunctionOverBoundary(xGrid, xAtoms, xWeights, fMatch);
  
  // Create a boundary point iterator and iterate over all points
  MedialBoundaryPointIterator *itBoundary = xGrid->NewBoundaryPointIterator();

  // Iterate through all boundary points
  while(!itBoundary->IsAtEnd())
    {
    // Evaluate the image gradient at this point
    unsigned int iBnd = itBoundary->GetIndex();
    SMLVec3d X = GetBoundaryPoint(itBoundary, xAtoms).X;
    fMatch->ComputeGradient(X, xOutGradient[iBnd]);
    xOutGradient[iBnd] *= xWeights[iBnd];

    // On to the next point
    ++(*itBoundary);
    }

  // Return match scaled by total weight
  return xMatch;
}

void TestTriangleAreaPartialDerivative()
{
  SMLVec3d X[3], Xp[3];
  for(size_t i=0;i<3;i++) for(size_t j=0;j<3;j++)
    {
    X[i][j] = rand() * 1.0 / RAND_MAX;
    Xp[i][j] = rand() * 1.0 / RAND_MAX;
    }

  double eps = 0.001;
  double A0 = TriangleArea(X[0], X[1], X[2]);
  double A1 = 
    TriangleArea(X[0] + eps * Xp[0], X[1] + eps * Xp[1], X[2] + eps * Xp[2]);
  double A2 = 
    TriangleArea(X[0] - eps * Xp[0], X[1] - eps * Xp[1], X[2] - eps * Xp[2]);

  cout << "Numerical derivative : " << 0.5 * (A1 - A2) / eps << endl;
  cout << "Analytical derivative : " << TriangleAreaPartialDerivative(
    X[0], X[1], X[2], Xp[0], Xp[1], Xp[2], A0) << endl;
}
  
void TestTetrahedronVolumePartialDerivative()
{
  double eps = 0.001;
  SMLVec3d X[4], Xp[4], X1[4], X2[4];
  for(size_t i=0;i<4;i++) for(size_t j=0;j<3;j++)
    {
    X[i][j] = rand() * 1.0 / RAND_MAX;
    Xp[i][j] = rand() * 1.0 / RAND_MAX;
    X1[i][j] = X[i][j] + eps * Xp[i][j];
    X2[i][j] = X[i][j] - eps * Xp[i][j];
    }

  double V0 = TetrahedronVolume( X[0], X[1], X[2], X[3] );
  double V1 = TetrahedronVolume( X1[0], X1[1], X1[2], X1[3] );
  double V2 = TetrahedronVolume( X2[0], X2[1], X2[2], X2[3] );

  cout << "Testing tetrahedron partial derivative computation" << endl;
  cout << "Numerical derivative : " << 0.5 * (V1 - V2) / eps << endl;
  cout << "Analytical derivative : " << 
    TetrahedronVolumePartialDerivative(
    X[0], X[1], X[2], X[3], 
    Xp[0], Xp[1], Xp[2], Xp[3]) << endl;
}
  



