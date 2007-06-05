#include "BruteForceSubdivisionMedialModel.h"

BruteForceSubdivisionMedialModel::BruteForceSubdivisionMedialModel() :
  SubdivisionMedialModel()
{
  xAtoms = NULL;
}

BruteForceSubdivisionMedialModel::~BruteForceSubdivisionMedialModel()
{
}

void
BruteForceSubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Call the parent method
  SubdivisionMedialModel::SetMesh(mesh, C, u, v, nAtomSubs, nCoeffSubs);

  // Pass the mesh to the loop scheme
  xLoopScheme.SetMeshLevel(&mlAtom);
}

void
BruteForceSubdivisionMedialModel
::ComputeAtoms(const double *xHint)
{
  size_t i;

  // The first step is to compute the X and rho of the atoms based on the
  // coefficients. This step is performed in this class, not in the solver
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i]; a.X.fill(0.0); a.R = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; double w = it.Value();
      a.X += w * xCoefficients.extract(3, j);
      a.R += w * xCoefficients[j+3];
      }

    // Set F = R^2
    a.F = a.R * a.R;
    }

  // We must now compute the derivatives of X and R to make real atoms
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];
    
    // Compute partial derivatives
    a.Fu = xLoopScheme.Fu(i, xAtoms);
    a.Fv = xLoopScheme.Fv(i, xAtoms);
    a.Xu = xLoopScheme.Xu(i, xAtoms);
    a.Xv = xLoopScheme.Xv(i, xAtoms);

    // Compute the differential geometry
    a.G.SetOneJet(a.X.data_block(), a.Xu.data_block(), a.Xv.data_block());

    // Compute the normal vector
    a.ComputeNormalVector();

    // Compute all the components of the atoms
    a.ComputeBoundaryAtoms(!mlAtom.IsVertexInternal(i));

    // Compute some form of strain tensor. Loop over all neighbors
    /*
    SMLVec3d XU, XV; XU.fill(0); XV.fill(0);
    EdgeWalkAroundVertex it(&mlAtom, i);

    // Flag all neighbors
    for(; !it.IsAtEnd(); ++it)
      {
      MedialAtom &aj = xAtoms[it.MovingVertexId()];
      XU += (aj.X - a.X) / (aj.u - a.u);
      XV += (aj.X - a.X) / (aj.v - a.v);
      }

    double E = dot_product(XU, XU);
    double F = dot_product(XU, XV);
    double G = dot_product(XV, XV);
    a.G.xChristoffelSecond[0][0][0] = E + G;
    */
    

    // Compute the Riemannian gradient of the parameterization (for stretch tensors)
    /*
    double udu = xLoopScheme.Uu(i, xAtoms);
    double udv = xLoopScheme.Uv(i, xAtoms);
    double vdu = xLoopScheme.Vu(i, xAtoms);
    double vdv = xLoopScheme.Vv(i, xAtoms);
    double xMagGradU = 
      udu * udu * a.G.xContravariantTensor[0][0] + 
      2 * udu * udv * a.G.xContravariantTensor[0][1] + 
      udv * udv * a.G.xContravariantTensor[1][1];
    double xMagGradV = 
      vdu * vdu * a.G.xContravariantTensor[0][0] + 
      2 * vdu * vdv * a.G.xContravariantTensor[0][1] + 
      vdv * vdv * a.G.xContravariantTensor[1][1];
    double xStretch = xMagGradV + xMagGradU;
    a.G.xChristoffelSecond[0][0][0] = xStretch;
    */
    }

  // On the third pass, we compute the second order partial derivatives
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // The second order derivatives are set to zero
    a.Xuu.fill(0.0); a.Xuv.fill(0.0); a.Xvv.fill(0.0);
    
    // Compute the derivatives of the normal vector
    SMLVec3d Nu = xLoopScheme.Nu(i, xAtoms);
    SMLVec3d Nv = xLoopScheme.Nv(i, xAtoms);

    // Compute the principal curvatures
    double ee = - dot_product(Nv, a.Xv);
    double ff = -0.5 * (dot_product(Nu, a.Xv) + dot_product(Nv, a.Xu));
    double gg = - dot_product(Nu, a.Xu);
    a.xMeanCurv = 0.5 * (
      ee * a.G.xContravariantTensor[0][0] +
      2 * ff * a.G.xContravariantTensor[1][0] + 
      gg * a.G.xContravariantTensor[1][1]);
    a.xGaussCurv = (ee * gg - ff * ff) * a.G.gInv;
    a.Xuu = Nu;
    a.Xvv = Nv;
    }
}

BruteForceSubdivisionMedialModel::Vec
BruteForceSubdivisionMedialModel::GetHintArray() const
{
  Vec xHint(1, 0.0);
  return xHint;
}

void
BruteForceSubdivisionMedialModel
::PrepareAtomsForVariationalDerivative(
  const Vec &xVariation, MedialAtom *dAtoms) const
{
  size_t i;

  // First, we must label the dependency structure. To begin with, for each atom
  // we must mark whether is is affected by the variation or not. At the same time
  // we can compute the derivative of X and R with respect to the variation 
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &da = dAtoms[i]; da.X.fill(0.0); da.R = 0.0;

    // Initialize the user_data flag to indicate no dependence on the variation
    da.user_data = 0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::ConstRowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; 
      double w = it.Value();

      SMLVec3d dXdVar = w * xVariation.extract(3, j);
      double dRdVar = w * xVariation[j+3];

      da.X += dXdVar;
      da.R += dRdVar;

      // Check if there is a dependency
      if(dXdVar.squared_magnitude() > 0.0 || dRdVar != 0.0)
        da.user_data = 1;
      }
    }

  // Propagate the dependency flag (user_data) to the 2-nd ring of neighbors. This
  // is because the neighbors are used in the computation of Xu and Xv
  for(size_t level = 0; level < 2; level++)
    {
    for(i = 0; i < mlAtom.nVertices; i++)
      {
      if(dAtoms[i].user_data == 1)
        {
        // Create a walk around the vertex
        EdgeWalkAroundVertex it(&mlAtom, i);

        // Flag all neighbors
        for(; !it.IsAtEnd(); ++it)
          dAtoms[it.MovingVertexId()].user_data = 1;
        }
      }
    }

  // Next, precompute the first partial derivatives of X and R
  for(i = 0; i < mlAtom.nVertices; i++) 
    {
    MedialAtom &da = dAtoms[i];

    // Set the values of Xu and Xv
    da.Xu = xLoopScheme.Xu(i, dAtoms);
    da.Xv = xLoopScheme.Xv(i, dAtoms);

    // Also compute Ru and Rv (although atoms don't have a place to stick this)
    da.Ru = xLoopScheme.Ru(i, dAtoms);
    da.Rv = xLoopScheme.Rv(i, dAtoms);


    // The second order derivatives are set to zero
    da.Xuu.fill(0.0); da.Xuv.fill(0.0); da.Xvv.fill(0.0);
    }

  // Set the derivative atoms unaffected by the computation to zeros
  for(i = 0; i < mlAtom.nVertices; i++) 
    {
    if(dAtoms[i].user_data == 0)
      dAtoms[i].SetAllDerivativeTermsToZero();
    }
}

void
BruteForceSubdivisionMedialModel
::ComputeAtomGradient(std::vector<MedialAtom *> &xVariations)
{
  size_t i, q;

  // Precompute common terms for atom derivatives
  MedialAtom::DerivativeTerms *dt = 
    new MedialAtom::DerivativeTerms[mlAtom.nVertices];
  for(i = 0; i < mlAtom.nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(dt[i]);

  // Compute the derivative for each variation
  for(q = 0; q < xVariations.size(); q++)
    {
    // Get the derivative atoms for this variation
    MedialAtom *dAtoms = xVariations[q];

    // Compute F at each atom
    for(i = 0; i < mlAtom.nVertices; i++) 
      if(dAtoms[i].user_data)
        dAtoms[i].F = 2.0 * dAtoms[i].R * xAtoms[i].R;

    // Compute dF/dv for each atom
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      // Only consider atoms affected by the variation
      if(dAtoms[i].user_data)
        {
        // Get the current atom and the derivative (which we are computing)
        MedialAtom &a = xAtoms[i];
        MedialAtom &da = dAtoms[i];

        // Set Fu and Fv
        da.Fu = xLoopScheme.Fu(i, dAtoms);
        da.Fv = xLoopScheme.Fv(i, dAtoms);
      
        // Compute the metric tensor derivatives of the atom
        xAtoms[i].ComputeMetricTensorDerivatives(da);

        // Compute the derivatives of the boundary nodes
        xAtoms[i].ComputeBoundaryAtomDerivatives(da, dt[i]);
        }
      }

    // Comput the curvatures
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      // Only consider atoms affected by the variation
      if(dAtoms[i].user_data)
        {
        // Get the current atom and the derivative (which we are computing)
        MedialAtom &a = xAtoms[i];
        MedialAtom &da = dAtoms[i];

        // Compute the derivative of mean and gauss curvatures
        SMLVec3d Nu = xLoopScheme.Nu(i, xAtoms);
        SMLVec3d Nv = xLoopScheme.Nv(i, xAtoms);
        SMLVec3d dNu = xLoopScheme.Nu(i, dAtoms);
        SMLVec3d dNv = xLoopScheme.Nv(i, dAtoms);

        // TODO: this should be precomputed, figure out a memory-sane solution
        double ee = - dot_product(Nv, a.Xv);
        double ff = -0.5 * (dot_product(Nu, a.Xv) + dot_product(Nv, a.Xu));
        double gg = - dot_product(Nu, a.Xu);

        // Compute the derivatives of 2nd fundamental form
        double dee = - (dot_product(dNv, a.Xv) + dot_product(Nv, da.Xv));
        double dff = -0.5 * 
          ( dot_product(dNv, a.Xu) + dot_product(Nv, da.Xu) + 
            dot_product(dNu, a.Xv) + dot_product(Nu, da.Xv) );
        double dgg = - (dot_product(dNu, a.Xu) + dot_product(Nu, da.Xu));

        // Compute the derivatives of the mean curvature
        da.xMeanCurv = 0.5 * (
          dee * a.G.xContravariantTensor[0][0] + 
          ee * da.G.xContravariantTensor[0][0] +
          2 * (
            dff * a.G.xContravariantTensor[1][0] +
            ff * da.G.xContravariantTensor[1][0] ) +
          dgg * a.G.xContravariantTensor[1][1] + 
          gg * da.G.xContravariantTensor[1][1]);

        // Compute the derivatives of the Gauss curvature
        da.xGaussCurv = 
          ( (dee * gg + ee * dgg - 2 * ff * dff) * a.G.g - 
            (ee * gg - ff * ff) * da.G.g ) / (a.G.g * a.G.g);

        da.Xuu = dNu;
        da.Xvv = dNv;
        }
      }
    }

  delete dt;
}


void
BruteForceSubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  SubdivisionMedialModel::WriteToRegistry(R);

  // Set the model subtype
  R["Grid.Model.SolverType"] << "BruteForce";
}


void
BruteForceSubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}

