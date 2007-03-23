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
    }

  // On the third pass, we compute the second order partial derivatives
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];
    
    // Compute partial derivatives
    a.Xuu = xLoopScheme.Xuu(i, xAtoms);
    a.Xuv = xLoopScheme.Xuv(i, xAtoms);
    a.Xvv = xLoopScheme.Xvv(i, xAtoms);

    // Compute the differential geometry
    a.ComputeDifferentialGeometry();
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
  // This method must compute the terms in the derivative atoms that will not
  // change as the coefficients themselves change. This simply means setting
  // the values of R and X.
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &da = dAtoms[i]; da.X.fill(0.0); da.R = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::ConstRowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; double w = it.Value();
      da.X += w * xVariation.extract(3, j);
      da.R += w * xVariation[j+3];
      }
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

    // Compute dF/dv for each atom
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      // Get the current atom and the derivative (which we are computing)
      MedialAtom &a = xAtoms[i];
      MedialAtom &da = dAtoms[i];

      // Set dF/dv = d(R^2)/dv = 2R dR/dV
      da.F = 2.0 * da.R * a.R;
      }

    // Compute the derivatives of Fu and Fv
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      MedialAtom &da = dAtoms[i];
      
      // Set the values of Fu and Fv
      da.Fu = xLoopScheme.Fu(i, dAtoms);
      da.Fv = xLoopScheme.Fv(i, dAtoms);
      da.Xu = xLoopScheme.Xu(i, dAtoms);
      da.Xv = xLoopScheme.Xv(i, dAtoms);
      }

    // Compute second derivatives
    for(i = 0; i < mlAtom.nVertices; i++) 
      {
      MedialAtom &da = dAtoms[i];

      da.Xuu = xLoopScheme.Xuu(i, dAtoms);
      da.Xuv = xLoopScheme.Xuv(i, dAtoms);
      da.Xvv = xLoopScheme.Xvv(i, dAtoms);

      // Compute the metric tensor derivatives of the atom
      xAtoms[i].ComputeMetricTensorDerivatives(da);

      // Compute the derivatives of the boundary nodes
      xAtoms[i].ComputeBoundaryAtomDerivatives(da, dt[i]);
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

