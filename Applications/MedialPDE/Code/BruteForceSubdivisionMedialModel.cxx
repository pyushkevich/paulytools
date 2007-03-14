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
    MedialAtom &a = xAtoms[i]; a.X.fill(0.0); a.xLapR = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; double w = it.Value();
      a.X += w * xCoefficients.extract(3, j);
      a.R += w * xCoefficients[j+3];
      }
    }

  // We must now compute the derivatives of X and R to make real atoms
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // Compute partial derivatives
    a.Fu = xLoopScheme.PartialPhi(0, i, xAtoms);
    a.Fv = xLoopScheme.PartialPhi(1, i, xAtoms);
    a.Xu = xLoopScheme.TangentX(0, i, xAtoms);
    a.Xv = xLoopScheme.TangentX(1, i, xAtoms);

    // Compute the differential geometry
    a.G.SetOneJet(a.X.data_block(), a.Xu.data_block(), a.Xv.data_block());

    // Compute the normal vector
    a.ComputeNormalVector();

    // Compute all the components of the atoms
    a.ComputeBoundaryAtoms(!mlAtom.IsVertexInternal(i));
    }

  // Finally we 'fix' the boundary atoms to have low grad R


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
  // the values of xLapR and X.
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &da = dAtoms[i]; da.X.fill(0.0); da.xLapR = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::ConstRowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; double w = it.Value();
      da.X += w * xVariation.extract(3, j);
      da.xLapR += w * xVariation[j+3];
      }
    }
}

void
BruteForceSubdivisionMedialModel
::ComputeAtomGradient(std::vector<MedialAtom *> &dAtoms)
{
  // TODO: Write this method!
}
