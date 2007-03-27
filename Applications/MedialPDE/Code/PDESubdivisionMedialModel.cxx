#include "PDESubdivisionMedialModel.h"
#include "SubdivisionSurfaceMedialIterationContext.h"

/**
 * Constructor does nothing, just calls parent constructor
 */
PDESubdivisionMedialModel::PDESubdivisionMedialModel() :
  SubdivisionMedialModel()
{
}

void
PDESubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Call the parent method
  SubdivisionMedialModel::SetMesh(mesh, C, u, v, nAtomSubs, nCoeffSubs);

  // Pass the mesh information to the solver
  xSolver.SetMeshTopology(&mlAtom, xAtoms);
}

void 
PDESubdivisionMedialModel
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
      a.xLapR += w * xCoefficients[j+3];
      }
    }

  // If the initial solution is specified, take the F values from it
  if(xHint != NULL) 
    {
    for(i = 0; i < mlAtom.nVertices; i++)
      xAtoms[i].F = xHint[i];
    }

  // Now have the solver solve the equation
  xSolver.SolveEquation(NULL, true);
}

PDESubdivisionMedialModel::Vec
PDESubdivisionMedialModel::GetHintArray() const
{
  Vec xHint(mlAtom.nVertices, 0.0);
  for(size_t i = 0; i < xHint.size(); i++)
    xHint[i] = xAtoms[i].F;
  return xHint;
}

void
PDESubdivisionMedialModel
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
PDESubdivisionMedialModel
::ComputeAtomGradient(std::vector<MedialAtom *> &dAtoms)
{
  xSolver.ComputeGradient(dAtoms);
}


void
PDESubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  SubdivisionMedialModel::WriteToRegistry(R);

  // Set the model subtype
  R["Grid.Model.SolverType"] << "PDE";

  // Save the values of R computed at the atom level
  Vec phi = GetPhi();
  R["Grid.PhiAvailable"] << true;
  R.Folder("Grid.Phi").PutArray(phi.size(), phi.data_block());
}

void
PDESubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}
