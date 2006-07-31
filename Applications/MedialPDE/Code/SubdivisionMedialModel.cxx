#include "SubdivisionMedialModel.h"
#include "SubdivisionSurfaceMedialIterationContext.h"

SubdivisionMedialModel::SubdivisionMedialModel() :
  xCTFDescriptor(0)
{
  xAtoms = NULL;
  xIterationContext = NULL;
  xSubdivisionLevel = 0;
}

void
SubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, const vnl_vector<double> &C,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Validity check
  assert(nAtomSubs >= nCoeffSubs);

  // Set the subdivision level
  xSubdivisionLevel = nAtomSubs - nCoeffSubs;

  // Create a vector of mesh levels
  vector<const MeshLevel *> xTempLevels; xTempLevels.push_back(&mesh);

  // Subdivide the input mesh into the coefficient-level mesh
  SubdivisionSurface::RecursiveSubdivide(&mesh, &mlCoefficient, nCoeffSubs);

  // Compute the coefficients from the input data
  if(nCoeffSubs > 0)
    {
    xCoefficients.set_size(mlCoefficient.nVertices * 4);
    SubdivisionSurface::ApplySubdivision(
      C.data_block(), xCoefficients.data_block(), 4, mlCoefficient);
    }
  else
    xCoefficients = C;
    
  // Set the coefficient-level mesh as the new root (forget the input)
  mlCoefficient.SetAsRoot();

  // Subdivide the coefficient-level mesh up to the atom-level mesh
  SubdivisionSurface::RecursiveSubdivide(&mlCoefficient, &mlAtom, nAtomSubs - nCoeffSubs);

  // Set up the iteration context
  if(this->xIterationContext != NULL) delete this->xIterationContext;
  this->xIterationContext = new SubdivisionSurfaceMedialIterationContext(&mlAtom);

  // Pass the mesh information to the solver
  xSolver.SetMeshTopology(&mlAtom);

  // The the atoms to point into the solver
  xAtoms = xSolver.GetAtomArray();

  // Set the coarse-to-fine descriptor
  xCTFDescriptor = SubdivisionSurfaceCoarseToFineMappingDescriptor(mlCoefficient.nVertices);
}

void
SubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, SMLVec3d *X, double *rho, 
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Flatten the input data into a coefficient array
  Vec C0(mesh.nVertices * 4, 0.0);
  for(size_t i = 0, j = 0; i < mesh.nVertices; i++)
    {
    C0[j++] = X[i][0]; C0[j++] = X[i][1]; C0[j++] = X[i][2]; 
    C0[j++] = rho[i];
    }

  // Pass in to the other SetMesh method
  this->SetMesh(mesh, C0, nAtomSubs, nCoeffSubs);
}

void 
SubdivisionMedialModel
::ComputeAtoms(MedialAtom *xInitialSolution)
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
  if(xInitialSolution != NULL)
    for(i = 0; i < mlAtom.nVertices; i++)
      xAtoms[i].F = xInitialSolution[i].F;

  // Now have the solver solve the equation
  xSolver.SolveEquation(NULL, true);
  xSolver.SolveEquation(NULL, true);
}

void
SubdivisionMedialModel
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
SubdivisionMedialModel
::ComputeAtomGradient(std::vector<MedialAtom *> &dAtoms)
{
  xSolver.ComputeGradient(dAtoms);
}

const CoarseToFineMappingDescriptor *
SubdivisionMedialModel
::GetCoarseToFineMappingDescriptor() const
{

}
