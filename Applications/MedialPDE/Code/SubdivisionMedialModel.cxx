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
::SetMesh(
  const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
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

  // Compute the coefficients and u/v arrays from the input data
  if(nCoeffSubs > 0)
    {
    // Initialize the arrays
    xCoefficients.set_size(mlCoefficient.nVertices * 4);
    uCoeff.set_size(mlCoefficient.nVertices);
    vCoeff.set_size(mlCoefficient.nVertices);

    // Apply the subdivision to the coefficients
    SubdivisionSurface::ApplySubdivision(
      C.data_block(), xCoefficients.data_block(), 4, mlCoefficient);

    // Apply to the u and v arrays
    SubdivisionSurface::ApplySubdivision(
      u.data_block(), uCoeff.data_block(), 1, mlCoefficient);
    SubdivisionSurface::ApplySubdivision(
      v.data_block(), vCoeff.data_block(), 1, mlCoefficient);
    }
  else
    {
    xCoefficients = C;
    uCoeff = u; vCoeff = v;
    }

  // Set the coefficient-level mesh as the new root (forget the input)
  mlCoefficient.SetAsRoot();

  // Subdivide the coefficient-level mesh up to the atom-level mesh
  cout << "Enter  SubdivisionSurface::RecursiveSubdivide" << endl;
  SubdivisionSurface::RecursiveSubdivide(
    &mlCoefficient, &mlAtom, nAtomSubs - nCoeffSubs);
  cout << "Got past SubdivisionSurface::RecursiveSubdivide" << endl;

  // Apply the subdivision to the u and v coordinates
  uAtom.set_size(mlAtom.nVertices); vAtom.set_size(mlAtom.nVertices);
  SubdivisionSurface::ApplySubdivision(
    uCoeff.data_block(), uAtom.data_block(), 1, mlAtom);
  SubdivisionSurface::ApplySubdivision(
    vCoeff.data_block(), vAtom.data_block(), 1, mlAtom);

  // Create the atoms array
  if(xAtoms) delete xAtoms;
  xAtoms = new MedialAtom[mlAtom.nVertices];

  // Copy the u, v values into the atoms
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    xAtoms[i].u = uAtom[i];
    xAtoms[i].v = vAtom[i];
    }

  // Set up the iteration context
  if(this->xIterationContext != NULL) delete this->xIterationContext;
  this->xIterationContext = new SubdivisionSurfaceMedialIterationContext(&mlAtom);

  // Set the coarse-to-fine descriptor
  xCTFDescriptor = SubdivisionSurfaceCoarseToFineMappingDescriptor(mlCoefficient.nVertices);
}

void
SubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  // Mesh type
  R["Grid.Type"] << "LoopSubdivision";

  // Save the subdivision level info
  R["Grid.Model.Atom.SubdivisionLevel"] << GetSubdivisionLevel();

}

void
SubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}

const CoarseToFineMappingDescriptor *
SubdivisionMedialModel
::GetCoarseToFineMappingDescriptor() const
{
  return &xCTFDescriptor;

}
