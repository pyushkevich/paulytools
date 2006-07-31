#ifndef __SubdivisionMedialModel_h_
#define __SubdivisionMedialModel_h_

#include "GenericMedialModel.h"
#include "MeshMedialPDESolver.h"
#include "SubdivisionSurface.h"

class SubdivisionMedialModel : public GenericMedialModel
{
public:
  typedef SubdivisionSurface::MeshLevel MeshLevel;

  SubdivisionMedialModel();
  
  /** 
   * This method associates the medial model with a mesh. The two parameters
   * that follow are levels of subdivision: the number of levels at which the
   * atoms are interpolated and a smaller number, at which the coefficients
   * are generated. By default, the nodes in the mesh itself serve as the
   * coefficients, but it is possible to request that the input mesh first be
   * subdivided and the nodes of that mesh be used as coefficients
   */
  void SetMesh(const MeshLevel &mesh, SMLVec3d *X, double *rho, 
    size_t nAtomSubs, size_t nCoeffSubs = 0);

  /**
   * A version of SetMesh that takes a coefficient array rather than X and rho
   */
  void SetMesh(const MeshLevel &mesh, const vnl_vector<double> &C,
    size_t nAtomSubs, size_t nCoeffSubs);

  /** Get the solver */
  MeshMedialPDESolver *GetSolver() 
    { return &xSolver; }

  /** Get the number of optimizable coefficients that define this model */
  size_t GetNumberOfCoefficients() const
    { return xCoefficients.size(); }

  /** Get the array of (optimizable) coefficients that define this model */
  const Vec GetCoefficientArray() const
    { return xCoefficients; }

  /** Get one of the coefficients defining this model */
  double GetCoefficient(size_t i) const
    { return xCoefficients[i]; }

  /** Set one of the coefficients defining this model */
  void SetCoefficient(size_t i, double x)
    { xCoefficients[i] = x; }

  /** Set the array of coefficients defining this model */
  void SetCoefficientArray(const Vec &xData)
    { xCoefficients = xData; }

  /** Compute the atoms from given a set of coefficients.  */
  void ComputeAtoms(MedialAtom *xInitialSolution = NULL);

  /**
   * This method is called before multiple calls to the ComputeJet routine. Given a  
   * of variation (direction in which directional derivatives will be computed), and 
   * an array of atoms corresponding to the variation, this method will precompute the
   * parts of the atoms that do not depend on the position in coefficient space where 
   * the gradient is taken. This method must be called for each variation in the Jet.
   */
  void PrepareAtomsForVariationalDerivative(
    const Vec &xVariation, MedialAtom *dAtoms) const;

  /**
   * This method computes the 'gradient' (actually, any set of directional derivatives)
   * for the current values of the coefficients. As in ComputeAtoms, a set of initial 
   * values can be passed in. By default, the previous solution is used to initialize.
   */
  void ComputeAtomGradient(std::vector<MedialAtom *> &dAtoms);

  /**
   * Load a medial model from the Registry.
   */
  void ReadFromRegistry(Registry &folder);

  /**
   * Save the model to registry folder
   */
  void WriteToRegistry(Registry &folder);

  /**
   * Get the subdivision level
   */
  size_t GetSubdivisionLevel() const { return xSubdivisionLevel; }

  /**
   * Get a pointer to the affine transform descriptor corresponding to this class. 
   * The descriptor is a lightweight object, and its allocation should be
   * managed by the child of GenericMedialModel.
   */
  const AffineTransformDescriptor *GetAffineTransformDescriptor() const
    { return &xAffineDescriptor; }

  /**
   * Get a pointer to the coarse-to-fine masking descriptor corresponding to
   * this class. This descriptor is managed internally by the child of
   * GenericMedialModel.
   */
  const CoarseToFineMappingDescriptor *
    GetCoarseToFineMappingDescriptor() const;

  /** Get a pointer to the coefficient mesh stored in this model */
  const MeshLevel &GetCoefficientMesh() const
    { return mlCoefficient; }

  /** Get a pointer to the atom mesh stored in this model */
  const MeshLevel &GetAtomMesh() const
    { return mlAtom; }

  /** Get the vector of phi-values for this model (equal to xAtoms.F) */
  Vec GetPhi() const 
    { 
    Vec phi(this->GetNumberOfAtoms());
    for(size_t i = 0; i < phi.size(); i++)
      phi[i] = xAtoms[i].F;
    return phi;
    }

  /** Set the vector of phi-values for the model */
  void SetPhi(const Vec &phi)
    {
    for(size_t i = 0; i < phi.size(); i++)
      xAtoms[i].F = phi[i];
    }

private:
  // Vector typedef
  typedef vnl_vector<double> Vec;

  // PDE solver which does most of the work in this method
  MeshMedialPDESolver xSolver;

  // Mesh levels corresponding to the coefficients and the atoms
  MeshLevel mlCoefficient, mlAtom;

  // The array of coefficients (4-tuples at each vertex)
  Vec xCoefficients;

  // Coarse-to-fine mapping descriptor (dummy)
  SubdivisionSurfaceCoarseToFineMappingDescriptor xCTFDescriptor;
  
  // Affine transform mapping descriptor
  PointArrayAffineTransformDescriptor xAffineDescriptor;

  // The subdivision level (number of subdivisions from mlCoefficient to
  // mlAtom)
  size_t xSubdivisionLevel;

  friend class SubdivisionMedialModelIO;
};


#endif // __SubdivisionMedialModel_h_
