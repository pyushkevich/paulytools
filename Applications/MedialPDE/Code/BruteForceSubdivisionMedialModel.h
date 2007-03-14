#ifndef __BruteForceSubdivisionMedialModel_h_
#define __BruteForceSubdivisionMedialModel_h_

#include "SubdivisionMedialModel.h"
#include "MedialAtomGrid.h"

/**
 * This is a type of subdivision medial model that does not solve any PDE.
 * Instead, in this model the subdivision surface represents the function
 * (xyzr) and no equality contraints are enforced. Who knows if that is going to
 * work, but at least we can try!
 */
class BruteForceSubdivisionMedialModel : public SubdivisionMedialModel
{
public:
  // Mesh level definition from parent
  typedef SubdivisionMedialModel::MeshLevel MeshLevel;

  BruteForceSubdivisionMedialModel();
  ~BruteForceSubdivisionMedialModel();

  void SetMesh(const MeshLevel &mesh, 
    const Vec &C, const Vec &u, const Vec &v,
    size_t nAtomSubs, size_t nCoeffSubs);

  /** Get the hint array. This returns a single double, it's a dummy method */
  Vec GetHintArray() const;

  /** Compute the atoms from given a set of coefficients.  */
  void ComputeAtoms(const double *xHint = NULL);

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

private:
  // Vector typedef
  typedef vnl_vector<double> Vec;

  // Loop scheme for computing tangents
  MedialAtomLoopScheme xLoopScheme;
};


#endif // __BruteForceSubdivisionMedialModel_h_
