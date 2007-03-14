#ifndef __PDESubdivisionMedialModel_h_
#define __PDESubdivisionMedialModel_h_

#include "SubdivisionMedialModel.h"

class PDESubdivisionMedialModel : public SubdivisionMedialModel
{
public:
  // Mesh level definition from parent
  typedef SubdivisionMedialModel::MeshLevel MeshLevel;

  PDESubdivisionMedialModel();

  /** Set the mesh and initial coefficient array */
  void SetMesh(
    const MeshLevel &mesh, 
    const Vec &C, const Vec &u, const Vec &v,
    size_t nAtomSubs, size_t nCoeffSubs);
  
  /** Get the solver */
  MeshMedialPDESolver *GetSolver() 
    { return &xSolver; }

  /** Get the hint array (solution phi) */
  virtual Vec GetHintArray() const;

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
};


#endif // __PDESubdivisionMedialModel_h_
