#include "CoefficientMapping.h"
#include "PrincipalComponents.h"

PCACoefficientMapping
::PCACoefficientMapping(PrincipalComponents *pca, size_t nModes)
{ 
  this->pca = pca;
  this->n = nModes;
  this->m = pca->GetMean().size(); 
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::Apply(const Vec &C, const Vec &P)
{
  return C + pca->MapToFeatureSpaceZeroMean(P);
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) 
{ 
  return pca->MapToFeatureSpaceZeroMean(varP); 
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
{ 
  return varC; 
}

PCAPlusAffineCoefficientMapping
::PCAPlusAffineCoefficientMapping(
  GenericMedialModel *model, PrincipalComponents *pca, size_t nModes) :
  CompositionCoefficientMapping(
    new AffineTransformCoefficientMapping(model),
    new PCACoefficientMapping(pca, nModes))
{
}

PCAPlusAffineCoefficientMapping
::~PCAPlusAffineCoefficientMapping()
{
  delete this->f;
  delete this->g;
}
  
