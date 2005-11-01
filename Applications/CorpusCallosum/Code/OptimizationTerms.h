#ifndef _OptimizationTerms_H
#define _OptimizationTerms_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "CMRep2D.h"
#include "TemplatedWaveletRep.h"

#include "optima.h"

class CMRep2DOptimizationProblem : public DifferentiableFunction {
 private:
  double imageBlurVariance;
  ImageInterpolator::Pointer originalImage;
  ImageInterpolator::Pointer blurImage;
  const int dim;
  CMRep2D *cm;
  double *phi;
  Vector TX;
  int tSize;
  int tJmax;
  Vector fixedX;
  int fixedSize;
  int fixedJmax;
  double stepEpsilon[3];
  double gradientScaleFactor[3];
 
 public:
  /* Initialize the problem */
  CMRep2DOptimizationProblem(const char *inputImageFile, double _imageBlurVariance, const int _dim, const Vector& TX);
  ~CMRep2DOptimizationProblem();

  void setStepEpsilon(double xEpsilon, double yEpsilon, double rhoEpsilon);
  void setGradientScaleFactor(double xGf, double yGf, double rhoGf);
  void setFixedX(const Vector& _fixedX);
  double evaluate(const Vector& X);
  double computeOneJet(const Vector& X, Vector& XGrad);
  bool CMRep2DOptimizationProblem::getBoundary(Vector& X, Vector& bx, Vector& by);
  void computeScaleFactors(const Vector& X, Vector& scale);
  double areaOverlapRatio(const Vector& X); 
};

#endif
