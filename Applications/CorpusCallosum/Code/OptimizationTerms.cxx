#include "OptimizationTerms.h"
#include <cmath>

using namespace std;

typedef float PixelType;
typedef itk::Image< PixelType, Dimension >	ImageType;
typedef itk::ImageFileReader< ImageType >	ImageReader;
typedef itk::LinearInterpolateImageFunction< ImageType, double > ImageInterpolator;
typedef itk::DiscreteGaussianImageFilter< ImageType, ImageType > ImageSmoother;

CMRep2DOptimizationProblem::CMRep2DOptimizationProblem(const char *inputImageFile, double _imageBlurVariance, const int _dim, const Vector &TX) : imageBlurVariance(_imageBlurVariance), dim(_dim){
  /* Initialize image */

  ImageReader::Pointer reader = ImageReader::New();
  ImageSmoother::Pointer smooth = ImageSmoother::New();
  reader->SetFileName(inputImageFile);
  reader->Update();
  originalImage = ImageInterpolator::New();
  originalImage-> SetInputImage(reader->GetOutput());
  reader->SetFileName(inputImageFile);
  reader->Update();
  smooth->SetVariance(imageBlurVariance);
  smooth->SetInput(reader->GetOutput());
  smooth->Update();
  blurImage = ImageInterpolator::New();
  blurImage->SetInputImage(smooth->GetOutput()); 
  /*
  image = ImageInterpolator::New();
  image->SetInputImage(reader->GetOutput());
  */
  cm = new CMRep2D(dim);
  
  phi = new double[dim + 3];
  for (int i = 0; i < dim + 3; ++i){
    phi[i] = 1.0;
  };

  this->TX = TX;
  tSize = TX.size()/3;
  int tmp = tSize - 1;
  int counter = 0;
  for (; tmp > 1; counter ++) {
     tmp /= 2;
  }
  tJmax = counter - 1;

  // the coeffs remain unchanged in optimization, useful in multiscale optimzation
  fixedX.setSize(0);
  fixedSize = 0;
  fixedJmax = 0;

  // std::cout << "tJmax = " << tJmax << std::endl;
  stepEpsilon[0] = 1.0e-04;
  stepEpsilon[1] = 1.0e-04;
  stepEpsilon[2] = 1.0e-05;

  gradientScaleFactor[0] = 1.0;
  gradientScaleFactor[1] = 1.0;
  gradientScaleFactor[2] = 1.0e-07;
}

CMRep2DOptimizationProblem::~CMRep2DOptimizationProblem(){
  delete cm;
  delete[] phi;
}

void CMRep2DOptimizationProblem::setStepEpsilon(double xEpsilon, double yEpsilon, double rhoEpsilon) {
  stepEpsilon[0] = xEpsilon;
  stepEpsilon[1] = yEpsilon;
  stepEpsilon[2] = rhoEpsilon; 
}


void CMRep2DOptimizationProblem::setGradientScaleFactor(double xGf, double yGf, double rhoGf) {
  gradientScaleFactor[0] = xGf;
  gradientScaleFactor[1] = yGf;
  gradientScaleFactor[2] = rhoGf;
}

void CMRep2DOptimizationProblem::setFixedX(const Vector& _fixedX) {
  this->fixedX = _fixedX;
  fixedSize = fixedX.size()/3;
  int counter = 0;
  int tmp = fixedSize - 1;
  for(; tmp >1; counter++) {
    tmp /=2;
  }
  fixedJmax = counter - 1;
}

double CMRep2DOptimizationProblem::evaluate(const Vector &X) {
  int dSize = X.size()/3 + fixedSize;
  int tmp = dSize - 1;
  int counter = 0;
  for (; tmp > 1; counter ++) {
     tmp /= 2;
  }
  int dJmax = counter - 1;
  Vector allX;
  allX.setSize(dSize*3);
  for (int i = 0; i < fixedSize; i++) {
    allX(i) = fixedX(i);
    allX(i+dSize) = fixedX(i+fixedSize);
    allX(i+2*dSize) = fixedX(i+2*fixedSize);
  }
  for (int i = 0; i < dSize - fixedSize; i++) {
    allX(fixedSize + i) = X(i);
    allX(fixedSize + i + dSize) = X(i + dSize - fixedSize);
    allX(fixedSize + i + 2*dSize) = X(i + 2*dSize - 2*fixedSize);
  }
  //cout << "dSize = " << dSize << endl;
  //cout << "dJmax = " << dJmax << endl;
  TemplatedWaveletRep fx(TX.getDataArray(), tJmax, allX.getDataArray(), dJmax);
  TemplatedWaveletRep fy(TX.getDataArray() + tSize, tJmax, allX.getDataArray() + dSize, dJmax);
  TemplatedWaveletRep frho(TX.getDataArray() + tSize*2, tJmax, allX.getDataArray() + 2*dSize, dJmax);
  
  double  minPhi = cm->buildCMRep2D(fx, fy, frho, phi);
  double foldJac = cm->checkBoundaryFold();
  double arcLenVar = cm->getArcLenVar();
  arcLenVar *= 0.5;
  //  cout << " arcLenVar = " << arcLenVar << endl;
  double overlap = cm->computeAreaOverlap(blurImage, 0.8, 10, 5);
  //cout << "overlap = " << overlap << endl;
  double alpha = -50.0;
  double eps = 0.1;

  return  -overlap + arcLenVar + exp(alpha*(minPhi-eps)) + exp(alpha*(foldJac-eps));
  
}

double CMRep2DOptimizationProblem::computeOneJet(const Vector &X, Vector &XGrad) {
  double fx = evaluate(X);
  int dSize = X.size()/3;
  XGrad.setSize(X.size());

  Vector x1 = X;
   for (int i = 0; i < dSize; ++i) {
    x1(i) += stepEpsilon[0];
    double fx1 = evaluate(x1);
    x1(i) -= 2*stepEpsilon[0];
    double fx2 = evaluate(x1);
    x1(i) += stepEpsilon[0];
    XGrad(i) = gradientScaleFactor[0]*(fx1 - fx2) /(2*stepEpsilon[0]);

    //    cout << " fx1 fx2: " << fx1 << " , " << fx2 << endl;
   
    x1(i+dSize) += stepEpsilon[1];
    fx1 = evaluate(x1);
    x1(i+dSize) -= 2*stepEpsilon[1];
    fx2 = evaluate(x1);
    x1(i+dSize) += stepEpsilon[1];
    XGrad(i+dSize) = gradientScaleFactor[1]*(fx1 - fx2) /(2*stepEpsilon[1]);
   
    x1(i+2*dSize) += stepEpsilon[2];
    fx1 = evaluate(x1);
    x1(i+2*dSize) -= 2*stepEpsilon[2];
    fx2 = evaluate(x1);
    x1(i+2*dSize) += stepEpsilon[2];
    XGrad(i+2*dSize) = gradientScaleFactor[2]*(fx1 - fx2) / (2*stepEpsilon[2]);
    //     cout << " fx1 fx2: " << fx1 << " , " << fx2 << endl;
  }
  // check the norm of XGrad, make sure it's not too wild
  float len = XGrad.infinityNorm();
  if (len > 1.0e30) {
    XGrad /= len;
  }

  return fx;
}

bool CMRep2DOptimizationProblem::getBoundary(Vector &X, Vector &bx, Vector &by) {
  int dSize = X.size()/3 + fixedSize;
  int tmp = dSize - 1;
  int counter = 0;
  for (; tmp > 1; counter ++) {
     tmp /= 2;
  }
  int dJmax = counter - 1;
  Vector allX;
  allX.setSize(dSize*3);
  for (int i = 0; i < fixedSize; i++) {
    allX(i) = fixedX(i);
    allX(i+dSize) = fixedX(i+fixedSize);
    allX(i+2*dSize) = fixedX(i+2*fixedSize);
  }
  for (int i = 0; i < dSize - fixedSize; i++) {
    allX(fixedSize + i) = X(i);
    allX(fixedSize + i + dSize) = X(i + dSize - fixedSize);
    allX(fixedSize + i + 2*dSize) = X(i + 2*dSize - 2*fixedSize);
  }
  TemplatedWaveletRep fx(TX.getDataArray(), tJmax, allX.getDataArray(), dJmax);
  TemplatedWaveletRep fy(TX.getDataArray() + tSize, tJmax, allX.getDataArray() + dSize, dJmax);
  TemplatedWaveletRep frho(TX.getDataArray() + tSize*2, tJmax, allX.getDataArray() + 2*dSize, dJmax);
 
  bool nonnegative = cm->buildCMRep2D(fx, fy, frho, phi);
  if (!nonnegative) {
    return false;
  }
  cm->getBoundary(bx,by);
  return true;
}

void CMRep2DOptimizationProblem::computeScaleFactors(const Vector &X, Vector &scale) {
  double fx = evaluate(X);
  int dSize = X.size()/3;
  scale.setSize(3*dSize);

  Vector x1 = X;
  for (int i = 0; i < dSize; ++i) {
    x1(i) += stepEpsilon[0];
    double fx1 = evaluate(x1);
    x1(i) -= 2*stepEpsilon[0];
    double fx2 = evaluate(x1);
    x1(i) += stepEpsilon[0];
    scale(i) = (fx1 + fx2 - 2*fx) / (stepEpsilon[0]*stepEpsilon[0]);
    scale(i) = sqrt(fabs(scale(i)));
   
    x1(i+dSize) += stepEpsilon[1];
    fx1 = evaluate(x1);
    x1(i+dSize) -= 2*stepEpsilon[1];
    fx2 = evaluate(x1);
    x1(i+dSize) += stepEpsilon[1];
    scale(i+dSize) = (fx1 + fx2 - 2*fx) / (stepEpsilon[1]*stepEpsilon[1]);
    scale(i+dSize) = sqrt(fabs(scale(i+dSize)));   

    x1(i+2*dSize) += stepEpsilon[2];
    fx1 = evaluate(x1);
    x1(i+2*dSize) -= 2*stepEpsilon[2];
    fx2 = evaluate(x1);
    x1(i+2*dSize) += stepEpsilon[2];
    scale(i+2*dSize) = (fx1 + fx2 - 2*fx) / (stepEpsilon[2]*stepEpsilon[2]);
    scale(i+2*dSize) = sqrt(fabs(scale(i+2*dSize)));     
  }

  cout << "scale = " ;
  for (int i=0; i< 3*dSize; ++i) {
    cout << scale(i) << " , ";
  }
  cout << endl;

}

double CMRep2DOptimizationProblem::areaOverlapRatio(const Vector& X) {
  int dSize = X.size()/3 + fixedSize;
  int tmp = dSize - 1;
  int counter = 0;
  for (; tmp > 1; counter ++) {
     tmp /= 2;
  }
  int dJmax = counter - 1; 
  Vector allX;
  allX.setSize(dSize*3);
  for (int i = 0; i < fixedSize; i++) {
    allX(i) = fixedX(i);
    allX(i+dSize) = fixedX(i+fixedSize);
    allX(i+2*dSize) = fixedX(i+2*fixedSize);
  }
  for (int i = 0; i < dSize - fixedSize; i++) {
    allX(fixedSize + i) = X(i);
    allX(fixedSize + i + dSize) = X(i + dSize - fixedSize);
    allX(fixedSize + i + 2*dSize) = X(i + 2*dSize - 2*fixedSize);
  }
  TemplatedWaveletRep fx(TX.getDataArray(), tJmax, allX.getDataArray(), dJmax);
  TemplatedWaveletRep fy(TX.getDataArray() + tSize, tJmax, allX.getDataArray() + dSize, dJmax);
  TemplatedWaveletRep frho(TX.getDataArray() + tSize*2, tJmax, allX.getDataArray() + 2*dSize, dJmax);
  
  double  minPhi = cm->buildCMRep2D(fx, fy, frho, phi);
  if (minPhi <= 0) {
    cout << "not a valid CMRep!" << endl;
  }
  double foldJac = cm->checkBoundaryFold();
  if (foldJac <= 0) {
    cout << "boundary folds!" << endl;
  }
  double overlap = cm->computeAreaOverlap(originalImage, 0.8, 10, 5);
  double areaOfCMRep = cm->getAreaOfCMRep();

  double areaOfImage = 0.0;
  for (int i = 6; i < 250; i++) {
    for (int j = 6; j< 250; j++) {
     double index[2] = {i,j};
      areaOfImage += originalImage->Evaluate(index);
    }
  }
  cout << "area of CMRep = " << areaOfCMRep << endl;
  cout << "area of Image = " << areaOfImage << endl;
  cout << "area of overlap = " << overlap << endl;
  return overlap/(areaOfCMRep+areaOfImage-overlap);
}
