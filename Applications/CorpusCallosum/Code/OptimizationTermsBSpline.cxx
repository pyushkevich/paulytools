#include "OptimizationTermsBSpline.h"
#include <cmath>

using namespace std;

typedef float PixelType;
typedef itk::Image< PixelType, Dimension >	ImageType;
typedef itk::ImageFileReader< ImageType >	ImageReader;
typedef itk::LinearInterpolateImageFunction< ImageType, double > ImageInterpolator;
typedef itk::DiscreteGaussianImageFilter< ImageType, ImageType > ImageSmoother;

CMRep2DOptimizationProblem::CMRep2DOptimizationProblem(char *_inputImageFile, double _imageBlurVariance, const int _dim) : inputImageFile(_inputImageFile), imageBlurVariance(_imageBlurVariance), dim(_dim){
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
  fixedX.setSize(3);
  for (int i = 0; i<3; i++) {
    fixedX(i) = 0.0;
  };
  phi = new double[dim + 3];
  phi[0] = -100.0;
  phi[dim+2] = -100.0;
  for (int i = 1; i < dim + 2; ++i){
    phi[i] = 1.0;
  };

  stepEpsilon[0] = 1.0e-04;
  stepEpsilon[1] = 1.0e-04;
  stepEpsilon[2] = 1.0e-05;

  gradientScaleFactor[0] = 1.0;
  gradientScaleFactor[1] = 1.0;
  gradientScaleFactor[2] = 1.0;
}

CMRep2DOptimizationProblem::~CMRep2DOptimizationProblem(){
  delete cm;
  delete[] phi;
}
void CMRep2DOptimizationProblem::setBlurVariance(double _imageBlurVariance) {
  imageBlurVariance = _imageBlurVariance;
  ImageReader::Pointer reader = ImageReader::New();
  ImageSmoother::Pointer smooth = ImageSmoother::New();
  reader->SetFileName(inputImageFile);
  reader->Update();
  if(imageBlurVariance > 0.0){
    smooth->SetVariance(_imageBlurVariance);
    smooth->SetInput(reader->GetOutput());
    smooth->Update();
    blurImage->SetInputImage(smooth->GetOutput()) ;
  }
  if(imageBlurVariance==0.0){
    blurImage->SetInputImage(reader->GetOutput());
  }
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

void  CMRep2DOptimizationProblem::setFixedX(const Vector& _fixedX) {
 this->fixedX = _fixedX;
}

double CMRep2DOptimizationProblem::evaluate(const Vector &X) {
  int fixedSize = fixedX.size()/3;
  int xSize = X.size()/3;

  SumBSplineCubUni fx(fixedX.getDataArray(), fixedSize, X.getDataArray(), xSize);
  SumBSplineCubUni fy(fixedX.getDataArray() + fixedSize, fixedSize, X.getDataArray() + xSize, xSize);
  SumBSplineCubUni frho(fixedX.getDataArray() + 2*fixedSize, fixedSize, X.getDataArray() + 2*xSize, xSize);

  double  minPhi = cm->buildCMRep2D(fx, fy, frho, phi);
  double foldJac = cm->checkBoundaryFold();
  double arcLenVar = cm->getArcLenVar();
  arcLenVar *= 1;
  // cout << " arcLenVar = " << arcLenVar << endl;
  double overlap = cm->computeAreaOverlap(blurImage, 0.8, 10, 5);
  // double overlap = cm->computeAreaOverlap(originalImage, 0.8, 10, 5);
  //  cout << "overlap = " << overlap << endl;
  //  cout << "minPhi = " << minPhi << endl;
  // cout << "foldJac = " << foldJac << endl;
  double alpha = -50.0;
  double eps1 = 0.03;
  double eps2 = 0.03;
  //   if( minPhi < 0) {
  //  return -minPhi;  
  // }
  // if( foldJac < 0) {
  // return -foldJac; 
  // }

  return  -overlap + arcLenVar + exp(alpha*(minPhi-eps1)) + exp(alpha*(foldJac-eps2));
  
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
    cout << " Gradient overflow!" << endl;
    XGrad /= len;
  }

  return fx;
}

bool CMRep2DOptimizationProblem::getBoundary(Vector &X, Vector &bx, Vector &by) {
  int fixedSize = fixedX.size()/3;
  int xSize = X.size()/3;

  SumBSplineCubUni fx(fixedX.getDataArray(), fixedSize, X.getDataArray(), xSize);
  SumBSplineCubUni fy(fixedX.getDataArray() + fixedSize, fixedSize, X.getDataArray() + xSize, xSize);
  SumBSplineCubUni frho(fixedX.getDataArray() + 2*fixedSize, fixedSize, X.getDataArray() + 2*xSize, xSize);
 
  double minPhi = cm->buildCMRep2D(fx, fy, frho, phi);
  if (minPhi <= 0) {
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
  int fixedSize = fixedX.size()/3;
  int xSize = X.size()/3;

  SumBSplineCubUni fx(fixedX.getDataArray(), fixedSize, X.getDataArray(), xSize);
  SumBSplineCubUni fy(fixedX.getDataArray() + fixedSize, fixedSize, X.getDataArray() + xSize, xSize);
  SumBSplineCubUni frho(fixedX.getDataArray() + 2*fixedSize, fixedSize, X.getDataArray() + 2*xSize, xSize);
  
  double  minPhi = cm->buildCMRep2D(fx, fy, frho, phi);
  if (minPhi <= 0) {
    cout << "not a valid CMRep!" << endl;
  }
  double foldJac = cm->checkBoundaryFold();
  if (foldJac <= 0) {
    cout << "boundary folds!" << endl;
  }
  double overlapOri = cm->computeAreaOverlap(originalImage, 0.8, 200, 500);
  double overlap = cm->overlapAndCMRep(originalImage, 0.8, 200, 500, 0.0);
  double areaOfCMRep = cm->getAreaOfCMRep();

  double areaOfImage = 0.0;
  for (int i = 6; i < 250; i++) {
    for (int j = 6; j< 250; j++) {
     double index[2] = {i,j};
     areaOfImage +=  originalImage->Evaluate(index);
    }
  }
  cout << "area of CMRep = " << areaOfCMRep << endl;
  cout << "area of Image = " << areaOfImage << endl;
  cout << "area of overlap before shift = " << overlapOri << endl;
  cout << "area of overlap = " << overlap << endl;
  return overlap/(areaOfCMRep+areaOfImage-overlap);
}
