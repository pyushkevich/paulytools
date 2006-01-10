#include "OptimizationTermsBSpline.h"
#include <iostream>
#include "optima.h"
#include "ConjugateGradientMethod.h"

using namespace std;

double optimize(CMRep2DOptimizationProblem& fcm, Vector& x, Vector& xGrad, double stepSize, double imageBlurVariance, char *inputFileName);
//void addToVector(Vector& currentX, Vector& x);

int main(int argc, char *argv[]) {

 
  if (argc < 4) {
    cerr <<"Not enough arguments" <<endl;
    cerr <<"Usage: " << argv[0] << "inputImageFile dim"<<endl;
    exit(1);
  }
  char *inputImageFile = argv[1];
  int dim = atoi(argv[2]);
  double scaleFactor = atof(argv[3]);
  double stepSize = atof(argv[4]);

  // cerr << "Begin to optimize coeffs for image " << inputImageFile << endl; 

  const double x0[8*3] = {
    /* 56.0068, 58.6106, 68.3392, 62.8632, 49.797, 37.2775,
    23.3805, 9.1962, -5.46414, -19.7429, -33.4751, -46.5852,
    -59.3993, -69.4148, -69.686, -68.479,

    -26.3054, -23.613, -17.3027, -2.97595, 4.27851, 11.3017, 
    16.0621,19.1998, 19.0956, 16.2316, 11.3509, 5.17381,
    -1.44601, -12.3681, -22.6114,-26.6637,

    6.40185, 0.462407, -1.44657, 0.461113, -0.0267031, 0.110174, -0.0507504, 
    0.0484413, -0.113689, 0.0806307, -0.0264758, 0.385522, -0.275604, 0.268138, 
    -3.51883, 14.1603
    };
  
    //  200.0068, 100.8632, 50.3805, -19.7429, -80.3993, -200.479,
     56.0068, 62.8632, 23.3805, -19.7429, -59.3993, -68.479,
        -26.3054, -2.97595, 16.0621, 16.2316, -1.44601, -26.6637,
     //-80.3054, -20.97595, 16.0621, 16.2316, -20.44601, -80.6637,
    1.0, -1.0, -0.1, 0.1,  -1.0, 0.5   
    };*/
   -26.038781346731948, 79.80877917668593, 57.216606149055146,20.572369752040505, -17.053366355024266, -52.75306161265423, -82.16512590218723, -32.796713791079426,
    -30.18894204154545, -28.927657211699714, 3.8075596689297297, 19.069372026847958, 19.434008249360303, 3.4943536404368625, -20.64086668985507, -55.1118916378337,
   -0.004,   -0.004,   -0.004,   -0.004,   -0.004,   -0.004,   -0.004,   -0.004};
  Vector x1;
  x1.setSize(8*3);
  for (int i = 0; i < x1.size(); ++i) {
    x1(i) = x0[i];
  }
  Vector x1Grad;
  

// generate a CMRep2DOptimizationProblem class object and test it

  CMRep2DOptimizationProblem fcm(inputImageFile, 1.0, dim);
  
  fcm.setGradientScaleFactor(1.0, 1.0, scaleFactor);
  fcm.setStepEpsilon(1.0e-3,1.0e-3,1.0e-4);    
  
  double  fx1 =  fcm.evaluate(x1);
    
  cerr <<"f(x1) = "<< fx1 << endl;
  //  for (int i =0; i < x1Grad.size(); ++i) {
  //  cerr << "x1Grad(" << i << ") = " << x1Grad(i) << endl;
  // }
      
  // start from optimize 8*3 coeffs
  //  cerr << "start from " << x1.size()/3 << "*3 coeffs ......" << endl;
  double iteration = optimize(fcm, x1, x1Grad,stepSize, 0.1, inputImageFile);
  double OverlapRatio1 = fcm.areaOverlapRatio(x1);
 

  // output the optimization results for 8*3 coeffs
  
   cerr << "x1Grad = {" ;
   for (int i = 0; i < x1Grad.size()-1; ++i) {
    cerr << x1Grad(i) << ",";
    }
    cerr << x1Grad(x1Grad.size()-1) << "}" << endl; 
  
  cerr << "x1" << inputImageFile << " = {" ;
  for(int i = 0; i < x1.size()-1; ++i) {
    cerr <<  x1(i) << ", ";
  }
  cerr << x1(x1.size()-1) << "};" << endl ;
  cerr << "overlapRatio1" << inputImageFile << " = " << OverlapRatio1 << ";" << endl ;

 fx1 =  fcm.computeOneJet(x1,x1Grad);
 cout << "fx(x1) = " << fx1 << endl;
 double dx = 1.0e-3;
 for (int i = 0; i < x1.size(); ++i) {
   x1(i) += dx;
   cout << "fx(" << i << "+dx) = " << fcm.evaluate(x1) << ",";
   x1(i) -= 2.0*dx;
   cout << "fx(" << i << "-dx) = " << fcm.evaluate(x1) << ";";
   x1(i) += dx;
 }
 cout << endl;
  Vector bx;
  Vector by;
  bool valid = fcm.getBoundary(x1, bx, by);
  if ( valid) {
    cerr << "boundaryPoints"<< x1.size()/3 << inputImageFile << " = {" ;
    for( int i = 0; i < bx.size() -1; ++i) {
      cerr << "{" << bx(i) << "," << by(i) << "}, " ;
    } 
    cerr << "{" << bx(bx.size() - 1) << "," << by(by.size() - 1) << "}};" << endl << endl ;
  }

  /*    
  // now fix the value we get from previous optimization, add (32)*3 coeffs to optimize
  fcm.setFixedX(x1);
  fcm.setGradientScaleFactor(1.0, 1.0, scaleFactor);
  Vector x2;
  x2.setSize(32*3);
  for (int i = 0; i < x2.size(); ++i) {
    x2(i) = 0.0;
  }
  Vector x2Grad;
  //  cerr << "then add " << x2.size()/3 << "*3 coeffs ......" << endl;
  iteration += optimize(fcm, x2, x2Grad, stepSize , 0.1, inputImageFile);
  double  OverlapRatio2 = fcm.areaOverlapRatio(x2);
  
  // output the optimization results for 8*3 + 32*3 coeffs

  /*
  cerr << "x2Grad = {" ;
  for (int i = 0; i <x2Grad.size()-1; ++i) {
    cerr << x2Grad(i) << ",";
  }
  cerr << x2Grad(x2Grad.size()-1) << "}" << endl; 

  cerr << "totalIteration" << inputImageFile << " = " << iteration << ";" << endl << endl;
  cerr << "x2" << inputImageFile << " = {" ;
  for(int i = 0; i < x2.size()-1; ++i) {
    cerr <<  x2(i) << ", ";
  }
  cerr << x2(x2.size()-1) << "};" << endl << endl;
  cerr << "overlapRatio2" << inputImageFile << " = " << OverlapRatio2 << ";" << endl << endl;
  
 //generate boundary point and output
  valid = fcm.getBoundary(x2, bx, by);
    if ( valid) {
      cerr << "boundaryPoints"<< x2.size()/3 << inputImageFile << " = {" ;
      for( int i = 0; i < bx.size() -1; ++i) {
	cerr << "{" << bx(i) << "," << by(i) << "}, " ;
      } 
      cerr << "{" << bx(bx.size() - 1) << "," << by(by.size() - 1) << "}};" << endl << endl ;
    }
*/ 
  return 0;
};


double optimize(CMRep2DOptimizationProblem& fcm, Vector& x, Vector& xGrad, double stepSize, double imageBlurVariance, char *inputFileName) {
  fcm.setBlurVariance(1.0);
  ConjugateGradientMethod CG(fcm,x);
  CG.setTolerance(1.0e-4);
  CG.setStepSize(stepSize);
  int xSize = x.size()/3;  
  double objectiveF;
  bool finished = false;
  int  counter = 0;
  for(counter;!finished; ++counter) {
    CG.performIteration();
    finished = CG.isFinished();
    objectiveF = CG.getBestEverValue();
    cerr << "optimizing " << xSize << "*3 coeffs, objective function value= " << objectiveF << endl;
  }
  x = CG.getBestEverX();
  fcm.setGradientScaleFactor(1.0e-4, 1.0e-4, 1.0e-4);
  fcm.setBlurVariance(imageBlurVariance);
  ConjugateGradientMethod CG2(fcm,x);
  CG2.setTolerance(1.0e-4);
  CG2.setStepSize(stepSize);
  finished = false;
  for(counter;!finished; ++counter) {
    CG2.performIteration();
    finished = CG2.isFinished();
    objectiveF = CG2.getBestEverValue();
    cerr << "optimizing " << xSize << "*3 coeffs, objective function value= " << objectiveF << endl;
  }
  x = CG2.getBestEverX();
  double fx =  fcm.computeOneJet(x,xGrad);
  cerr << "iteration" << xSize << inputFileName << " = " << counter << ";" << endl << endl; 
  cerr << "objectiveFunctionValue" << xSize << inputFileName << " = " << fx << ";" << endl << endl; 

  return counter;

} 
