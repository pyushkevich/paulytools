#include "OptimizationTerms.h"
#include <iostream>
#include "optima.h"
#include "ConjugateGradientMethod.h"

using namespace std;

double optimize(CMRep2DOptimizationProblem& fcm, Vector& x, Vector& xGrad);
void addToVector(Vector& currentX, Vector& x);

int main(int argc, char *argv[]) {

  if (argc < 3) {
    cerr <<"Not enough arguments" <<endl;
    cerr <<"Usage: " << argv[0] << "inputImageFile imageBlurVariance dim"<<endl;
    exit(1);
  }
  const double p0[51] = {
	-28.2364,
	60.1971, -0.377662,
	-4.78251, 4.24725,
	4.21383, 1.00549, -0.33139, -1.87177,
	0.568533, -0.603949, 0.154689, -0.1541, 0.072736, -0.123597, 0.104211, -0.666733,
	
	-26.1098,
	29.7462, -0.00113348,
	2.65794, 2.91456,
	1.34698, 0.12066, 0.22392, 2.53658,
	-0.783693, 0.131142, 0.28515, 0.0791865, -0.0151148, 0.130991, -0.0664204, 0.20153,


	-0.2780468875676537,
	0.8430449632231002, -0.4481581004039171,
    -0.3168084491489501, 0.4390913581310029,
	-0.5411271156693166, -0.09268973864155224, -0.09077243805937582, 0.1103090808314945,
	-0.6752267991252191, -0.02694546649788497, -0.05732010382287213, -0.03807655204148059, -0.1325252333762741, -0.21076001234790967, -0.31639484874831914, -0.39152479130559137
};
/*
	-0.025,
	0.00, 0.00,
	0.00, 0.00,
	0.00, 0.00, 0.00, 0.00,
	0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00
	};
*/

 const double x0[27] = {
   0.0,
   0.0, 0.0,
   0.0, 0.0,
   0.0, 0.0, 0.0, 0.0,

   0.0,
   0.0, 0.0,
   0.0, 0.0,
   0.0, 0.0, 0.0, 0.0,

   0.0,
   0.0, 0.0,
   0.0, 0.0,
   0.0, 0.0, 0.0, 0.0,

 };
  // {-5.08345, -0.864421, -1.93286, 1.27896, 0.29162, 1.28341, -0.130888, -0.434299, -1.0902, 1.27523, -1.90642, 0.294781, -0.824298, -2.75131, -0.922801, -1.89178, 0.309253, 1.5994, -0.0189866, 0.020699, 0.00368445, 0.032449, -0.0166652, -0.105217, -0.0727629, 0.0110231, 0.0820438};

  Vector tx;
  tx.setSize(51);
  
   for (int i = 0; i < tx.size(); ++i) {
   tx(i) = p0[i];
  // cout << "tx(" << i << ") = " << tx(i) << endl;
   }
  
  Vector x;
  x.setSize(9*3);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = 0.0;
  }
  Vector xGrad;
 
  const char *inputImageFile = argv[1];
  double imageBlurVariance = atof(argv[2]);
  int dim = atoi(argv[3]);

  // generate a CMRep2DOptimizationProblem class object and test it
  CMRep2DOptimizationProblem fcm(inputImageFile, imageBlurVariance, dim, tx);
  fcm.setGradientScaleFactor(1.0, 1.0, 1.0);
  fcm.setStepEpsilon(1.0e-5,1.0e-5,1.0e-6);    
  double  fx =  fcm.computeOneJet(x,xGrad);
  cout <<"computeOneJet(x) = "<< fx << endl;
  for (int i =0; i < xGrad.size(); ++i) {
    cout << "xGrad(" << i << ") = " << xGrad(i) << endl;
  }
  
  // start from optimize 9*3 coeffs
  double iteration = optimize(fcm, x, xGrad);
  double currentOverlapRatio = fcm.areaOverlapRatio(x);
  Vector  currentX = x;
  Vector currentXGrad = xGrad;

  // output the optimization results for 9*3 coeffs
  cout << " total iteration = " << iteration << endl;
  cout << "xGrad = {" ;
  for (int i = 0; i < currentXGrad.size()-1; ++i) {
    cout << currentXGrad(i) << ",";
  }
  cout << currentXGrad(currentXGrad.size()-1) << "}" << endl; 
  cout << " x = {" ;
  for(int i = 0; i < currentX.size()-1; ++i) {
    cout <<  currentX(i) << ", ";
  }
  cout << currentX(currentX.size()-1) << "}" << endl;
  cout << "overlap Ratio = " << currentOverlapRatio << endl;

  // now fix the value we get from previous optimization, add 8*3 coeffs to optimize
  fcm.setFixedX(currentX);
  fcm.setGradientScaleFactor(1.0, 1.0, 1.0);
  x.setSize(8*3);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = 0.0;
  }
  iteration += optimize(fcm, x, xGrad);
  currentOverlapRatio = fcm.areaOverlapRatio(x);
  addToVector(currentX, x);
  addToVector(currentXGrad, xGrad);
  
  // output the optimization results for 8*3 + 9*3 coeffs
  cout << " total iteration = " << iteration << endl;
  cout << "xGrad = {" ;
  for (int i = 0; i < currentXGrad.size()-1; ++i) {
    cout << currentXGrad(i) << ",";
  }
  cout << currentXGrad(currentXGrad.size()-1) << "}" << endl; 
  cout << " x = {" ;
  for(int i = 0; i < currentX.size()-1; ++i) {
    cout <<  currentX(i) << ", ";
  }
  cout << currentX(currentX.size()-1) << "}" << endl;
  cout << "overlap Ratio = " << currentOverlapRatio << endl;

  // now fix the value we get from previous optimization, add 16*3 coeffs to optimize
  fcm.setFixedX(currentX);
  fcm.setGradientScaleFactor(1.0, 1.0, 1.0);
  x.setSize(16*3);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = 0.0;
  }
  fx =  fcm.computeOneJet(x,xGrad);
  cout << "fx = " << fx << endl;  
  iteration += optimize(fcm, x, xGrad);
  currentOverlapRatio = fcm.areaOverlapRatio(x);
  addToVector(currentX, x);
  addToVector(currentXGrad, xGrad);
  
  // output the optimization results for 16*3 + 8*3 + 9*3 coeffs
  cout << " total iteration = " << iteration << endl;
  cout << "xGrad = {" ;
  for (int i = 0; i < currentXGrad.size()-1; ++i) {
    cout << currentXGrad(i) << ",";
  }
  cout << currentXGrad(currentXGrad.size()-1) << "}" << endl; 
  cout << " x = {" ;
  for(int i = 0; i < currentX.size()-1; ++i) {
    cout <<  currentX(i) << ", ";
  }
  cout << currentX(currentX.size()-1) << "}" << endl;
  cout << "overlap Ratio = " << currentOverlapRatio << endl;
 
  
 /* //generate boundary point and output
    Vector bx;
    Vector by;
    bool valid = fcm.getBoundary(x,bx, by);
    if ( valid) {
      cout << " boundaryPoint = {" ;
      for( int i = 0; i < bx.size() - 1; ++i) {
	cout << "{" << bx(i) << "," << by(i) << "}, " ;
      } 
      cout << "{" << bx(bx.size() - 1) << "," << by(by.size() - 1) << "}}" << endl ;
    }
    */
  return 0;
};


double optimize(CMRep2DOptimizationProblem& fcm, Vector& x, Vector& xGrad) {
  ConjugateGradientMethod CG(fcm,x);
  CG.setTolerance(1.0e-6);
  CG.setStepSize(1.0e-4);
  int xSize = x.size()/3;  
  double objectiveF;
  bool finished = false;
  int  counter = 0;
  for(counter;!finished; ++counter) {
    CG.performIteration();
    finished = CG.isFinished();
    objectiveF = CG.getBestEverValue();
    cout << "optimizing " << xSize << "*3 coeffs, objective function value= " << objectiveF << endl;
  }
  x = CG.getBestEverX();  
  double fx =  fcm.computeOneJet(x,xGrad);
  cout << "optimized " << xSize << "*3 coeffs, iteration = " << counter << endl; 
  cout << "optimized " << xSize << "*3 coeffs, objective Function value = " << fx << endl;
   
  return counter;

} 

void addToVector(Vector& currentX, Vector& x) {
  Vector tmpV = currentX;
  currentX.setSize(tmpV.size() + x.size());
  for (int i = 0; i < tmpV.size()/3; ++i) {
    currentX(i) = tmpV(i);
    currentX(i + currentX.size()/3) = tmpV(i + tmpV.size()/3);
    currentX(i + 2*currentX.size()/3) = tmpV(i + 2*tmpV.size()/3);
  }
  for (int i = 0; i < x.size()/3; ++i) {    
    currentX(tmpV.size()/3 + i) = x(i);
    currentX(tmpV.size()/3 + currentX.size()/3 + i) = x(i + x.size()/3);
    currentX(tmpV.size()/3 + 2*currentX.size()/3 + i) = x(i + 2*x.size()/3);
  }
}
