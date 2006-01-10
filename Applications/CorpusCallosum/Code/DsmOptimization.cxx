#include "OptimizationTermsBSpline.h"
#include <iostream>
#include "/mnt/data1/huiz/research/eclipse/workspace/Mala_3D/src/numerics/numerics.h"
#include "/mnt/data1/huiz/research/eclipse/workspace/Mala_3D/src/numerics/Function.h"
using namespace std;

class numerics::Function;
double compute8 (const double p[]);
double compute32 (const double p[]);
CMRep2DOptimizationProblem *fcm;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr <<"Not enough arguments" <<endl;
    cerr <<"Usage: " << argv[0] << "inputImageFile dim"<<endl;
    exit(1);
  }
  char *inputImageFile = argv[1];
  int dim = atoi(argv[2]);
  //fcm is defined golabally
  fcm = new  CMRep2DOptimizationProblem(inputImageFile, 1.0, dim);
  double x1[8*3] = {
   -26.038781346731948, 79.80877917668593, 57.216606149055146,20.572369752040505, -17.053366355024266, -52.75306161265423, -82.16512590218723, -32.796713791079426,
    -30.18894204154545, -28.927657211699714, 3.8075596689297297, 19.069372026847958, 19.434008249360303, 3.4943536404368625, -20.64086668985507, -55.1118916378337,
   -0.04,   -0.04,   -0.04,   -0.04,   -0.04,   -0.04,   -0.04,   -0.04};

  double ftol = 1.0e-4;
  numerics::Function func8(8*3, compute8);
  numerics::DirectionSetMinimizer dsm;
  fcm->setBlurVariance(1.0);
  double min8 = dsm.run(x1, ftol, func8);
  fcm->setBlurVariance(0.1);
  min8 = dsm.run(x1, ftol, func8);
 
 //get the result
  Vector v1;
  v1.setSize(8*3);
  for (int i = 0; i < v1.size(); ++i) {
    v1(i) = x1[i];
  }
  double OverlapRatio1 = fcm->areaOverlapRatio(v1);
  Vector bx1, by1;
  bool valid1 = fcm->getBoundary(v1, bx1, by1);

  // now fix the value we get from previous optimization, add (32)*3 coeffs to optimize
  fcm->setFixedX(v1);
  double x2[32*3] = {
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,

    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  numerics::Function func32(32*3, compute32);
  fcm->setBlurVariance(1.0); 
  double min32 = dsm.run(x2, ftol, func32);
  fcm->setBlurVariance(0.1);
  min32 = dsm.run(x2, ftol, func32);
 
  //get the result of 32*3 coeffs
  Vector v2;
  v2.setSize(32*3);
  for (int i = 0; i < v2.size(); ++i) {
    v2(i) = x2[i];
  }
  double  OverlapRatio2 = fcm->areaOverlapRatio(v2);
  Vector bx2, by2;
  bool valid2 = fcm->getBoundary(v2, bx2, by2);


  // output the optimization results for 8*3 coeffs
  cerr << "v1[[" << inputImageFile << "]] = {" ;
  for(int i = 0; i < v1.size()-1; ++i) {
    cerr <<  v1(i) << ", ";
  }
  cerr << v1(v1.size()-1) << "};" << endl ;
  cerr << "overlapRatio1[[" << inputImageFile << "]] = " << OverlapRatio1 << ";" << endl ;

  if ( valid1 ) {
    cerr << "boundaryPoints"<< v1.size()/3 << "[[" << inputImageFile << "]] = {" ;
    for( int i = 0; i < bx1.size() -1; ++i) {
      cerr << "{" << bx1(i) << "," << by1(i) << "}, " ;
    } 
    cerr << "{" << bx1(bx1.size() - 1) << "," << by1(by1.size() - 1) << "}};" << endl << endl ;
  }


  // output the optimization results for 8*3 + 32*3 coeffs
  cerr << "v2[[" << inputImageFile << "]] = {" ;
  for(int i = 0; i < v2.size()-1; ++i) {
    cerr <<  v2(i) << ", ";
  }
  cerr << v2(v2.size()-1) << "};" << endl << endl;
  cerr << "overlapRatio2[[" << inputImageFile << "]] = " << OverlapRatio2 << ";" << endl << endl;
  
  if ( valid2) {
    cerr << "boundaryPoints"<< v2.size()/3 << "[[" << inputImageFile << "]] = {" ;
    for( int i = 0; i < bx2.size() -1; ++i) {
      cerr << "{" << bx2(i) << "," << by2(i) << "}, " ;
    } 
    cerr << "{" << bx2(bx2.size() - 1) << "," << by2(by2.size() - 1) << "}};" << endl << endl ;
  }
 
  return 0;
};


double compute8 (const double p[]) {
  Vector x;
  x.setSize(8*3);
  for (int i = 0; i < 8*3; ++i){
    x(i) = p[i];
  }
  return fcm->evaluate(x);
}


double compute32 (const double p[]) {
  Vector x;
  x.setSize(32*3);
  for (int i = 0; i < 32*3; ++i){
    x(i) = p[i];
  }
  return fcm->evaluate(x);
}
