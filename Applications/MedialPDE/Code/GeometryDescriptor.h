#ifndef __GeometryDescriptor_h_
#define __GeometryDescriptor_h_

#include <iostream>

struct GeometryDescriptor
{
  double xCovariantTensor[2][2];
  double xContravariantTensor[2][2];
  double xChristoffelFirst[2][2][2];
  double xChristoffelSecond[2][2][2];

  // The determinant of the covariant tensor and its inverse
  double g, gInv;

  // Initialize the descriptor using a Jet
  void SetJet(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv);

  // Dump information out
  void PrintSelf(std::ostream &str);
};

#endif
