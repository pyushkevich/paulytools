#include "BSplineCurve.h"
#include <iostream>

using namespace std;

int main(void)
{
  // Pick the type of b-spline
  typedef BSplineCurve <3,double,3> CubicCurve3D;
  typedef BSplineCurve <2,float,2> QuadraticCurve2D;
  
  // Create a 3D b-spline object
  CubicCurve3D curve3d;
  
  // Initialize the control points for the curve
  curve3d.SetNumberOfControlPoints(5);
  curve3d.SetControlPoint(0, CubicCurve3D::Point(3.0, 2.0,  5.0));
  curve3d.SetControlPoint(1, CubicCurve3D::Point(4.0, 4.0,  4.0));
  curve3d.SetControlPoint(2, CubicCurve3D::Point(5.0, 8.0,  3.0));
  curve3d.SetControlPoint(3, CubicCurve3D::Point(6.0, 16.0, 2.0));
  curve3d.SetControlPoint(4, CubicCurve3D::Point(7.0, 32.0, 1.0));

  // Interpolate the curve at a given position
  cout << "Position at u = 0.5 \t" << curve3d.Evaluate(0.5,0) << endl;

  // Create an interpolation grid for the curve
  CubicCurve3D::EvaluationGrid grid3d;
  curve3d.CreateUniformEvaluationGrid(20, grid3d);

  // Evaluate the grid
  CubicCurve3D::Point *xPoints = new CubicCurve3D::Point[grid3d.size()];
  for(unsigned int i=0;i<grid3d.size();i++)
    {
    xPoints[i] = curve3d.EvaluateGridPoint(grid3d,i,0);
    cout << "Position at u = " << grid3d[i].u << "\t" << xPoints[i] << endl;
    }
  
  // Refit the spline using a different number of control points
  curve3d.SetNumberOfControlPoints(4);
  curve3d.FitToPoints(grid3d.size(),xPoints);

  // Print the new control points
  for(unsigned int j=0;j<curve3d.GetNumberOfControlPoints();j++)
    {
    cout << "New control point " << j << "\t" << curve3d.GetControlPoint(j);
    }
}
