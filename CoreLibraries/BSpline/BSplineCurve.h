#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include <cmath>
#include <vector>
#include <matrix.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

using namespace std;

template <class TReal = double, unsigned int VOrder = 3, unsigned int VJetOrder = 2>
class BSplineKnotList
{
public:
  // Type representing the basis function multiplier for each spline patch
  typedef vnl_vector_fixed<TReal, (VOrder+1)> BasisVectorType;
  typedef vnl_vector<TReal> VectorType;

  /**
   * Create a new vector of knots
   * z: Number of zero-knots at each end, between 2 and VOrder+1 
   */
  BSplineKnotList(unsigned int m);

  /** Destructor */
  ~BSplineKnotList() {}

  /**
   * This method evaluates the basis function and its derivatives
   * for a patch starting at a given knot index
   * 
   * i: Knot index (first relevant knot starts at VOrder + 1)
   * u: Value of the parameter, between 0 and 1
   * W: Array of VJetOrder vectors of size vOrder + 1 that will receive the derivatives
   * This function depends on the knot sequence.
   * When there are only two zero knots at ends, the 
   * function is matrix multiplication. 
   */
  void ComputeBasisJet(unsigned int i,TReal u,BasisVectorType W[VJetOrder+1]);

  /**
   * Return the u or v value for a given knot index
   */
  TReal GetParameterValueAtKnot(unsigned int knotIndex) const 
    {
    return m_Knots[knotIndex];
    }

  /**
   * Get knot value for a given t
   */
  unsigned int GetKnotAtParameterValue(TReal u);

  /**
   * Create a knot index sequence for an ordered sequence of u values
   */
  vnl_vector<unsigned int> GetKnotSequence(unsigned int nPoints, const TReal *points);

private:
  // Number of knots
  unsigned int m_NumberOfKnots;

  // Number of control points
  unsigned int m_NumberOfControlPoints;

  // Number of patches on the spline
  unsigned int m_NumberOfPatches;
 
  // Piegl's m (n. controls - 1)
  unsigned int m_PieglM;

  // A vector of knots
  vnl_vector<TReal> m_Knots;
};

/**
 * A definition of a BSplineCurve in N-space. Each control point is a vector
 * of TReal of dimension VDimension. The spline is of order 3 and can be
 * evaluated up to the VJetOrder'th derivative 
 */
template <unsigned int VDimension, class TReal = double, 
  unsigned int VOrder = 3, unsigned int VJetOrder = 2>
class BSplineCurve
{
public:
  typedef BSplineKnotList<TReal, VOrder, VJetOrder> KnotList;
  typedef typename KnotList::BasisVectorType BasisVector;

  // Vector used to represent a control point's location in N-space
  typedef vnl_vector_fixed<TReal,VDimension> Point;

  // A matrix used to represent the span of VOrder+1 control points
  typedef vnl_matrix_fixed<TReal, VDimension, VOrder+1> ControlMatrix;

  /** 
   * Create a new B-Spline Curve. By default this generates a minimal
   * possible curve of one dimension 
   */
  BSplineCurve() : m_KnotList(NULL)
    {
    SetNumberOfControlPoints(VOrder + 1);
    }

  /** Destructor */
  ~BSplineCurve()
    {
    delete m_KnotList;
    }

  /** Set the number of control points for the curve. This method discards
   * previously set control points */
  void SetNumberOfControlPoints(unsigned int nPoints)
    {
    // Resize the control point array
    m_Controls.resize(nPoints, Point(0.0));
    m_ControlGroups.resize(nPoints,ControlMatrix(0.0));

    // Reinitialize the knot vector
    if(m_KnotList) delete m_KnotList;
    m_KnotList = new KnotList(nPoints - 1);

    // Store the number of contol points
    m_NumberOfControlPoints = nPoints;
    }

  /** Set the N-th control point */
  void SetControlPoint(unsigned int iPoint, const Point &point)
    {
    // Check input
    assert(iPoint < m_NumberOfControlPoints);
    
    // Set the control point's value
    m_Controls[iPoint] = point;

    // Set the neighbor information for adjacent control points
    for(unsigned int d = 0; d <= VOrder && d <= iPoint; d++)
      {
      m_ControlGroups[iPoint - d].set_column(d, point);
      }
    }

  /** Get the N-th control point */
  const Point &GetControlPoint(unsigned int iPoint)
    {
    return m_Controls[iPoint];
    }

  /** Get the number of control points */
  unsigned int GetNumberOfControlPoints()
    {
    return m_Controls.size();
    }

  /** Interpolate the spline or its derivative at a parameter value u. This is a 
   * convenience method that is a bit too slow to use in practice */
  Point Evaluate(TReal u, unsigned int order)
    {
    // Get the knot index corresponding to u
    unsigned int iKnot = m_KnotList->GetKnotAtParameterValue(u);

    // Compute the basis for the knot
    BasisVector basis[VJetOrder];
    m_KnotList->ComputeBasisJet(iKnot,u,basis);

    // Evaluate the basis for the requested order
    return m_ControlGroups[iKnot] * basis[order];
    }

  /** A structure used in the interpolation grid */
  struct GridPoint 
    { 
    TReal u;
    unsigned int iKnot;
    BasisVector basis[VJetOrder];
    };

  // Interpolation grid typedef 
  typedef std::vector<GridPoint> EvaluationGrid;
   
  /** Create an 'interpolation grid' that can be used to quickly interpolate
   * the spline at a number of points. The interpolation grid becomes useless
   * if the number of control points changes after it's creation */
  void CreateEvaluationGrid(unsigned int nPoints, const TReal *points, EvaluationGrid &grid);

  /** Create a uniformely spaced evaluation grid */
  void CreateUniformEvaluationGrid(unsigned int nPoints, EvaluationGrid &grid);

  /** Fit the spline to a data matrix */
  void FitToPoints(unsigned int nPoints, const Point *xPoints);
    
  /** Evaluate the spline at a given order over the evaluation grid */
  Point EvaluateGridPoint(const EvaluationGrid &grid, unsigned int iPoint, unsigned int iJetOrder)
    {
    // Check the index
    assert(iPoint < grid.size());
    
    // Get the grid point
    const GridPoint &point = grid[iPoint];

    // Compute the interpolation
    return m_ControlGroups[point.iKnot] * point.basis[iJetOrder];
    }

private:
  // Number of the control points
  unsigned int m_NumberOfControlPoints;

  // The knot list
  KnotList *m_KnotList;

  // List of control points
  vector<Point> m_Controls;

  // List of control point neighborhood matrices
  vector<ControlMatrix> m_ControlGroups;
};

#include "BSplineCurve.txx"

#endif
