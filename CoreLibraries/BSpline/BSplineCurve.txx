#include "bspline.h"
#include <vnl/algo/vnl_svd.h>
#include <iostream>

template <class TReal, unsigned int VOrder, unsigned int VJetOrder>
BSplineKnotList<TReal,VOrder,VJetOrder>
::BSplineKnotList(unsigned int m) 
{
  unsigned int i;

  // Save the m
  m_PieglM = m; 
  m_NumberOfControlPoints = m + 1;
  m_NumberOfKnots = m_NumberOfControlPoints + VOrder + 1;
  m_NumberOfPatches = m_NumberOfControlPoints - VOrder;

  // Create the knot array (by uniform interpolation)
  m_Knots.set_size(m_NumberOfKnots);

  // Set the knot sequence
  TReal step = 1.0 / m_NumberOfPatches;
  for(unsigned int i = 0;i < m_NumberOfKnots; i++)
    {
    if(i <= VOrder)
      m_Knots[i] = 0.0;
    else if (i >= m_NumberOfKnots - (VOrder + 1))
      m_Knots[i] = 1.0;
    else
      m_Knots[i] = step * (i - VOrder);
    }
}

template <class TReal, unsigned int VOrder, unsigned int VJetOrder>
void
BSplineKnotList<TReal,VOrder,VJetOrder>
::ComputeBasisJet(unsigned int i,TReal u,BasisVectorType N[]) 
{
  // We use the variables p and n because we want the code to match the code in Piegl
  const unsigned int p = VOrder;
  const unsigned int n = VJetOrder;
  TReal *U = m_Knots.data_block();

  // Left and right arrays are not dynallocated for effifiency.  Don't pass p >= 10
  TReal ndu[p+1][p+1],a[2][p+1],left[p+1],right[p+1];

  TReal saved,temp;
  int j,k,j1,j2,r;

  ndu[0][0]=1.0f;
  for (j=1; j<=p; j++ )
    {
    left[j]=u-U[i+1-j];
    right[j]=U[i+j]-u;
    saved=0.0f;

    for ( r=0; r<j; r++ )
      {
      ndu[j][r]=right[r+1]+left[j-r];
      temp=ndu[r][j-1]/ndu[j][r];

      ndu[r][j]=saved+right[r+1]*temp;
      saved=left[j-r]*temp;
      }
    ndu[j][j]=saved;
    }

  for ( j=0; j<=p; j++ )
    N[0][j]=ndu[j][p];  // basis function

  // compute derivatives
  for ( r=0; r<=p; r++ )
    {
    int s1=0, s2=1;  // alternate rows in a
    a[0][0]=1.;

    // compute k'th derivative
    for ( k=1; k<=n; k++ )
      {
      TReal d=0.0;
      int rk=r-k, pk=p-k;
      if ( r>=k )
        {
        a[s2][0]=a[s1][0]/ndu[pk+1][rk];
        d=a[s2][0]*ndu[rk][pk];
        }

      j1 = rk >= -1 ? 1 : -rk;
      j2 = (r-1<=pk) ? k-1 : p-r;

      for ( j=j1; j<=j2; j++ )
        {
        a[s2][j]=(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
        d+=a[s2][j]*ndu[rk+j][pk];
        }
      if ( r<=pk )
        {
        a[s2][k]= -a[s1][k-1]/ndu[pk+1][r];
        d+=a[s2][k]*ndu[r][pk];
        }
      N[k][r]=d;
      j=s1; s1=s2; s2=j;  // switch rows
      }
    }

  r=p;
  for ( k=1; k<=n; k++ )
    {
    for ( j=0; j<=p; j++ )
      N[k][j]*=r;
    r*=p-k;
    }    
}

template <class TReal, unsigned int VOrder, unsigned int VJetOrder>
vnl_vector<unsigned int>
BSplineKnotList<TReal,VOrder,VJetOrder>
::GetKnotSequence(unsigned int nPoints, const TReal *x) 
{
  vnl_vector<unsigned int> z(nPoints);
  
  int ki = 0;
  for (unsigned int i=0;i<nPoints;i++)
    {
    while (GetParameterValueAtKnot(ki+1) <= x[i] && x[i] < 1.0)
      ki++;
    z[i] = ki;
    }

  return z;
}

template <class TReal, unsigned int VOrder, unsigned int VJetOrder>
unsigned int
BSplineKnotList<TReal,VOrder,VJetOrder>
::GetKnotAtParameterValue(TReal u) 
{
  unsigned int n = m_NumberOfControlPoints - 1;
  if (u >= m_Knots[n+1]) 
    return n;
  int lo = 3,hi = n+1,mid = (lo+hi)/2;
  while (u < m_Knots[mid] || u>=m_Knots[mid+1])
    {
    if (u < m_Knots[mid]) hi = mid;
    else lo = mid;
    mid = (lo+hi)/2;
    }
  return mid;
}

template <unsigned int VDimension, class TReal, unsigned int VOrder, unsigned int VJetOrder>
void 
BSplineCurve<VDimension,TReal,VOrder,VJetOrder>
::FitToPoints(unsigned int nPoints, const Point *points) 
{
  unsigned int i;

  // M is the index of the last point 
  unsigned int m = nPoints - 1;

  // N is the index of the last control point
  unsigned int n = m_NumberOfControlPoints - 1;

  // D is the length of the curve Q
  TReal len = 0;
  for (i=1;i<=m;i++)
    len += (points[i] - points[i-1]).two_norm();

  // Construct an array of u-values for the Q points
  typename KnotList::VectorType uk(m+1);
  uk(0) = 0; uk(m) = 1.0;
  for (i=1;i<m;i++)
    {
    uk(i) = uk(i-1) + (points[i] - points[i-1]).two_norm() / len;
    }

  // Find knot indices of each uk
  vnl_vector<unsigned int> uki = m_KnotList->GetKnotSequence(m+1,uk.data_block());

  // Constuct a 'band' matrix of N_i,p(u_j)
  typedef vnl_matrix<TReal> MatrixType;
  MatrixType BN(m+1,n+1,0.0);
  for (int r=0;r<=m;r++)
    {
    BasisVector W[VJetOrder];
    m_KnotList->ComputeBasisJet(uki[r],uk[r],W);
    unsigned int c0 = uki[r] - VOrder;
    unsigned int c1 = c0 + VOrder < n ? c0 + VOrder : n;
    for (unsigned int c=c0;c<=c1;c++)
      {
      BN(r,c) = W[0][c-c0];
      }
    }

  // Construct matrices N,R
  MatrixType R(n-1,VDimension,0.0);

  // R is computed as a sum
  typedef vnl_vector<double> VectorType;
  for (unsigned int c=1;c<=n-1;c++)
    {
    for (unsigned int r=1;r<=m-1;r++)
      {
      VectorType L = points[r] - points[0] * BN(r,0) - points[m] * BN(r,n);
      R.set_row(c-1,R.get_row(c-1) + L * BN(r,c));
      }
    }

  // N is a chunk of BN
  MatrixType N = BN.extract(m-1,n-1,1,1);  
  MatrixType NTN = N.transpose() * N;

  // Solve for P
  vnl_svd<double> svd(NTN);
  MatrixType R1 = svd.solve(R);

  // Set the control points
  SetControlPoint(0, points[0]);
  for (unsigned int i=1;i<n;i++)
    SetControlPoint(i, R1.get_row(i-1));
  SetControlPoint(n, points[m]);

  // Done!
  std::cout << "fit " << std::endl;
}

template <unsigned int VDimension, class TReal, unsigned int VOrder, unsigned int VJetOrder>
void
BSplineCurve<VDimension,TReal,VOrder,VJetOrder>
::CreateEvaluationGrid(unsigned int nPoints, const TReal *points, EvaluationGrid &grid)
{
  // Get the knot indices for the points
  vnl_vector<unsigned int> uk = m_KnotList->GetKnotSequence(nPoints,points);

  // Create the grid
  grid.resize(nPoints);
  for (unsigned int i=0;i<nPoints;i++)
    {
    grid[i].u = points[i];
    grid[i].iKnot = uk[i];
    m_KnotList->ComputeBasisJet(grid[i].iKnot,grid[i].u,grid[i].basis);
    }
}

template <unsigned int VDimension, class TReal, unsigned int VOrder, unsigned int VJetOrder>
void
BSplineCurve<VDimension,TReal,VOrder,VJetOrder>
::CreateUniformEvaluationGrid(unsigned int nPoints, EvaluationGrid &grid)
{
  double *values = new double[nPoints];
  double step = 1.0 / (nPoints - 1);
  double val = 0.0;
  for (unsigned int i=0;i<nPoints;i++)
    {
    values[i] = val < 1.0 ? val : 1.0;
    val += step;
    }
  CreateEvaluationGrid(nPoints,values,grid);
  delete values;
}

