#ifndef __METISTools_h_
#define __METISTools_h_

#include "ImageToGraphFilter.h"
#include <itkSingleValuedCostFunction.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkNormalVariateGenerator.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_powell.h>

using namespace itk;

typedef int idxtype;
extern "C" {
  void METIS_WPartGraphKway(
  	int *, idxtype *, idxtype *, idxtype *, idxtype *, 
  	int *, int *, int *, float *, int *, int *, idxtype *);

  void METIS_WPartGraphRecursive(
  	int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, 
    int *, int *, float *, int *, int *, idxtype *);
}

/** Function to run METIS using ImageToGraphFilter */
template< class TImage >
int RunMETISPartition(
  ImageToGraphFilter<TImage> *fltGraph,
  int nParts,
  float *xPartWeights,
  int *outPartition,
  bool useRecursiveAlgorithm = true)
{
  // Variables used to call METIS
  int nVertices = fltGraph->GetNumberOfVertices();
  int wgtflag = 3;
  int numflag = 0;
  int options[] = {0,0,0,0,0};
  int edgecut = 0;

  if( useRecursiveAlgorithm )
    {
    METIS_WPartGraphRecursive(
      &nVertices,
      fltGraph->GetAdjacencyIndex(),
      fltGraph->GetAdjacency(),
      fltGraph->GetVertexWeights(),
      fltGraph->GetEdgeWeights(),
      &wgtflag,
      &numflag,
      &nParts,
      xPartWeights,
      options,
      &edgecut,
      outPartition);
    }
  else
    {
    METIS_WPartGraphKway(
      &nVertices,
      fltGraph->GetAdjacencyIndex(),
      fltGraph->GetAdjacency(),
      fltGraph->GetVertexWeights(),
      fltGraph->GetEdgeWeights(),
      &wgtflag,
      &numflag,
      &nParts,
      xPartWeights,
      options,
      &edgecut,
      outPartition);
    }  

  return edgecut;
}

/*
 * METIS PARTITION OPTIMIZATION
 *
 * This experimental code is used to optimize the METIS result over the
 * relative weights of the partitions. This optimization calls METIS in
 * the inner loop and therefore can be a little slow.
 */

class MetisPartitionProblem : public SingleValuedCostFunction
{
public:
  typedef MetisPartitionProblem Self;
  typedef SingleValuedCostFunction Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  itkTypeMacro(MetisPartitionProblem, SingleValuedCostFunction);

  itkNewMacro(Self);
  
  typedef itk::Image<short,3> ImageType; 
  typedef ImageToGraphFilter<ImageType> GraphFilter;

  typedef Superclass::MeasureType MeasureType;
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::DerivativeType DerivativeType;

  /** Set the problem parameters */
  void SetProblem(GraphFilter *fltGraph, unsigned int nParts);

  /** Return the number of parameters */
  unsigned int GetNumberOfParameters() const
    { return m_NumberOfParameters; }
  
  /** Virtual method from parent class */
  MeasureType GetValue(const ParametersType &x) const;

  /** Not used, since there are no derivatives to evaluate */
  void GetDerivative(const ParametersType &, DerivativeType &) const {}

  /** Get the result of the last partition */
  idxtype *GetLastPartition() { return m_Partition; }
  
private:
  /** The stored graph information */
  GraphFilter *m_Graph;

  /** The partition array */
  idxtype *m_Partition;

  /** Problem size */
  unsigned int m_NumberOfParameters;
};


#endif // __METISTools_h_
