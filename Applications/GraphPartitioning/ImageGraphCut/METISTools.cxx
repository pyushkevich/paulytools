#include "METISTools.h"
#include <iostream>

using namespace std;

void
MetisPartitionProblem
::SetProblem(GraphFilter *fltGraph, unsigned int nParts)
{
  m_Graph = fltGraph;
  m_Partition = new idxtype[m_Graph->GetNumberOfVertices()];
  m_NumberOfParameters = nParts;
}

MetisPartitionProblem::MeasureType
MetisPartitionProblem
::GetValue(const ParametersType &x) const 
{
  // Convert the vector to floating array
  vnl_vector<float> wgt(x.size() + 1);
  float sum = 0;
  for(unsigned int i=0;i<x.size();i++)
    {
    wgt[i] = (float) x[i];
    sum += wgt[i];
    }
  wgt[x.size()] = 1.0f - sum;

  // Run the METIS code
  cout << " Running METIS iteration [ x = " << x << "] " << endl;
  int edgecut = 
    RunMETISPartition(m_Graph, wgt.size(), wgt.data_block(), m_Partition); 

  // Report the edge cut
  cout << "  - done - edge cut " << edgecut << endl;

  return (MeasureType) edgecut;
}