#include <vector>
#include <queue>

using namespace std;

template <class TVertex, class TDistance>
class DijkstraShortestPath 
{
public:
  /** Constructor, creates shortest path computer for given graph
   * specified in METIS form */
  DijkstraShortestPath(unsigned int nVertices, 
    TVertex *xAdjacencyIndex,TVertex *xAdjacency, TDistance *xEdgeLen);

  /** Compute the shortest paths for given source vertex */
  void ComputePathsFromSource(TVertex iSource);

private:
  vector<unsigned int> m_QueBase;    
  vector<TDistance> m_Distance;
  vector<TVertex> m_Predecessor;  
};

template<class TVertex, class TDistance>
DijkstraShortestPath<TVertex, TDistance>
::DijkstraShortestPath(unsigned int nVertices, 
    TVertex *xAdjacencyIndex,TVertex *xAdjacency, TDistance *xEdgeLen)
: m_QueBase(nVertices), m_Distance(nVertices), m_Predecessor(nVertices)
{
  
}

  