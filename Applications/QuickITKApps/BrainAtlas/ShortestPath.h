#ifndef __ShortestPath_h_
#define __ShortestPath_h_

#include "BinaryHeap.h"
#include <limits>

/**
 * This class implements the classic shortest path algorithm by the
 * legendary Dijkstra. It uses the binary heap implementation of
 * the priority queue that is more flexible than the implementation 
 * in STL
 *
 * The graph must be represented in the format that I lifted from the 
 * METIS program. It's very simple. Let's say we have a graph G(V,E). 
 * You must supply two arrays of integers: the first called the 
 * adjacency index array (AI) of length |V|+1 and the second called the 
 * adjacency (A) array of length |E|. The elements of AI are indices into
 * A, and elements of A represent adjacencies. 
 *
 * The element i in the array AI points to a position in array A. Starting
 * at that position are listed the vertices adjacent to i. The last element
 * in AI points past the end of A, i.e., AI[ |V| ] = |E|.
 *
 * To determine the vertices adjacent to vertex i, one looks at the values
 * A[AI[i]], ... ,A[AI[i+1] - 1]. 
 */
template <class TWeight>
class DijkstraShortestPath 
{
public:
  /** Constant representing infinite distance to soruce vertex, ie, no
   * path to source */
  static const TWeight INFINITE_WEIGHT;

  /** Constant representing no path to source vertex */
  static const unsigned int NO_PATH;

  /** 
   * Constructor, creates shortest path computer for given graph
   * specified in METIS format. In addition pass in the weights
   * associated with the edges. */
  DijkstraShortestPath(
    unsigned int nVertices, unsigned int *xAdjacencyIndex,
    unsigned int *xAdjacency, TWeight *xEdgeLen)
    {
    // Store the sizes
    m_NumberOfVertices = nVertices;
    m_NumberOfEdges = xAdjacencyIndex[nVertices];

    // Store the adjacency data
    m_Adjacency = xAdjacency;
    m_AdjacencyIndex = xAdjacencyIndex;
    m_EdgeWeight = xEdgeLen;

    // Allocate the array of distances
    m_Distance = new TWeight[nVertices];
    m_Predecessor = new unsigned int[nVertices];

    // Create the Binary heap (priority que)
    m_Heap = new BinaryHeap<TWeight>(m_NumberOfVertices, m_Distance);
    }

  /** Destructor, cleans up pointers */
  ~DijkstraShortestPath()
    {
    delete m_Heap;
    delete m_Distance;
    delete m_Predecessor;
    }

  /** Compute the shortest paths for given source vertex */
  void ComputePathsFromSource(unsigned int iSource)
    {
    unsigned int i;

    // Initialize the predecessor array (necessary?)
    for(i = 0; i < m_NumberOfVertices; i++)
      { m_Predecessor[i] = NO_PATH; } 

    // Reset the binary heap and weights
    m_Heap->InsertAllElementsWithEqualWeights(INFINITE_WEIGHT);

    // Change the distance for the first weight to 0
    m_Heap->DecreaseElementWeight(iSource, 0);

    // Set the predecessor of iSource to itself
    m_Predecessor[iSource] = iSource;

    // Change the distance to the adjacent vertices of the source
    for(i = m_AdjacencyIndex[iSource]; i < m_AdjacencyIndex[iSource+1]; i++)
      {
      // Get the neighbor of i
      unsigned int iNbr = m_Adjacency[i];

      // Get the edge weight associated with it and update it in the queue
      m_Heap->DecreaseElementWeight(iNbr, m_EdgeWeight[i]);

      // Update the predecessor as well 
      m_Predecessor[iNbr] = iSource;
      }

    // Continue while the heap is not empty
    while(m_Heap->GetSize())
      {
      // Pop off the closest vertex
      unsigned int w = m_Heap->PopMinimum();

      // Relax the vertices remaining in the heap
      for(i = m_AdjacencyIndex[w]; i < m_AdjacencyIndex[w+1]; i++)
        {
        // Get the neighbor of i
        unsigned int iNbr = m_Adjacency[i];

        // Get the edge weight associated with it and update it in the queue
        if(m_Heap->ContainsElement(iNbr))
          {
          // If the distance to iNbr more than distance thru w, update it
          TWeight dTest = m_Distance[w] + m_EdgeWeight[i];

          if(dTest < m_Distance[iNbr])
            {
            // Update the weight
            m_Heap->DecreaseElementWeight(iNbr, dTest);

            // Set the predecessor
            m_Predecessor[iNbr] = w;
            }
          }
        }
      } // while heap not empty
    }

  /** Get the predecessor array */
  const unsigned int *GetPredecessorArray()
    { return m_Predecessor; }

  /** Get the distance array */
  const TWeight * GetDistanceArray()
    { return m_Distance; }

private:
  BinaryHeap<TWeight> *m_Heap;
  TWeight *m_Distance;
  TWeight *m_EdgeWeight;
  unsigned int *m_Predecessor;
  unsigned int *m_AdjacencyIndex, *m_Adjacency;
  unsigned int m_NumberOfVertices, m_NumberOfEdges;
};

template<class TWeight>
const TWeight
DijkstraShortestPath<TWeight>
::INFINITE_WEIGHT = numeric_limits<TWeight>::max();

template<class TWeight>
const unsigned int
DijkstraShortestPath<TWeight>
::NO_PATH = numeric_limits<unsigned int>::max();

#endif
