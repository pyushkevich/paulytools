#ifndef __ShortestPath_h_
#define __ShortestPath_h_

#include "BinaryHeap.h"
#include <limits>
#include <iostream>

using namespace std;

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
template <typename TWeight>
class DijkstraShortestPath 
{
public:
  /** Constant representing infinite distance to soruce vertex, ie, no
   * path to source */
  static const TWeight INFINITE_WEIGHT;

  /** Constant representing no path to source vertex */
  static const size_t NO_PATH;

  /** 
   * Constructor, creates shortest path computer for given graph
   * specified in METIS format. In addition pass in the weights
   * associated with the edges. */
  DijkstraShortestPath(
    size_t nVertices, size_t *xAdjacencyIndex,
    size_t *xAdjacency, TWeight *xEdgeLen)
    {
    // Store the sizes
    this->m_NumberOfVertices = nVertices;
    this->m_NumberOfEdges = xAdjacencyIndex[nVertices];

    // Store the adjacency data
    m_Adjacency = xAdjacency;
    m_AdjacencyIndex = xAdjacencyIndex;
    m_EdgeWeight = xEdgeLen;

    // Allocate the array of distances
    m_Distance = new TWeight[nVertices];
    m_Predecessor = new size_t[nVertices];

    // Create the Binary heap (priority que)
    m_Heap = new BinaryHeap<TWeight>(m_NumberOfVertices, m_Distance);
    }

  /** Destructor, cleans up pointers */
  virtual ~DijkstraShortestPath()
    {
    delete m_Heap;
    delete m_Distance;
    delete m_Predecessor;
    }

  /** Compute the shortest paths for given source vertex */
  void ComputePathsFromSource(size_t iSource, size_t iTarget = 0xffffffff)
    {
    size_t i;

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
      size_t iNbr = m_Adjacency[i];

      // Get the edge weight associated with it and update it in the queue
      m_Heap->DecreaseElementWeight(iNbr, m_EdgeWeight[i]);

      // Update the predecessor as well 
      m_Predecessor[iNbr] = iSource;
      }

    // Continue while the heap is not empty
    while(m_Heap->GetSize())
      {
      // Pop off the closest vertex
      size_t w = m_Heap->PopMinimum();

      // If w is the target vertex, we are done
      if(w == iTarget) break;

      // Relax the vertices remaining in the heap
      for(i = m_AdjacencyIndex[w]; i < m_AdjacencyIndex[w+1]; i++)
        {
        // Get the neighbor of i
        size_t iNbr = m_Adjacency[i];

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
  const size_t *GetPredecessorArray()
    { return m_Predecessor; }

  /** Get the distance array */
  const TWeight * GetDistanceArray()
    { return m_Distance; }

protected:
  BinaryHeap<TWeight> *m_Heap;
  TWeight *m_Distance;
  TWeight *m_EdgeWeight;
  size_t *m_Predecessor;
  size_t *m_AdjacencyIndex, *m_Adjacency;
  size_t m_NumberOfVertices, m_NumberOfEdges;
};

template<class TWeight>
class GraphVoronoiDiagram : public DijkstraShortestPath<TWeight>
{
public:
  typedef DijkstraShortestPath<TWeight> Superclass;
  
  GraphVoronoiDiagram(
    size_t nVertices, size_t *xAdjacencyIndex,
    size_t *xAdjacency, TWeight *xEdgeLen) : 
    Superclass(nVertices, xAdjacencyIndex, xAdjacency, xEdgeLen)
    {
      m_Source = new size_t[nVertices];
      cout << "GVD : " << nVertices << " " << this->m_NumberOfEdges << endl;
    }

  virtual ~GraphVoronoiDiagram()
    { delete m_Source; }

  /** Compute paths from multiple sources. Use this method to construct a
   * sort of a Voronoi diagram of the graph */ 
  void ComputePathsFromManySources(size_t nSources, size_t *lSources)
    {
    size_t i;

    cout << "Number of sources: " << nSources << endl;
    cout << "Number of vertices: " << this->m_NumberOfVertices << endl;
    cout << "Number of edges: " << this->m_NumberOfEdges << endl;

    // Initialize the predecessor array and the source array 
    for(i = 0; i < this->m_NumberOfVertices; i++)
      { 
      this->m_Predecessor[i] = DijkstraShortestPath<TWeight>::NO_PATH; 
      this->m_Source[i] = DijkstraShortestPath<TWeight>::NO_PATH;
      } 

    // Reset the binary heap and weights
    this->m_Heap->InsertAllElementsWithEqualWeights(this->INFINITE_WEIGHT);

    // Initialize the heap with the sources
    for(size_t iSource = 0; iSource < nSources; iSource++)
      {
      // Get the source vertex
      size_t id = lSources[iSource];

      // Change the weight for all the sources to 0
      this->m_Heap->DecreaseElementWeight(id, 0);

      // Set the predecessor of iSource to itself
      m_Source[id] = this->m_Predecessor[id] = id;
      }

    // Continue while the heap is not empty
    while(this->m_Heap->GetSize())
      {
      // Pop off the closest vertex
      size_t w = this->m_Heap->PopMinimum();

      // Relax the vertices remaining in the heap
      for(i = this->m_AdjacencyIndex[w]; i < this->m_AdjacencyIndex[w+1]; i++)
        {
        // Get the neighbor of i
        size_t iNbr = this->m_Adjacency[i];

        // Get the edge weight associated with it and update it in the queue
        if(this->m_Heap->ContainsElement(iNbr))
          {
          // If the distance to iNbr more than distance thru w, update it
          TWeight dTest = this->m_Distance[w] + this->m_EdgeWeight[i];

          if(dTest < this->m_Distance[iNbr])
            {
            // Update the weight
            this->m_Heap->DecreaseElementWeight(iNbr, dTest);

            // Set the predecessor
            this->m_Predecessor[iNbr] = w;

            // Set the source of the vertex
            this->m_Source[iNbr] = this->m_Source[w];
            }
          }
        }
      } // while heap not empty
    }

  /** Get the source of a vertex, i.e., the source vertex that is closest to
    the given vertex, or NO_PATH if there is no path to the given vertex from
    any of the sources */
  size_t GetVertexSource(size_t iVertex) const
    { return m_Source[iVertex]; }

private:
  size_t *m_Source;
};

template<class TWeight>
const TWeight
DijkstraShortestPath<TWeight>
::INFINITE_WEIGHT = std::numeric_limits<TWeight>::max();

template<class TWeight>
const size_t
DijkstraShortestPath<TWeight>
::NO_PATH = std::numeric_limits<size_t>::max();

#endif
