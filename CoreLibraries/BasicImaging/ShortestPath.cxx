#include "ShortestPath.h"

using namespace std;

void TestShortestPath()
{
  // Create a graph (taken from website below)
  // http://ciips.ee.uwa.edu.au/~morris/Year2/PLDS210/dij-op.html
  size_t AI[] = {0,2,4,5,7,10};
  size_t A[] =  {1,  4,  2,  4,  3,  2,  0,  1,  2,  3};
  size_t w[] =  {10, 5,  1,  2,  4,  6,  7,  3,  9,  2};
  
  // Create shortest pather
  DijkstraShortestPath<size_t> sp(5, AI, A, w);

  // Compute shortest path from starting point
  sp.ComputePathsFromSource(0);

  // Reconstruct the shortest path to each vertex
  for(size_t i=1; i<5; i++)
    {
    cout << "Path to 0 from " << i << " : ";
    size_t j = i;
    while(j != 0)
      {
      cout << j << " - ";
      j = sp.GetPredecessorArray()[j];
      }
    cout << "0" << endl;
    }
}

/* 
int main(int argc, char *argv[])
{
  TestShortestPath();
}
*/ 


