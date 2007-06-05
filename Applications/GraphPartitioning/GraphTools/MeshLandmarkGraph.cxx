#include "ShortestPath.h"
#include "ReadWriteVTK.h"
#include "VTKMeshShortestDistance.h"
#include <vtkCleanPolyData.h>
#include "SparseMatrix.h"
#include <iostream>
#include <algorithm>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matlab_filewrite.h>

using namespace std;

int usage()
{
  cerr << "Program: Extract random landmarks from a VTK mesh" << endl;
  cerr << "Usage: meshlandgraph mesh.vtk xyz.mat dist.mat n" << endl;
  cerr << "Parameters: " << endl;
  cerr << "  mesh.vtk  Input binary image" << endl;
  cerr << "  xyz.mat    Output, coordinates of selected landmarks" << endl;
  cerr << "  dist.mat   Squared geodesic distance matrix between landmarks" << endl;
  cerr << "  n          Number of landmarks to generate" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  size_t i, j;
  if(argc < 5)
    return usage();

  // Load the mesh
  vtkPolyData *xMeshRaw = ReadVTKData(argv[1]);

  // Convert the mesh to triangles
  vtkTriangleFilter *fltTri = vtkTriangleFilter::New();
  fltTri->SetInput(xMeshRaw);

  // Clean the mesh
  vtkCleanPolyData *fltClean = vtkCleanPolyData::New();
  fltClean->SetInput(fltTri->GetOutput());
  fltClean->Update();
  vtkPolyData *xMesh = fltClean->GetOutput();
  xMesh->BuildCells();
  xMesh->BuildLinks();
  
  // Compute all the edges
  typedef ImmutableSparseArray<double> ImmutableArr;
  typedef ImmutableArr::VNLSourceType MutableArr;

  // Initialize the adjacency matrix
  size_t np = xMesh->GetNumberOfPoints();
  MutableArr M0(np, np);
  
  // Traverse all the cells in the VTK mesh, recording all the available
  // edges in the mesh. 
  for(unsigned int iCell = 0; iCell < xMesh->GetNumberOfCells(); iCell++)
    {
    // Get the points for this cell
    vtkIdType nPoints, *xPoints;
    xMesh->GetCellPoints(iCell, nPoints, xPoints);

    // Walk around the list of points
    for(unsigned int j = 0; j < nPoints; j++)
      {
      // Get the head and the tail of the current half-edge
      unsigned int iTail = xPoints[j], iHead = xPoints[(j+1) % nPoints];

      // Compute the distance
      vnl_vector<double> x1(xMesh->GetPoint(iTail), 3);
      vnl_vector<double> x2(xMesh->GetPoint(iHead), 3);
      double dist = (x1 - x2).magnitude();

      // Add the edge to the list
      M0(iTail,iHead) = dist;
      M0(iHead,iTail) = dist;
      }
    }

  // Build a sparse matrix from the edge data
  ImmutableArr M;
  M.SetFromVNL(M0); 
  
  // Generate a random subset of landmarks in the image
  size_t n = atoi(argv[4]);
  size_t *xLandmarks = new size_t[np];
  for(i = 0; i < np; i++) xLandmarks[i] = i;
  random_shuffle(xLandmarks, xLandmarks+np);

  // Compute distances between landmarks
  typedef DijkstraShortestPath<double> Dijkstra;
  Dijkstra dijk(np, M.GetRowIndex(), M.GetColIndex(), M.GetSparseData());

  // Initialize the distance matrix
  vnl_matrix<double> mDist(n, n, 0.0);
  
  // Using the brute force approach
  for(i = 0; i < n; i++)
    {
    // Compute all pairs shortest paths
    dijk.ComputePathsFromSource(xLandmarks[i]);
    const double *xDist = dijk.GetDistanceArray();

    // Get the landmark-to-landmark distances
    for(j = 0; j < n; j++)
      {
      double d = xDist[xLandmarks[j]];
      mDist[i][j] = d * d;
      }

    cout << ".";
    if( ((i+1) % 64) == 0 || i + 1 == n )
      cout << " n = " << i+1 << endl;
    else
      cout << flush;
    }

  // Save the distance matrix
  vnl_matlab_filewrite exporter(argv[3]);
  exporter.write(mDist, "dist");

  // Save the coordinates of the landmarks
  vnl_matrix<double> mCoord(n, 3, 0.0);
  for(i = 0; i < n; i++)
    {
    double *xPoint = xMesh->GetPoint(xLandmarks[i]);
    mCoord[i][0] = xPoint[0];
    mCoord[i][1] = xPoint[1];
    mCoord[i][2] = xPoint[2];
    }

  // Save the distance matrix
  vnl_matlab_filewrite exporter2(argv[2]);
  exporter2.write(mCoord, "xyz");

  delete xLandmarks;
}
