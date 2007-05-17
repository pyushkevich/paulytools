#include "ImageToGraphFilter.h"
#include "ShortestPath.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include <iostream>
#include <algorithm>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matlab_filewrite.h>

using namespace std;

int usage()
{
  cerr << "Program: Extract random landmarks from binary image" << endl;
  cerr << "Usage: imglandgraph image.img xyz.mat dist.mat n" << endl;
  cerr << "Parameters: " << endl;
  cerr << "  image.img  Input binary image" << endl;
  cerr << "  xyz.mat    Output, coordinates of selected landmarks" << endl;
  cerr << "  dist.mat   Squared geodesic distance matrix between landmarks" << endl;
  cerr << "  n          Number of landmarks to generate" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 5)
    return usage();

  // Load the image
  typedef itk::Image<short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[1]);
  try {
    fltReader->Update();
  } catch(itk::ExceptionObject &exc) {
    cerr << "Unable to read image" << endl;
    cerr << exc.what() << endl;
  }
  ImageType::Pointer img = fltReader->GetOutput();

  // Convert the image to a graph
  typedef ImageToGraphFilter<ImageType, size_t, double> GraphFilter;
  GraphFilter::Pointer fltGraph = GraphFilter::New();
  fltGraph->SetInput(img);
  fltGraph->Update();

  // Get the number of landmarks
  size_t i, j;
  size_t n = atoi(argv[4]);

  // Generate a random subset of landmarks in the image
  size_t *xLandmarks = new size_t[fltGraph->GetNumberOfVertices()];
  for(i = 0; i < fltGraph->GetNumberOfVertices(); i++)
    xLandmarks[i] = i;
  random_shuffle(xLandmarks, xLandmarks+fltGraph->GetNumberOfVertices());

  // Compute distances between landmarks
  typedef DijkstraShortestPath<double> Dijkstra;
  Dijkstra dijk(fltGraph->GetNumberOfVertices(),
    fltGraph->GetAdjacencyIndex(),
    fltGraph->GetAdjacency(),
    fltGraph->GetEdgeWeights());

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
    }

  // Save the distance matrix
  vnl_matlab_filewrite exporter(argv[3]);
  exporter.write(mDist, "dist");

  // Save the coordinates of the landmarks
  vnl_matrix<double> mCoord(n, 3, 0.0);
  for(i = 0; i < n; i++)
    {
    ImageType::IndexType idx = fltGraph->GetVertexImageIndex(xLandmarks[i]);
    mCoord[i][0] = idx[0];
    mCoord[i][1] = idx[1];
    mCoord[i][2] = idx[2];
    }

  // Save the distance matrix
  vnl_matlab_filewrite exporter2(argv[2]);
  exporter2.write(mCoord, "xyz");

  delete xLandmarks;
}
