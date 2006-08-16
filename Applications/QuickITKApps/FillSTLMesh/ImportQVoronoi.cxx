#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCell.h>
#include <vtkCellLocator.h>

#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImage.h>

#include "VTKMeshShortestDistance.h"

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_cross.h>

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

using namespace std;

typedef itk::Image<float, 3> ImageType;

typedef std::pair<vtkIdType, vtkIdType> VertexPair;
typedef std::set< std::pair<vtkIdType, vtkIdType> > VertexPairSet;
typedef std::vector<VertexPair> VertexPairArray;

void ComputeExactMeshToMeshDistance(vtkPolyData *source, vtkPolyData *target, vnl_vector<vtkFloatingPointType> &dist)
{
  // Create a cell locator from the target mesh
  vtkCellLocator *xLocator = vtkCellLocator::New();
  xLocator->SetDataSet(target);
  xLocator->BuildLocator();

  // Expand the distance array to match the number of points in the source
  dist.set_size(source->GetNumberOfPoints());

  // Compute each distance
  for(size_t i = 0; i < dist.size(); i++)
    {
    vtkFloatingPointType xClosest[3];
    vtkIdType cellId;
    int subId;
    xLocator->FindClosestPoint(source->GetPoint(i), xClosest, cellId, subId, dist[i]);
    }

  // Take square root
  dist = dist.apply(sqrt);
}

/** 
 * Compute the distances between vertex pairs. The first parameter is the mesh,
 * the second is the list of vertex pairs, the third is the list of 'maximal' 
 * distances (after which the geodesic is not computed) and the last parameter 
 * is the output vector
 */
void ComputeGeodesics(
  vtkPolyData *mesh,                  // Boundary mesh (VTK)
  VertexPairArray &pairs,             // Pairs of points between which to compute geodesics
  vector<double> &dmax,               // Threshold for geodesic distance for each point pair
  double dmaxScale,                   // Constant scaling factor for above threshold
  vector<double> &d)                  // Output vector (geodesic distance)
{
  size_t i, n = mesh->GetNumberOfPoints();

  // For each vertex, compute its radius (right now these are only for pairs)
  vector<double> rmax(n, 0.0);
  for(i = 0; i < dmax.size(); i++)
    {
    rmax[pairs[i].first] = std::max(rmax[pairs[i].first], dmax[i]);
    rmax[pairs[i].second] = std::max(rmax[pairs[i].second], dmax[i]);
    }

  // Fill the output vector with zero
  d.resize(pairs.size(), 0.0);

  // We need to construct a data structure where for each vertex we have the list
  // of vertex pairs that pass through it. Right now we use brute force...
  
  // Wrap the mesh in a half-edge structure
  VTKMeshHalfEdgeWrapper hewrap(mesh);

  // Create a Dijkstra object
  VTKMeshShortestDistance dijkstra;
  dijkstra.SetInputMesh(&hewrap);

  // Set the weight function to euclidean distance
  EuclideanDistanceMeshEdgeWeightFunction wfunc;
  dijkstra.SetEdgeWeightFunction(&wfunc);

  // Compute the edge graph
  dijkstra.ComputeGraph();

  // For each vertex, compute the points around it
  for(i = 0; i < n; i++)
    {
    // Compute distances in a geodesic ball of radius rmax * dscale
    dijkstra.ComputeDistances(i, rmax[i] * dmaxScale);

    // Compute distances for all pairs that involve i (brute force again)
    size_t j = 0;
    for(VertexPairArray::iterator it = pairs.begin(); it != pairs.end(); ++it, ++j)
      {
      if(it->first == i)
        d[j] = dijkstra.GetVertexDistance(it->second);
      else if(it->second == i)
        d[j] = dijkstra.GetVertexDistance(it->first);
      }
    }
}

vtkPolyData *ReadVoronoiOutput(
  string fn,                // Filename from which to read the points
  ImageType *mask,          // The mask image (to tell if diagram is inside the object)
  float threshold,          // The threshold for the mask image
  VertexPairArray &src)       // Output: for each cell in the VD, the pair of generators
{
  // Build an interpolator for the image
  typedef itk::LinearInterpolateImageFunction<ImageType, double> FuncType;
  FuncType::Pointer fMask = FuncType::New();
  fMask->SetInputImage(mask);
  
  // Load the file
  ifstream fin(fn.c_str());

  // Load the numbers
  size_t nv, np, junk;
  
  // First two lines
  fin >> junk;
  fin >> nv; 

  // Create an array of points
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(nv);
  for(size_t i = 0; i < nv; i++)
    {
    double x,y,z;
    fin >> x;
    fin >> y;
    fin >> z;
    pts->SetPoint(i,x,y,z);
    }

  // Read the number of cells
  fin >> np;

  // Clear the src and r arrays
  src.clear();

  // Create the polygons 
  vtkCellArray *cells = vtkCellArray::New();
  for(size_t j = 0; j < np; j++)
    {
    bool isinf = false;
    bool isout = false;
    
    size_t m; fin >> m; m-=2;
    vtkIdType ip1, ip2; fin >> ip1; fin >> ip2;
    
    vtkIdType *ids = new vtkIdType[m];
    for(size_t k = 0; k < m; k++)
      {
      fin >> ids[k];

      // Is this point at infinity?
      if(ids[k] == 0) isinf = true; else ids[k]--;

      // Is this point in the image?
      FuncType::PointType P;
      P[0] = pts->GetPoint(ids[k])[0];
      P[1] = pts->GetPoint(ids[k])[1];
      P[2] = pts->GetPoint(ids[k])[2];
      if(!fMask->IsInsideBuffer(P) || fMask->Evaluate(P) > threshold)
        isout = true;
      }

    if(!isinf && !isout)
      {
      cells->InsertNextCell(m, ids);
      src.push_back(make_pair(ip1, ip2)); // TODO: is this numbering 0-based?
      }

    delete ids;
    }

  // Create the vtk poly data
  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(pts);
  poly->SetPolys(cells);
  return poly;
}

inline vtkFloatingPointType TriangleArea(
  const vnl_vector_fixed<vtkFloatingPointType,3> &A, 
  const vnl_vector_fixed<vtkFloatingPointType,3> &B, 
  const vnl_vector_fixed<vtkFloatingPointType,3> &C)
{
  return 0.5 * vnl_cross_3d(B - A, C - A).magnitude();
}

double ComputeAverageEdgeLength(vtkPolyData *poly)
{
  double l = 0.0;
  size_t n = 0;
  
  vtkIdType nCells = poly->GetNumberOfCells();
  for(vtkIdType iCell = 0; iCell < nCells; iCell++)
    {
    // Get the points in this cell
    vtkIdType nPoints, *xPoints;
    poly->GetCellPoints(iCell, nPoints, xPoints);

    for(vtkIdType j = 0; j < nPoints; j++)
      {
      vtkIdType k = (j + 1) % nPoints;
      vnl_vector_fixed<vtkFloatingPointType,3> x1(poly->GetPoint(xPoints[j]));
      vnl_vector_fixed<vtkFloatingPointType,3> x2(poly->GetPoint(xPoints[k]));
      l += sqrt(dot_product(x1-x2,x1-x2));
      n++;
      }
    }

  return l / n;
}

void ComputeAreaElement(vtkPolyData *poly, vnl_vector<vtkFloatingPointType> &elt)
{
  // For each triangle in the polydata compute its area
  vtkIdType nCells = poly->GetNumberOfCells();
  
  // Initialize the area elt array
  elt.set_size(poly->GetNumberOfPoints());
  elt.fill(0.0);
  
  for(vtkIdType iCell = 0; iCell < nCells; iCell++)
    {
    // Get the points in this cell
    vtkIdType nPoints, *xPoints;
    poly->GetCellPoints(iCell, nPoints, xPoints);
    
    // Only triangles are admitted
    if(nPoints != 3)
      { 
      cerr << "Irregular face (n = " << nPoints << ") detected!" << endl;
      break;
      }

    // Get the three points
    vnl_vector_fixed<vtkFloatingPointType, 3> X0(poly->GetPoint(xPoints[0]));
    vnl_vector_fixed<vtkFloatingPointType, 3> X1(poly->GetPoint(xPoints[1]));
    vnl_vector_fixed<vtkFloatingPointType, 3> X2(poly->GetPoint(xPoints[2]));

    // Compute the area
    double xArea = TriangleArea(X0, X1, X2);
    if(xArea < 0)
      {
      cerr << "Negative area returned at cell " << iCell << endl;
      break;
      }

    // Add the area to all points
    elt[xPoints[0]] += xArea; elt[xPoints[1]] += xArea; elt[xPoints[2]] += xArea;
    }
}

int main(int argc, char *argv[])
{
  size_t k;

  if(argc != 7)
    {
    cerr << "Usage: importqvoronoi voronoi.txt mask.img mask_thresh bnd.vtk med.vtk output.vtk" << endl;
    cerr << "Diagrams must be gerenated with \"qvoronoi p Fv\"" << endl;
    return -1;
    }

  // Load the image
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  try {
    reader->Update();
  } catch (itk::ExceptionObject &exc) {
    cerr << exc << endl;
    return -1;
  }

  // Load the input mesh
  vtkPolyDataReader *pdr = vtkPolyDataReader::New();
  pdr->SetFileName(argv[4]);
  pdr->Update();
  vtkPolyData *bnd = pdr->GetOutput();
  bnd->BuildCells();

  // Load the input medial mesh
  vtkPolyDataReader *pdr2 = vtkPolyDataReader::New();
  pdr2->SetFileName(argv[5]);
  pdr2->Update();
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInput(pdr2->GetOutput());
  tri->Update();
  vtkPolyData *med = tri->GetOutput();
  med->BuildCells();
  med->BuildLinks();

  // Run the program 
  VertexPairArray pgen;
  vtkPolyData *poly = ReadVoronoiOutput(
    argv[1], reader->GetOutput(), atof(argv[3]), pgen);
  poly->BuildCells();
  poly->BuildLinks();

  // Compute the distances between pairs of associated boundary points
  std::vector<double> rad(pgen.size());
  for(k = 0; k < pgen.size(); k++)
    {
    vnl_vector_fixed<vtkFloatingPointType, 3> p1(bnd->GetPoint(pgen[k].first));
    vnl_vector_fixed<vtkFloatingPointType, 3> p2(bnd->GetPoint(pgen[k].second));
    rad[k] = (p1 - p2).magnitude();
    }

  // Compute the geodesic distances for the pairs of associated boundary points
  // TODO: ignore geodesics longer than k * r
  std::vector<double> geod(pgen.size());
  ComputeGeodesics(bnd, pgen, rad, 4.0, geod);

  // Create a new cell array for the pruning
  vtkCellArray *prune = vtkCellArray::New();

  // Prune the voronoi diagram
  for(k = 0; k < pgen.size(); k++)
    if(geod[k] >= 4.0 * rad[k])
      prune->InsertNextCell(poly->GetCell(k));
     
  // Insert the new pruned cell list
  poly->SetPolys(prune);

  // Compute the distance from the medial axis to the voronoi mesh
  vnl_vector<vtkFloatingPointType> dist;
  vnl_vector<vtkFloatingPointType> area;
  ComputeAreaElement(med, area);
  ComputeExactMeshToMeshDistance(med, poly, dist);
  cout << "Avg, Max, Edge: " << dot_product(dist, area) / area.one_norm() 
    << " " << dist.inf_norm() << " " << ComputeAverageEdgeLength(med) << endl;
  
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(poly);
  writer->SetFileName(argv[argc-1]);
  writer->Update();
  return 0;
}


