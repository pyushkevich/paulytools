#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>

#include <vtkCellArray.h>
#include <vtkCellDataToPointData.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include "ReadWriteVTK.h"
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkCleanPolyData.h>

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
 * Compute the mesh distances (#edges) between generator points in a mesh. Only
 * distances smaller than a threshold are considered
 */
void ComputeMeshDistances(
  vtkPolyData *mesh,
  VertexPairArray &pairs,
  size_t dMax,
  vector<double> &d)
{
  size_t i, n = mesh->GetNumberOfPoints();
  cout << "Computing Simple Geodesics on the Boundary" << endl;

  // Fill the output vector with zero
  d.resize(pairs.size(), 0.0);

  // Wrap the mesh in a half-edge structure
  VTKMeshHalfEdgeWrapper hewrap(mesh);

  // Create a Dijkstra object
  VTKMeshShortestDistance dijkstra;
  dijkstra.SetInputMesh(&hewrap);

  // Set the weight function to euclidean distance
  UnitLengthMeshEdgeWeightFunction wfunc;
  dijkstra.SetEdgeWeightFunction(&wfunc);

  // Compute the edge graph
  dijkstra.ComputeGraph();

  // For each vertex, compute the points around it
  for(i = 0; i < n; i++)
    {
    // Compute distances in a geodesic ball of radius rmax * dscale
    dijkstra.ComputeDistances(i, dMax);

    // Compute distances for all pairs that involve i (brute force again)
    size_t j = 0;
    for(VertexPairArray::iterator it = pairs.begin(); it != pairs.end(); ++it, ++j)
      {
      if(it->first == i)
        {
        d[j] = dijkstra.GetVertexDistance(it->second);
        }
      else if(it->second == i)
        {
        d[j] = dijkstra.GetVertexDistance(it->first);
        }
      }

    // Track progress
    if((i % 1000) == 999) cout << "." << flush;
    }
  cout << endl;
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
  cout << "Computing Geodesics on the Boundary" << endl;

  // For each vertex, compute its radius (right now these are only for pairs)
  vector<double> rmax(n, 0.0);
  for(i = 0; i < pairs.size(); i++)
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
        {
        d[j] = dijkstra.GetVertexDistance(it->second);
        }
      else if(it->second == i)
        {
        d[j] = dijkstra.GetVertexDistance(it->first);
        }
      }

    // Track progress
    if((i % 1000) == 999) cout << "." << flush;
    }
  cout << endl;
}

vtkPolyData *ReadVoronoiOutput(
  string fn,                  // Filename from which to read the points
  ImageType *mask,            // The mask image (to tell if diagram is inside the object)
  float threshold,            // The threshold for the mask image
  VertexPairArray &src)       // Output: for each cell in the VD, the pair of generators
{
  // Build an interpolator for the image
  typedef itk::LinearInterpolateImageFunction<ImageType, double> FuncType;
  FuncType::Pointer fMask = NULL;
  if(mask)
    {
    fMask = FuncType::New();
    fMask->SetInputImage(mask);
    }

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

  // Clear the src array
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
      if(!fMask.IsNull())
        {
        FuncType::PointType P;
        P[0] = pts->GetPoint(ids[k])[0];
        P[1] = pts->GetPoint(ids[k])[1];
        P[2] = pts->GetPoint(ids[k])[2];

        if(!fMask->IsInsideBuffer(P) || fMask->Evaluate(P) < threshold)
          isout = true;
        }
      }

    if(!isinf && !isout)
      {
      // add the cell
      cells->InsertNextCell(m, ids);

      // add the pair of generators
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

int usage()
{
  cout << "Usage: " << endl;
  cout << "    importqvoronoi [options] voronoi.txt mesh.vtk output.vtk" << endl;
  cout << "Parameters: " << endl;
  cout << "    voronoi.txt          Output produced by \"qvoronoi p Fv\"" << endl;
  cout << "    mesh.vtk             Mesh from which the voronoi diag was made" << endl;
  cout << "    output.vtk           Mesh where to save the skeleton" << endl;
  cout << "Options:" << endl;
  cout << "    -e N                 Minimal number of mesh edges separating two generator" << endl;
  cout << "                         points of a VD face for it to be considered (try 2, 3)" << endl;
  cout << "    -p X.XX              Prune the mesh using factor X.XX (try 2.0). The " << endl;
  cout << "                         pruning algorithm deletes faces in the VD for " << endl;
  cout << "                         which the ratio of the geodesic distance between " << endl;
  cout << "                         the generating points and the euclidean distance " << endl;
  cout << "                         between these points is less than X.XX" << endl;
  cout << "    -m mask.img X.XX     Use use-supplied binary image as a mask. Faces in the VD" << endl;
  cout << "                         that have a vertex where the mask is below the threshhold" << endl;
  cout << "                         X.XX will be eliminated" << endl;
  cout << "    -c N                 Take at most N connected components of the skeleton" << endl;
  cout << "    -s mesh.vtk          Load a skeleton from mesh.vtk and compare to the output skeleton" << endl;
  cout << "    -g                   Compute full geodesic information. This is only useful for" << endl;
  cout << "                         debugging the pruning code." << endl;
  return -1;
}
  

int main(int argc, char *argv[])
{
  size_t k;

  // Command line arguments
  string fnVoronoi, fnMesh, fnOutput, fnMask, fnSkel;
  double xPrune = 2.0, xThresh = 0.0;
  int nComp = 0, nDegrees = 0;
  bool flagGeodFull = false;

  // Check that there are at least three command line arguments
  if(argc < 4) return usage();
  fnVoronoi = argv[argc-3];
  fnMesh = argv[argc-2];
  fnOutput = argv[argc-1];
  
  // Parse the command line for options
  for(size_t iArg = 1; iArg < argc - 3; ++iArg)
    {
    string arg = argv[iArg];
    if(arg == "-p")
      {
      xPrune = atof(argv[++iArg]);
      }
    else if(arg == "-m")
      {
      fnMask = argv[++iArg];
      xThresh = atof(argv[++iArg]);
      }
    else if(arg == "-c")
      {
      nComp = atoi(argv[++iArg]);
      }
    else if(arg == "-e")
      {
      nDegrees = atoi(argv[++iArg]);
      }
    else if(arg == "-s")
      {
      fnSkel = argv[++iArg];
      }
    else if(arg == "-g")
      {
      flagGeodFull = true;
      }
    else return usage();
    }

  // Load the mask image if specified
  ImageType::Pointer imgMask = NULL;
  if(fnMask.length())
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fnMask.c_str());
    try 
      {
      reader->Update();
      imgMask = reader->GetOutput();
      } 
    catch (itk::ExceptionObject &exc) 
      {
      cerr << exc << endl;
      return -1;
      }
    }

  // Load the input mesh
  vtkPolyData *bnd = ReadVTKData(fnMesh);
  bnd->BuildCells();

  // Run the program 
  VertexPairArray pgen;
  vtkPolyData *poly = ReadVoronoiOutput(fnVoronoi, imgMask, xThresh, pgen);
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

  // Compute the number of mesh edges between generator points
  std::vector<double> meshdist(pgen.size());
  if(nDegrees)
    ComputeMeshDistances(bnd, pgen, nDegrees, meshdist);

  // Compute the geodesic distances for the pairs of associated boundary points
  std::vector<double> geod(pgen.size());
  ComputeGeodesics(bnd, pgen, rad, flagGeodFull ? 1e100 : xThresh, geod);

  // Create a new cell array for the pruning
  vtkCellArray *prune = vtkCellArray::New();

  // Allocate the output radius data array
  vtkDoubleArray *daRad = vtkDoubleArray::New();
  daRad->SetNumberOfComponents(1);
  daRad->SetName("Radius");

  // Another array for prune strength
  vtkDoubleArray *daPrune = vtkDoubleArray::New();
  daPrune->SetNumberOfComponents(1);
  daPrune->SetName("Pruning Ratio");

  // Another array for prune strength
  vtkDoubleArray *daGeod = vtkDoubleArray::New();
  daGeod->SetNumberOfComponents(1);
  daGeod->SetName("Geodesic");

  // Prune the voronoi diagram
  for(k = 0; k < pgen.size(); k++)
    {
    double xRad = rad[k];
    double xGeod = geod[k];
    double xRatio = xGeod / xRad;
    double xMeshDist = nDegrees ? meshdist[k] : 0.0;
    if(xMeshDist >= nDegrees && xRatio > xPrune)
      {
      prune->InsertNextCell(poly->GetCell(k));
      daRad->InsertNextTuple(&xRad);
      daGeod->InsertNextTuple(&xGeod);
      daPrune->InsertNextTuple(&xRatio);
      }
    }
     
  // Insert the new pruned cell list
  poly->SetPolys(prune);
  poly->GetCellData()->AddArray(daRad);
  poly->GetCellData()->AddArray(daPrune);
  poly->GetCellData()->AddArray(daGeod);
  poly->BuildCells();
  poly->BuildLinks();

  // Compute the distance from the medial axis to the voronoi mesh
  if(fnSkel != "")
    {
    // Load the input medial mesh
    vtkPolyData *pd = ReadVTKData(fnSkel);
    vtkTriangleFilter *tri = vtkTriangleFilter::New();
    tri->SetInput(pd);
    tri->Update();
    vtkPolyData *med = tri->GetOutput();
    med->BuildCells();
    med->BuildLinks();

    vnl_vector<vtkFloatingPointType> dist;
    vnl_vector<vtkFloatingPointType> area;
    ComputeAreaElement(med, area);
    ComputeExactMeshToMeshDistance(med, poly, dist);
    cout << "Avg, Max, Edge: " << dot_product(dist, area) / area.one_norm() 
      << " " << dist.inf_norm() << " " << ComputeAverageEdgeLength(med) << endl;
    }

  // Drop the singleton points from the diagram (points not in any cell)
  vtkCleanPolyData *fltClean = vtkCleanPolyData::New();
  fltClean->SetInput(poly);
  fltClean->Update();

  // The output from the next branch
  vtkPolyData *polySave = fltClean->GetOutput();

  // Compute the connected components
  if(nComp > 0)
    {
    vtkPolyDataConnectivityFilter *fltConnect = vtkPolyDataConnectivityFilter::New();
    fltConnect->SetInput(fltClean->GetOutput());

    if(nComp == 1)
      fltConnect->SetExtractionModeToLargestRegion();
    else 
      {
      fltConnect->SetExractionModeToSelectedRegions();
      fltConnect->InitializeSelectedRegionList();
      for(size_t rr = 0; rr < nComp; rr++)
        fltConnect->AddSelectedRegion(rr);
      }

    fltConnect->ScalarConnectivityOff();
    fltConnect->Update();
    polySave = fltConnect->GetOutput();
    }

  // Convert the cell data to point data
  vtkCellDataToPointData *c2p = vtkCellDataToPointData::New();
  c2p->SetInput(polySave);
  c2p->PassCellDataOn();
  c2p->Update();
    
  WriteVTKData(c2p->GetPolyDataOutput(), fnOutput);
  return 0;
}


