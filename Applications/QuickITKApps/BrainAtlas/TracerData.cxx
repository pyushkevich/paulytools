#include "TracerData.h"
#include <fstream>

const int
TracerData
::NO_FOCUS = -1;

const vtkIdType
TracerData
::NO_MARKER = -1;

TracerData
::TracerData()
{
  // Initialize the mesh edge weight function
  m_EdgeWeightFunction = new EuclideanDistanceMeshEdgeWeightFunction();

  // Initialize the shortest distance computer
  m_DistanceMapper = new VTKMeshShortestDistance();

  // Initialize the Voronoi diagram
  m_VoronoiDiagram = new VTKMeshVoronoiDiagram();
  
  // Create the necessary filters
  m_DataReader = vtkPolyDataReader::New();
  // m_Stripper = vtkStripper::New();
  m_Triangulator = vtkTriangleFilter::New();
  m_NormalsFilter = vtkPolyDataNormals::New();
  m_CleanFilter = vtkCleanPolyData::New();
  
  // m_DisplayMesh = 
  m_Mesh = NULL;
  m_FocusPoint = NO_FOCUS;
  m_FocusCurve = NO_FOCUS;
}

TracerData
::~TracerData()
{
  delete m_EdgeWeightFunction;
  delete m_DistanceMapper;
  delete m_VoronoiDiagram;
  
  m_DataReader->Delete();
  // m_Stripper->Delete();
  m_Triangulator->Delete();
  m_NormalsFilter->Delete();
  m_CleanFilter->Delete();

  // Clean up marker data
  RemoveMarkerData();
}

void
TracerData
::LoadInputMesh(const char *file)
{
  // Read the data from disk
  m_DataReader->SetFileName(file);
  m_DataReader->Update();
  m_Mesh = m_DataReader->GetOutput();

  // Compute all normals 
  /*
  m_NormalsFilter->SetInput(m_Mesh);
  m_NormalsFilter->ConsistencyOn();
  m_NormalsFilter->AutoOrientNormalsOn();
  m_NormalsFilter->NonManifoldTraversalOn();
  m_NormalsFilter->Update();
  m_Mesh = m_NormalsFilter->GetOutput();
  */

  // Clean the input data
  //m_CleanFilter->SetInput(m_Mesh);
  //m_CleanFilter->SetTolerance(0);
  //m_CleanFilter->Update();
  //m_Mesh = m_CleanFilter->GetOutput();

  // Convert the input to triangles
  m_Triangulator->PassLinesOff();
  m_Triangulator->PassVertsOff();
  m_Triangulator->SetInput(m_Mesh);
  m_Triangulator->Update();
  m_Mesh = m_Triangulator->GetOutput();

  // Build the cells and links
  m_Mesh->BuildCells();
  m_Mesh->BuildLinks();

  // Create a marker array and assign null marker to all cells
  vtkIdTypeArray *xArray = vtkIdTypeArray::New();
  xArray->SetNumberOfValues(m_Mesh->GetNumberOfCells());
  for(vtkIdType iCell = 0; iCell < m_Mesh->GetNumberOfCells(); iCell++)
    xArray->SetValue(iCell, NO_MARKER);
  xArray->SetName("markers");
  m_Mesh->GetCellData()->AddArray(xArray);

  // Clean up the markers too
  RemoveMarkerData();

  // Create a clear label marker and update it's mesh
  AddMarker(NO_MARKER, "no label", 0.5, 0.5, 0.5);
  m_Markers[NO_MARKER]->m_NumberOfCells = m_Mesh->GetNumberOfCells();
  m_Markers[NO_MARKER]->UpdateMesh(m_Mesh);

  // Set up the stripper
  // m_Stripper->SetInput(m_Mesh);
  // m_Stripper->Update();
  // m_DisplayMesh = m_Stripper->GetOutput();

  // Send the data to the distance mapper
  m_DistanceMapper->SetInputMesh(m_Mesh);
  m_DistanceMapper->ComputeGraph();

  // Send the data to the Voronoi diagram
  m_VoronoiDiagram->SetInputMesh(m_Mesh);
  m_VoronoiDiagram->ComputeGraph();
  
  // Clean up the curves, because they are no longer valid
  m_Curves.RemoveAllData();

  // Unset current curve and current point
  SetFocusPoint(NO_FOCUS); SetFocusCurve(NO_FOCUS);

  // Notify listeners that the curves have changed
  TracerDataEvent evt(this);
  BroadcastOnCurveListChange(&evt);

  // Nofity of the change in the mesh
  BroadcastOnMeshChange(&evt);
}

void
TracerData
::RemoveMarkerData()
{
  list<vtkIdType> lMarkers;
  this->GetMarkerIds(lMarkers);
  list<vtkIdType>::iterator it = lMarkers.begin();
  while(it != lMarkers.end())
    DeleteMarker(*it++);
}

void
TracerData
::SaveCurves(const char *file)
{
  ofstream fout(file,ios_base::out);
  m_Curves.SaveAsText(fout);
  fout.close();
}

bool
TracerData
::LoadCurves(const char *file)
{
  ifstream fin(file,ios_base::in);
  bool rc = m_Curves.LoadAsText(fin);
  fin.close();

  if(rc)
    {
    // Reset the current point and curve
    SetFocusPoint(NO_FOCUS);
    SetFocusCurve(NO_FOCUS);

    // Notify listeners that the curves have changed
    TracerDataEvent evt(this);
    BroadcastOnCurveListChange(&evt);
    }

  return rc;
}

void
TracerData
::SetFocusCurve(int inFocusCurve)
{
  if(inFocusCurve != m_FocusCurve)
    {
    // Set the new focus curve
    m_FocusCurve = inFocusCurve;

    // Fire the update event
    TracerDataEvent evt(this);
    BroadcastOnFocusCurveChange(&evt);
    }
}

void 
TracerData
::ComputeDistances(int inFocusPoint)
{
  // Get the associated mesh vertex
  vtkIdType iVertex = m_Curves.GetControlPointVertex(inFocusPoint);

  // Time the computation
  double tStart = (double) clock();

  // Compute shortest distances to that point
  m_DistanceMapper->ComputeDistances(iVertex);

  // Get the elapsed time
  double tElapsed = (clock() - tStart) * 1000.0 / CLOCKS_PER_SEC;    
  cout << "Shortest paths to source " << iVertex 
    << " computed in " << tElapsed << " ms." << endl;
}

void
TracerData
::ComputeMarkerSegmentation()
{
  // Create a list of edges
  list<pair<vtkIdType, vtkIdType> > lEdges;

  // Pass the separating edges as boundaries for the segmentation
  TracerCurves::IdList lCurves;
  m_Curves.GetCurveIdList(lCurves);
  TracerCurves::IdList::iterator itCurve = lCurves.begin();
  while(itCurve != lCurves.end())
    {
    // Get the links associated with this curve
    const TracerCurves::MeshCurve &curve = m_Curves.GetCurveVertices(*itCurve);
    if(curve.size() > 1)
      {
      TracerCurves::MeshCurve::const_iterator it1 = curve.begin();
      TracerCurves::MeshCurve::const_iterator it2 = it1; ++it2;
      while(it2 != curve.end())
        {
        if(*it1 < *it2)
          lEdges.push_back(make_pair(*it1, *it2));
        else
          lEdges.push_back(make_pair(*it2, *it1));
        ++it1; ++it2;
        }
      }
    
    ++itCurve;
    }

  // Pass the curves to the partitioner
  m_VoronoiDiagram->SetBarrierEdges(lEdges);

  // Time the computation
  double tStart = (double) clock();

  // Compute shortest distances to that point
  list<vtkIdType> lMarkers;
  GetMarkerIds(lMarkers);
  m_VoronoiDiagram->ComputeVoronoiDiagram(lMarkers);

  // Get the elapsed time
  double tElapsed = (clock() - tStart) * 1000.0 / CLOCKS_PER_SEC;    
  cout << "Voronoi segmentation computed in " << tElapsed << " ms." << endl;

  // Clear the number-of-cells counters for each marker, including the NOMARKER one
  MarkerMap::iterator itMarker = m_Markers.begin();
  while(itMarker != m_Markers.end())
    itMarker++->second->m_NumberOfCells = 0;

  // Assign scalar values to the points in the mesh
  vtkIdTypeArray *xArray = (vtkIdTypeArray *)
    m_Mesh->GetCellData()->GetArray("markers");
  for(vtkIdType iCell=0; iCell<m_Mesh->GetNumberOfCells(); iCell++)
    {
    // Get the marker for the given cell
    vtkIdType iSource = m_VoronoiDiagram->GetVertexSource(iCell);
    xArray->SetValue(iCell, iSource);
    m_Markers[iSource]->m_NumberOfCells++;
    }

  // Update the individual meshes
  itMarker = m_Markers.begin();
  while(itMarker != m_Markers.end())
    {
    cout << "Marker " << itMarker->first << " has " << itMarker->second->m_NumberOfCells << " cells" << endl;
    itMarker++->second->UpdateMesh(m_Mesh);
    }

  // Fire event: segmentation changed
  TracerDataEvent evt(this);
  BroadcastOnSegmentationChange(&evt);
}

void 
TracerData
::SetFocusPoint(int inFocusPoint)
{
  if(inFocusPoint != m_FocusPoint)
    {
    // Compute distance to the new focus point
    if(inFocusPoint != NO_FOCUS)
      ComputeDistances(inFocusPoint);

    // Set the new focus point
    m_FocusPoint = inFocusPoint;

    // Fire the update event
    TracerDataEvent evt(this);
    BroadcastOnFocusPointChange(&evt);
    }
}

void 
TracerData
::UpdateEdgeWeightFunction(MeshEdgeWeightFunction *fnNew)
{
  // Pass on the new function
  m_DistanceMapper->SetEdgeWeightFunction(fnNew);

  // Recompute the graph
  m_DistanceMapper->ComputeGraph();
      
  // If there is a focus point, compute distances to it
  if(m_FocusPoint != NO_FOCUS)
    ComputeDistances(m_FocusPoint);

  // Delete the old edge weight function
  delete m_EdgeWeightFunction;

  // Assign the edge weight function 
  m_EdgeWeightFunction = fnNew;

  // Fire the appropriate event
  TracerDataEvent evt(this);
  BroadcastOnEdgeWeightsUpdate(&evt);
}
