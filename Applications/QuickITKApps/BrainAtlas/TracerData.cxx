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
  m_CleanFilter->SetInput(m_Mesh);
  m_CleanFilter->SetTolerance(0.0);
  m_CleanFilter->Update();
  m_Mesh = m_CleanFilter->GetOutput();

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
  xArray->SetName("markers");
  m_Mesh->GetCellData()->AddArray(xArray);

  // Clean up the markers too
  ResetMarkerState();

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
  // Get all the marker data in the mesh
  MarkerIterator it = m_Markers.begin();
  while(it != m_Markers.end())
    {
    if(it->first == NO_MARKER)
      {
      delete it->second;
      m_Markers.erase(it->first);
      }
    else
      DeleteMarker(it->first);
    ++it;
    }
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
    // For old-version meshes / deformed meshes - compute the coordinates 
    // of the control points
    TracerCurves::IdList lCurveIds;
    m_Curves.GetCurveIdList(lCurveIds);
    TracerCurves::IdList::iterator itCurve = lCurveIds.begin();
    while(itCurve != lCurveIds.end())
      {
      // Get all the control points in the curve
      const TracerCurves::IdList &lControls = m_Curves.GetCurveControls(*itCurve++);
      TracerCurves::IdList::const_iterator itControls = lControls.begin();
      while(itControls != lControls.end())
        {
        // Get the coordinates of the control point
        vtkIdType id = m_Curves.GetControlPointVertexId(*itControls);
        Vec x(m_Mesh->GetPoint(id));
        m_Curves.SetControlPointPosition(*itControls++, x);
        }
      }

    // Remove all marker data
    ResetMarkerState();

    // Register each of the markers that has been loaded
    TracerCurves::IdList lMarkers;
    m_Curves.GetMarkerIdList(lMarkers);
    TracerCurves::IdList::iterator itMarker = lMarkers.begin();
    while(itMarker != lMarkers.end())
      RegisterMarker(*itMarker++);

    // Reset the current point and curve
    SetFocusPoint(NO_FOCUS);
    SetFocusCurve(NO_FOCUS);
    SetCurrentMarker(NO_FOCUS);

    // Notify listeners that the curves have changed
    TracerDataEvent evt(this);
    BroadcastOnCurveListChange(&evt);
    BroadcastOnMarkerListChange(&evt);
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
  vtkIdType iVertex = m_Curves.GetControlPointVertexId(inFocusPoint);

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

  // Get the list of all markers
  TracerCurves::IdList lMarkerIds;
  m_Curves.GetMarkerIdList(lMarkerIds);

  // Create a list of marker cell faces
  list<vtkIdType> lMarkerCells;
  TracerCurves::IdList::iterator itMarkerId = lMarkerIds.begin();
  while(itMarkerId != lMarkerIds.end())
    lMarkerCells.push_back(m_Curves.GetMarkerFace(*itMarkerId++));

  // Compute shortest distances to that point
  m_VoronoiDiagram->ComputeVoronoiDiagram(lMarkerCells);

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
    TracerCurves::IdType iSourceMarker = m_CellToMarker[iSource];

    // Assign the marker to the cell as additional data
    xArray->SetValue(iCell, iSourceMarker);

    // Increase the number of cells assigned to the source marker
    m_Markers[iSourceMarker]->m_NumberOfCells++;
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

void
TracerData
::ResetMarkerState()
{
  // Remove all markers except marker -1
  RemoveMarkerData();

  // Reset the label associated with each cell to -1
  vtkIdTypeArray *xArray = 
    (vtkIdTypeArray *) m_Mesh->GetCellData()->GetArray("markers");
  for(vtkIdType iCell = 0; iCell < m_Mesh->GetNumberOfCells(); iCell++)
    xArray->SetValue(iCell, NO_MARKER);

  // Associate mesh data with the 'clear' label
  m_Markers[NO_MARKER] = new SegmentationMarkerData();
  m_Markers[NO_MARKER]->m_Id = NO_MARKER;
  m_CellToMarker[NO_MARKER] = NO_MARKER;
  m_Markers[NO_MARKER]->m_NumberOfCells = m_Mesh->GetNumberOfCells();
  m_Markers[NO_MARKER]->UpdateMesh(m_Mesh);
}

void
TracerData
::RegisterMarker(TracerCurves::IdType idMarker)
{
  // Disassociate information previously associated with this id
  DeleteMarker(idMarker);

  // Add the mesh information associated with this marker
  m_Markers[idMarker] = new SegmentationMarkerData();
  m_Markers[idMarker]->m_Id = idMarker;
  m_Markers[idMarker]->m_NumberOfCells = 0;

  // Update the mesh for the marker
  m_Markers[idMarker]->UpdateMesh(m_Mesh);

  // Add a mapping from the cell to the marker
  vtkIdType idCell = m_Curves.GetMarkerFace(idMarker);
  m_CellToMarker[idCell] = idMarker;
}

bool
TracerData
::ImportCurves(const char *file)
{
  // Curve import works like a load, except that none of the id's into the
  // mesh can be trusted. Therefore, they are reconstructed by finding closest
  // vertices in the actual mesh, and finding shortest paths between the 
  // consequently added points.
  
  // Load a new set of curves
  TracerCurves xImportCurves;
  ifstream fin(file,ios_base::in);
  bool rc = xImportCurves.LoadAsText(fin);
  fin.close();

  // Return false if nothing happened
  if(!rc) return false;

  // Clear the current set of curves
  m_Curves.RemoveAllData();

  // Import each curve
  TracerCurves::IdList lCurveIds;
  xImportCurves.GetCurveIdList(lCurveIds);
  TracerCurves::IdList::iterator itCurve = lCurveIds.begin();
  while(itCurve != lCurveIds.end())
    {
    // Give the user some messages, since this takes a while
    cout << "Importing curve " << *itCurve << ": \t"
      << xImportCurves.GetCurveName(*itCurve) << endl;
    
    // Add the curve by name to the data
    this->AddNewCurve(xImportCurves.GetCurveName(*itCurve));

    // Get all the control points in the curve
    const TracerCurves::IdList &lControls = xImportCurves.GetCurveControls(*itCurve);
    TracerCurves::IdList::const_iterator itControls = lControls.begin();
    while(itControls != lControls.end())
      {
      // Get the coordinates of the control point
      Vec xPoint = xImportCurves.GetControlPointPosition(*itControls);

      // Find the closest point to the coordinate
      vtkIdType iFound = m_DistanceMapper->FindClosestVertexInSpace(xPoint);
      Vec xFound(m_Mesh->GetPoint(iFound));

      // Report the mapping
      cout << "-- Vertex " << xPoint << " mapped to " << xFound << endl;

      // Add the control point to the curve using this class, so that the 
      // shortest distances will be computed
      this->AddNewPoint(iFound);

      ++itControls;
      }

    ++itCurve;
    }

  // Remove the markers currently there
  ResetMarkerState();

  // Import the markers
  TracerCurves::IdList lMarkers;
  xImportCurves.GetMarkerIdList(lMarkers);
  TracerCurves::IdList::iterator itMarker = lMarkers.begin();
  while(itMarker != lMarkers.end())
    {
    // Get the coordinate for the marker
    Vec xCenter = xImportCurves.GetMarkerPosition(*itMarker);

    // Get the color of the marker
    Vec xColor = xImportCurves.GetMarkerColor(*itMarker);

    // Find the cell closest to the marker
    vtkIdType iCell = m_DistanceMapper->FindClosestCellInSpace(xCenter);

    // Add the marker using the normal process
    this->AddMarker(
      iCell, xImportCurves.GetMarkerName(*itMarker), 
      xColor[0], xColor[1], xColor[2]);

    ++itMarker;
    }

  // Reset the current point and curve
  SetFocusPoint(NO_FOCUS);
  SetFocusCurve(NO_FOCUS);

  // Notify listeners that the curves have changed
  TracerDataEvent evt(this);
  BroadcastOnCurveListChange(&evt);
  BroadcastOnMarkerListChange(&evt);
  
  return true;
}
