#include "TracerData.h"
#include <fstream>

TracerData
::TracerData()
{
  // Initialize the mesh edge weight function
  m_EdgeWeightFunction = new EuclideanDistanceMeshEdgeWeightFunction();

  // Initialize the shortest distance computer
  m_DistanceMapper = new VTKMeshShortestDistance();
  
  // Create the necessary filters
  m_DataReader = vtkPolyDataReader::New();
  m_Stripper = vtkStripper::New();
  m_Triangulator = vtkTriangleFilter::New();
  m_NormalsFilter = vtkPolyDataNormals::New();
  m_CleanFilter = vtkCleanPolyData::New();
  
  m_DisplayMesh = m_Mesh = NULL;
  m_FocusPoint = -1;
  m_FocusCurve = -1;
}

TracerData
::~TracerData()
{
  delete m_EdgeWeightFunction;
  delete m_DistanceMapper;
  
  m_DataReader->Delete();
  m_Stripper->Delete();
  m_Triangulator->Delete();
  m_NormalsFilter->Delete();
  m_CleanFilter->Delete();
  
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

  // Set up the stripper
  m_Stripper->SetInput(m_Mesh);
  m_Stripper->Update();
  m_DisplayMesh = m_Stripper->GetOutput();

  // Send the data to the distance mapper
  m_DistanceMapper->SetInputMesh(m_Mesh);
  m_DistanceMapper->ComputeGraph();
  
  // Clean up the curves, because they are no longer valid
  m_Curves.RemoveAllData();

  // Notify listeners that the curves have changed
  TracerDataEvent evt(this);
  BroadcastOnCurveListChange(&evt);

  // Nofity of the change in the mesh
  BroadcastOnMeshChange(&evt);
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
    SetFocusPoint(-1);
    SetFocusCurve(-1);

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
::SetFocusPoint(int inFocusPoint)
{
  if(inFocusPoint != m_FocusPoint)
    {
    // Compute distance to the new focus point
    if(inFocusPoint != -1)
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
  if(m_FocusPoint != -1)
    m_DistanceMapper->ComputeDistances(m_FocusPoint);

  // Delete the old edge weight function
  delete m_EdgeWeightFunction;

  // Assign the edge weight function 
  m_EdgeWeightFunction = fnNew;

  // Fire the appropriate event
  TracerDataEvent evt(this);
  BroadcastOnEdgeWeightsUpdate(&evt);
}
