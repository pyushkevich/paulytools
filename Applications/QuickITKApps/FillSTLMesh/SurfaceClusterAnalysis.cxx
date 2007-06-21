#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "ReadWriteVTK.h"
#include "vtkTriangleFilter.h"
#include "vtkBandedPolyDataContourFilter.h"
#include "vtkClipPolyData.h"
#include "vtkThresholdPoints.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCell.h"
#include "vtkTriangle.h"

#include "Registry.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int usage()
{
  cout << "This program performs cluster analysis on a VTK mesh" << endl;
  cout << "Usage: " << endl;
  cout << "    meshcluster configfile.txt DataArray Threshold Suffix" << endl;
  return -1;
}

struct Cluster
{
  // Surface area of the cluster
  double area;

  // Power of the thresholded data over the cluster
  double power;

  // Cluster size in nodes (?)
  size_t n;

  // p-values
  double pArea, pPower;

  // Dummy constructor
  Cluster() : n(0), area(0.0), power(0.0) {};
};

typedef std::vector<Cluster> ClusterArray;

/*
class ClusterGenerator
{
public:
  // Cluster definition
  typedef vector<Cluster> ClusterArray;

  // Constructor
  ClusterGenerator(vtkPolyData *mesh);

  // Compute clusters
  ClusterArray ComputeClusters(const char *data, const char *aelt, double thresh);

  // Get the output mesh
  vtkPolyData *GetOutputMesh()
    { return fConnect->GetOutput(); }

private:
  // Pipeline elements
  vtkClipPolyData *fContour;
  vtkThreshold *fThresh;
  vtkPolyDataConnectivityFilter *fConnect;

  // Meshes
  vtkPolyData *mesh;
};

ClusterGenerator::ClusterGenerator(vtkPolyData *inmesh)
{
  // Store the mesh
  mesh = inmesh;

  // Generate the pipeline
  fContour = vtkClipPolyData::New();
  fContour->SetInput(mesh);

  // fThresh = vtkThreshold::New();
  // fThresh->SetInput(fContour->GetOutput());
  // fThresh->ThresholdByLower(0);
  // fThresh->SetAttributeModeToUsePointData();

  static double srange[] = {-0.5, 0.5};
  
  fConnect = vtkPolyDataConnectivityFilter::New();
  fConnect->SetInput(fContour->GetOutput());
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  // fConnect->ScalarConnectivityOn();
  fConnect->SetScalarRange(srange);
}

ClusterGenerator::ClusterArray
ClusterGenerator::ComputeClusters(const char *data, const char *aelt, double thresh)
{
  // Update the mesh
  // mesh->GetPointData()->SetActiveScalars(data);
  // mesh->GetPointData()->CopyAllOn();

  // Compute the clusters
  fContour->SetValue(thresh);
  fConnect->Update();
  vtkPolyData *f = fContour->GetOutput();
  vtkPolyData *p = fConnect->GetOutput();

  vtkDataArray *daRegion = p->GetPointData()->GetScalars();
  vtkDataArray *daData = f->GetPointData()->GetArray(data);
  vtkDataArray *daArea = f->GetPointData()->GetArray(aelt);

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(size_t i = 0; i < p->GetNumberOfPoints(); i++)
    {
    size_t region = (size_t) (daRegion->GetTuple1(i));
    double x = daData->GetTuple1(i);
    double a = daArea->GetTuple1(i);
    ca[region].n++;
    ca[region].area += a;
    ca[region].power += a * x;
    }

  
  //for(size_t c = 0; c < ca.size(); c++)
  //  {
  //  printf("Cluster %d: n = %d, area = %f, power = %f, mean_t = %f\n",
  //    c, ca[c].n, ca[c].area, ca[c].power, ca[c].power / ca[c].area);
  //  }
    


  // Return the cluster array
  return ca;
}*/

ClusterArray ComputeClusters(
  vtkPolyData *mesh,
  const char *data, 
  double thresh, 
  vtkPolyData **mout = NULL)
{
  // Initialize mesh
  mesh->GetPointData()->SetActiveScalars(data);

  // Clip the data field at the threshold
  vtkClipPolyData *fContour = vtkClipPolyData::New();
  fContour->SetInput(mesh);
  fContour->SetValue(thresh);
  vtkPolyData *f = fContour->GetOutput();

  // Get the connected components
  vtkPolyDataConnectivityFilter * fConnect = vtkPolyDataConnectivityFilter::New();
  fConnect->SetInput(fContour->GetOutput());
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();
  vtkPolyData *p = fConnect->GetOutput();

  // Create output data arrays for computing area element
  vtkFloatArray *daArea = vtkFloatArray::New();
  daArea->SetName("area_element");
  daArea->SetNumberOfComponents(1);
  daArea->SetNumberOfTuples(p->GetNumberOfPoints());
  daArea->FillComponent(0, 0.0);

  // Compute the area of each triangle in the cluster set
  for(size_t k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
      throw("Wrong cell type");

    // Compute the area of the triangle
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    double p0[3], p1[2], p2[3];
    p->GetPoint(a0, p0);
    p->GetPoint(a1, p1);
    p->GetPoint(a2, p2);

    double area = vtkTriangle::TriangleArea(p0, p1, p2);

    // Split the area between neighbors
    daArea->SetTuple1(a0, area / 3.0 + daArea->GetTuple1(a0));
    daArea->SetTuple1(a1, area / 3.0 + daArea->GetTuple1(a1));
    daArea->SetTuple1(a2, area / 3.0 + daArea->GetTuple1(a2));
    }

  // The the important arrays in the resulting meshes
  vtkDataArray *daRegion = p->GetPointData()->GetScalars();
  vtkDataArray *daData = f->GetPointData()->GetArray(data);

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(size_t i = 0; i < p->GetNumberOfPoints(); i++)
    {
    size_t region = (size_t) (daRegion->GetTuple1(i));
    double x = daData->GetTuple1(i);
    double a = daArea->GetTuple1(i);
    ca[region].n++;
    ca[region].area += a;
    ca[region].power += a * x;
    }

  // Get the output if needed
  if(mout != NULL)
    {
    p->GetPointData()->AddArray(daArea);
    *mout = p;
    }
  else
    {
    // Delete the intermediates
    daArea->Delete();
    fConnect->Delete();
    fContour->Delete();
    }

  // Return the cluster array
  return ca;
}


void ComputeTTest(
  vtkPolyData *pd, 
  const char *var, 
  const vector<int> &labels,
  int l1, int l2)
{
  // Get the pointdata
  vtkDataArray *data = pd->GetPointData()->GetArray(var);
  vtkDataArray *ttest = pd->GetPointData()->GetArray("ttest");

  // Get the sizes of the cohorts
  int n1 = 0, n2 = 0;
  for(size_t j = 0; j < labels.size(); j++)
    {
    if(labels[j] == l1) 
      n1++;
    else    
      if(labels[j] == l2) 
        n2++;
    }

  // Loop over all the points
  for(size_t i = 0; i < pd->GetNumberOfPoints(); i++)
    {
    // Compute class-wise sums and sums of squares
    double s1 = 0, s2 = 0, ss1 = 0, ss2 = 0;
    for(size_t j = 0; j < labels.size(); j++)
      {
      double x = data->GetComponent(i, j);
      if(labels[j] == l1)
        { s1 += x; ss1 += x * x; }
      else if(labels[j] == l2)
        { s2 += x; ss2 += x * x; } 
      }

    // Compute t2 directly from sums and sums of squares
    double m1 = s1 / n1, m2 = s2 / n2;
    double v1 = ss1 - s1 * m1, v2 = ss2 - s2 * m2;
    double den = sqrt((v1+v2) * (1./n1 + 1./n2) / (n1 + n2 - 2.));
    double t = den > 0 ? (m1 - m2) / den : 0;

    // Add the t2 value to the array
    ttest->SetTuple1(i, t);
    }
}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  // Read the registry file
  Registry registry(argv[1]);

  // Get the number of permutations
  size_t np = registry["Analysis.NumberOfPermutations"][1000];

  // Get the cohort labels
  vector<int> labels = registry.Folder("Cohort").GetArray(-1);
  vector<int> true_labels = labels;

  // Get the variables of interest
  string sVOI = registry["Analysis.TestVariable"][""];
  if(argc > 2) sVOI = argv[2];
  cout << "Variable of Interest: " << sVOI << endl;

  // Get the clustering threshold
  double thresh = registry["Analysis.ClusteringThreshold"][0.0];
  if(argc > 3) thresh = atof(argv[3]);
  cout << "Threshold: " << thresh << endl;

  // If the threshold is negative, simply change the direction of the test
  int l1 = 0, l2 = 1;
  if(thresh < 0)
    { l1 = 1; l2 = 0; thresh = -thresh; }

  // Load each of the meshes and create a cluster analyzer
  vector<string> fnMeshes, fnOutMeshes;
  vector<vtkPolyData *> mesh;
  vector<ClusterArray> claTrue;
  fnMeshes = registry.Folder("Mesh").GetArray(string(""));
  fnOutMeshes = registry.Folder("OutputMesh").GetArray(string(""));
  if(fnMeshes.size() == 0)
    { cerr << "Missing mesh specification" << endl; return -1; }


  // If there is a suffix, add it to the output meshes
  if(argc > 4) 
    {
    string suffix = argv[4];
    for(size_t i = 0; i < fnOutMeshes.size(); i++)
      fnOutMeshes[i] = fnOutMeshes[i] + suffix;
    }
  
  // Read the meshes
  for(size_t i = 0; i < fnMeshes.size(); i++)
    {
    // Read mesh
    cout << "Reading mesh " << fnMeshes[i] << endl;
    mesh.push_back(ReadVTKData(fnMeshes[i].c_str()));

    // Create a t-test array
    vtkFloatArray *array = vtkFloatArray::New();
    array->SetName("ttest");
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(mesh[i]->GetNumberOfPoints());
    mesh[i]->GetPointData()->AddArray(array);
    mesh[i]->GetPointData()->SetActiveScalars("ttest");
    }

  // Run permutation analysis
  vector<double> hArea(np), hPower(np);
  for(size_t ip = 0; ip < np; ip++)
    {
    // Apply a random permutation to the labels array
    random_shuffle(labels.begin(), labels.end());
    hArea[ip] = 0; hPower[ip] = 0;

    // Build up the histogram of cluster areas (and powers)
    for(size_t i = 0; i < fnMeshes.size(); i++)
      {
      // For each mesh, compute the t-test and the clusters
      ComputeTTest(mesh[i], sVOI.c_str(), labels, l1, l2);
      ClusterArray ca = ComputeClusters(mesh[i], "ttest", thresh);

      // Now find the largest cluster
      for(size_t c = 0; c < ca.size(); c++)
        {
        if(ca[c].area > hArea[ip]) hArea[ip] = ca[c].area;
        if(ca[c].power > hPower[ip]) hPower[ip] = ca[c].power;
        }
      }
    cout << "." << flush;
    if((ip+1) % 80 == 0 || (ip+1) == np) 
      cout << " " << ip+1 << endl;
    }

  // Sort the histograms
  sort(hArea.begin(), hArea.end());
  sort(hPower.begin(), hPower.end());

  // Going back to the original meshes, assign a cluster p-value to each mesh
  for(size_t i = 0; i < fnMeshes.size(); i++)
    {
    vtkPolyData *mout;

    // Compute the t-test for the mesh with the correct labels
    ComputeTTest(mesh[i], sVOI.c_str(), true_labels, l1, l2);

    // Compute the clusters for the mesh with the correct labels
    ClusterArray ca = ComputeClusters(mesh[i], "ttest", thresh, &mout);
    printf("MESH %s HAS %d CLUSTERS: \n", fnMeshes[i].c_str(), ca.size());

    // Assign a p-value to each cluster
    for(size_t c = 0; c < ca.size(); c++)
      {
      // Brute force search in the histogram :(
      size_t zArea = 0, zPower = 0;
      while(zArea < np && hArea[zArea] < ca[c].area) zArea++;
      while(zPower < np && hPower[zPower] < ca[c].power) zPower++;
      ca[c].pArea = 1.0 - zArea * 1.0 / np;
      ca[c].pPower = 1.0 - zPower * 1.0 / np;
      bool sig = (ca[c].pArea <= 0.05 || ca[c].pPower <= 0.05);
      printf("Cluster %03d:  Area = %6f (p = %6f);  Power = %6f (p = %6f); %s\n",
        c, ca[c].area, ca[c].pArea, ca[c].power, ca[c].pPower, 
        sig ? "***" : "");
      }

    // Create output mesh arrays for p-values
    string snArea = sVOI + string(" p-cluster-area");
    vtkFloatArray *aArea = vtkFloatArray::New();
    aArea->SetName(snArea.c_str());
    aArea->SetNumberOfComponents(1);
    aArea->SetNumberOfTuples(mout->GetNumberOfPoints());
    mout->GetPointData()->AddArray(aArea);

    string snPower = sVOI + string(" p-cluster-power");
    vtkFloatArray *aPower = vtkFloatArray::New();
    aPower->SetName(snPower.c_str());
    aPower->SetNumberOfComponents(1);
    aPower->SetNumberOfTuples(mout->GetNumberOfPoints());
    mout->GetPointData()->AddArray(aPower);

    // Set the mesh arrays' p-values
    for(size_t p = 0; p < mout->GetNumberOfPoints(); p++)
      {
      size_t r = (size_t) mout->GetPointData()->GetScalars()->GetTuple1(p);
      aArea->SetTuple1(p, ca[r].pArea);
      aPower->SetTuple1(p, ca[r].pPower);
      }

    // Threshold out the non-significant clusters
    /*
    mout->GetPointData()->SetActiveScalars(snPower.c_str());
    vtkClipPolyData *fClip = vtkClipPolyData::New();
    fClip->SetInput(mout);
    fClip->SetValue(0.05);
    fClip->InsideOutOn();
    fClip->Update();
    */

    // Save the output mesh
    WriteVTKData(mout, fnOutMeshes[i].c_str());
    }
}
