#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkClipDataSet.h>
#include <vtkConnectivityFilter.h>
#include "vtkPointData.h"
#include "ReadWriteVTK.h"
#include "vtkTriangleFilter.h"
#include "vtkBandedPolyDataContourFilter.h"
#include "vtkClipPolyData.h"
#include "vtkThresholdPoints.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkTriangle.h"

#include "Registry.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int usage()
{
  cout << "This program performs cluster analysis on a VTK mesh (both surface and volume meshes)" << endl;
  cout << "Usage: " << endl;
  cout << "    meshcluster configfile.txt DataArray Threshold Suffix TypeOfTest" << endl;
  cout << "    TypeOfTest = 0 or none means two sample t-test" << endl;
  cout << "    TypeOfTest = 1 means correlation t-test" << endl;
  cout << "    TypeOfTest = 2 means paired t-test" << endl;
  cout << "    TypeOfTest = 3 means correlation with a summary measure" << endl;
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

  // t-value
  double tvalue;
 
  // Dummy constructor
  Cluster() : n(0), area(0.0), power(0.0), tvalue(0.0) {};
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

template <class TMeshType>
TMeshType * ReadMesh(const char *fname)
{ return NULL; }

template <>
vtkUnstructuredGrid *ReadMesh<>(const char *fname)
{
  vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <>
vtkPolyData *ReadMesh<>(const char *fname)
{
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}


template <class TMeshType>
void WriteMesh(TMeshType *mesh, const char *fname)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}


template <class TMeshType>
ClusterArray ComputeClusters(
  TMeshType *mesh,
  const char *data, 
  double thresh, 
  TMeshType **mout = NULL)
{ ClusterArray ca(1) ; return ca; }


template <>
ClusterArray ComputeClusters(
  vtkPolyData *mesh,
  const char *data, 
  double thresh, 
  vtkPolyData **mout )
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
    ca[region].tvalue += x;
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

template <>
ClusterArray ComputeClusters(
  vtkUnstructuredGrid *mesh,
  const char *data, 
  double thresh, 
  vtkUnstructuredGrid **mout )
{
  // Initialize mesh
  mesh->GetPointData()->SetActiveScalars(data);

  // Clip the data field at the threshold
  vtkClipDataSet *fContour = vtkClipDataSet::New();
  fContour->SetInput(mesh);
  fContour->SetValue(thresh);
  vtkUnstructuredGrid *f = fContour->GetOutput();

  // Get the connected components
  vtkConnectivityFilter * fConnect = vtkConnectivityFilter::New();
  fConnect->SetInput(fContour->GetOutput());
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();
  vtkUnstructuredGrid *p = fConnect->GetOutput();

  // Create output data arrays for computing volume element
  vtkFloatArray *daArea = vtkFloatArray::New();
  daArea->SetName("volume_element");
  daArea->SetNumberOfComponents(1);
  daArea->SetNumberOfTuples(p->GetNumberOfPoints());
  daArea->FillComponent(0, 0.0);

  // Compute the area of each tetra in the cluster set
  for(size_t k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
    if(cell->GetCellType() != VTK_TETRA)
      throw("Wrong cell type");

    // Compute the area of the tetra
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    vtkIdType a3 = cell->GetPointId(3);
    double p0[3], p1[2], p2[3], p3[3];
    p->GetPoint(a0, p0);
    p->GetPoint(a1, p1);
    p->GetPoint(a2, p2);
    p->GetPoint(a3, p3);

    double area = vtkTetra::ComputeVolume(p0, p1, p2, p3);

    // Split the area between neighbors
    daArea->SetTuple1(a0, area / 4.0 + daArea->GetTuple1(a0));
    daArea->SetTuple1(a1, area / 4.0 + daArea->GetTuple1(a1));
    daArea->SetTuple1(a2, area / 4.0 + daArea->GetTuple1(a2));
    daArea->SetTuple1(a3, area / 4.0 + daArea->GetTuple1(a3));
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
    ca[region].tvalue += x;
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


template <class TMeshType>
void ComputeTTest(
  TMeshType *pd, 
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

template <class TMeshType>
void ComputeCorrTest(
  TMeshType *pd, 
  const char *var, 
  const vector<int> &labels,
  const vector<int> &indiv_labels,
  int l1, int l2, int ispaired, const vector<float> &corrVar)
{
  // Get the pointdata
  vtkDataArray *data = pd->GetPointData()->GetArray(var);
  vtkDataArray *ttest = pd->GetPointData()->GetArray("ttest");

  // Create an r (correlation coefficient) array
  vtkFloatArray *array = vtkFloatArray::New();
  vtkDataArray *corrcoef ;
  if (ispaired == 1 || ispaired == 3)
  {
  	array->SetName("R");
  	array->SetNumberOfComponents(1);
  	array->SetNumberOfTuples(pd->GetNumberOfPoints());
  	pd->GetPointData()->AddArray(array);
  	corrcoef = pd->GetPointData()->GetArray("R");
  }

  int n;
  if (ispaired == 1 || ispaired == 2)
     n = indiv_labels.size()/2;
  else if (ispaired == 3)
     n = indiv_labels.size();
  // Loop over all the points
  for(size_t i = 0; i < pd->GetNumberOfPoints(); i++)
    {
    // Compute class-wise sums and sums of squares
    double s1 = 0, s2 = 0, ss1 = 0, ss2 = 0, s1s2 = 0, s12 = 0, ss12 = 0;
    for(size_t j = 0; j < n; j++)
      {
      double x = data->GetComponent(i, j);
      double y;
      if (ispaired == 1 || ispaired == 2)
         y = data->GetComponent(i, indiv_labels[j] + n);
      else if (ispaired == 3)
         y = corrVar[ indiv_labels[j] ];
        { s1 += x; ss1 += x * x; }
        { s2 += y; ss2 += y * y; } 
        s1s2 += x * y;
        if(labels[j] == l1)
	  { s12 += x - y; ss12 += (x - y)*(x - y); }
        else if(labels[j] == l2)
	  { s12 += y - x; ss12 += (x - y)*(x - y); }
	
      }

    double r = 0,t = 0;
    // Compute t2 directly from sums and sums of squares
    if (ispaired == 2) // paired t-test
    {
	double numerator = sqrt(n )*s12;
	double denominator = sqrt(n*ss12 - s12*s12);
	t = numerator/denominator;
    }
    else // paired correlation
    {
    	double numerator = n*s1s2 - s1 * s2;
    	double denominator = sqrt(n*ss1 - s1 * s1) * sqrt(n*ss2 - s2 * s2);
    	r = numerator/denominator;
    	t = r * sqrt((n-2) / (1 - r*r)); 
    	// If testing for negative correlation flip sign of t
    	if (l1 == 1) t = -t;
    }

    // Add the t2 value to the array
    ttest->SetTuple1(i, t);
    if (ispaired == 1 || ispaired == 3)
    	corrcoef->SetTuple1(i, r);
    }
}


template <class TMeshType>
int meshcluster(int argc, char *argv[], Registry registry, bool isPolyData)
{

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

  // Paired test desired. TODO registry ?
  int  ispaired = 0;
  if (argc > 5) ispaired = atoi(argv[5]);


  // SD If a correlation is desired, assume that the true labels are in the form of 000.. 111..
  // and that in each group the corresponding positions are paired, as in 123.. 123..
  // In this case, the permutation needs to be done only on one of the group's individual labels 123.. 
  // first create the individual labels with the above assumption
  
  
  int groupSize, Nlabels;
  vector<int> indiv_labels, true_indiv_labels ;
  vector<float> corrVar;
  if (ispaired == 1 ) // correlation between two variables 
  {
     Nlabels = (int)labels.size();
     if (Nlabels % 2 !=0)
        { cerr << "Total array size is odd, must be even for paired tests" << endl; return -1; }
     else
        groupSize = Nlabels/2;
     cout << "Generating individual labels for correlation, group size is " << groupSize << endl;
     for (int cohort = 0; cohort < 2; cohort++)
         for (int i=0; i< groupSize; i++) 
             indiv_labels.push_back( i ); 
     true_indiv_labels = indiv_labels;
  }
  else if (ispaired == 2) // paired t-test
  {
     Nlabels = (int)labels.size();
     if (Nlabels % 2 !=0)
        { cerr << "Total array size is odd, must be even for paired tests" << endl; return -1; }
     else
        groupSize = Nlabels/2;
     cout << "Generating individual labels for paired test, group size is " << groupSize << endl;
     for (int cohort = 0; cohort < 2; cohort++)
         for (int i=0; i< groupSize; i++)
             indiv_labels.push_back( i );
     true_indiv_labels = indiv_labels;
  }
  else if (ispaired == 3) // correlation with summary measurement
  {
     FILE *fd;
     fd = fopen("variable.txt","r");
     Nlabels = (int)labels.size();
     groupSize = Nlabels;
     float val;
     cout << "Variable values: " ;
     for (int i=0; i< groupSize; i++) 
     {
         indiv_labels.push_back( i );
         fscanf(fd, "%f\n", &val);
         corrVar.push_back( val ); 
         cout << val << " " ; 
     }
     true_indiv_labels = indiv_labels;
     cout << endl;
  }

  // If the threshold is negative, simply change the direction of the test
  string posneg("_pos");
  int l1 = 0, l2 = 1;
  if(thresh < 0)
    { l1 = 1; l2 = 0; thresh = -thresh; posneg = "_neg";}

  // Load each of the meshes and create a cluster analyzer
  vector<string> fnMeshes, fnOutMeshes;
  vector<ClusterArray> claTrue;
  fnMeshes = registry.Folder("Mesh").GetArray(string(""));
  fnOutMeshes = registry.Folder("OutputMesh").GetArray(string(""));
  if(fnMeshes.size() == 0)
    { cerr << "Missing mesh specification" << endl; return -1; }


  // If there is a suffix, add it to the output meshes
  if(argc > 4) 
    {
    string suffix = argv[4];
    string teststring;
    std::stringstream ss;
    ss << thresh;
    if (ispaired ==0)
    	teststring = "ttst";
    else if (ispaired == 1)
        teststring = "corr";
    else if (ispaired == 2)
    	teststring = "ptst";
    else if (ispaired == 3)
    	teststring = "vcorr";
    else
	teststring = "unkn"; 
    for(size_t i = 0; i < fnOutMeshes.size(); i++)
      fnOutMeshes[i] = fnOutMeshes[i] + teststring + ss.str() + posneg + suffix;
    }
  
  // Read the meshes
  vector<TMeshType *> mesh;
  for(size_t i = 0; i < fnMeshes.size(); i++)
    {
    // Read mesh
    cout << "Reading mesh " << fnMeshes[i] << endl;
    mesh.push_back(ReadMesh<TMeshType>(fnMeshes[i].c_str()));

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

    // SD shuffle individual member labels within group for correlations
    if (ispaired == 1 || ispaired == 3)
    {
       vector<int>::iterator it;
       random_shuffle(indiv_labels.begin(), indiv_labels.begin()+groupSize );
       /*
       for (it=indiv_labels.begin(); it!=indiv_labels.end(); ++it)
           cout << *it << " " << flush ;
       cout << endl;
       return 0;
       */
    }    

    // Build up the histogram of cluster areas (and powers)
    for(size_t i = 0; i < fnMeshes.size(); i++)
      {
      // SD paired tests
      if (ispaired > 0)
         ComputeCorrTest<TMeshType>(mesh[i], sVOI.c_str(), labels, indiv_labels, l1, l2, ispaired, corrVar);
      // For each mesh, compute the t-test and the clusters
      else
         ComputeTTest<TMeshType>(mesh[i], sVOI.c_str(), labels, l1, l2);
      ClusterArray ca = ComputeClusters<TMeshType>(mesh[i], "ttest", thresh);

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
    TMeshType *mout;

    // SD paired tests
    if (ispaired > 0)
       ComputeCorrTest<TMeshType>(mesh[i], sVOI.c_str(), true_labels, true_indiv_labels, l1, l2, ispaired, corrVar);
    else
    // Compute the t-test for the mesh with the correct labels
       ComputeTTest<TMeshType>(mesh[i], sVOI.c_str(), true_labels, l1, l2);

    // Compute the clusters for the mesh with the correct labels
    ClusterArray ca = ComputeClusters<TMeshType>(mesh[i], "ttest", thresh, &mout);
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
      printf("Cluster %03d:  AvgT = %6f; Area = %6f (p = %6f);  Power = %6f (p = %6f); %s\n",
        c, ca[c].tvalue/ca[c].n, ca[c].area, ca[c].pArea, ca[c].power, ca[c].pPower, 
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
    cout << "save the output ************ TODO" << endl;
    WriteMesh<TMeshType>(mout, fnOutMeshes[i].c_str());
    WriteMesh<TMeshType>(mesh[i], fnMeshes[i].c_str());
    }


}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();
  // Read the registry file
  Registry registry(argv[1]);

  vector<string> fnMeshes;
  fnMeshes = registry.Folder("Mesh").GetArray(string(""));
  if(fnMeshes.size() == 0)
    { cerr << "Missing mesh specification" << endl; return -1; }


  
  // Read the meshes
  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(fnMeshes[0].c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();

  bool isPolyData = true;
  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    isPolyData = false;
    return meshcluster<vtkUnstructuredGrid>( argc, argv, registry, isPolyData);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return meshcluster<vtkPolyData>( argc, argv, registry, isPolyData);

    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
