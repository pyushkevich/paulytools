#include <iostream>
#include "METISTools.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

using namespace std;
using namespace itk;

int usage()
{
  const char *usage = 
  "usage: metisseg [options] input.img output.img num_part" 
  "\n   uses METIS to segment a binary image into num_part partitions" 
  "\noptions: " 
  "\n   -w X.X X.X          Specify relative weights of the partitions" 
  "\n   -t file.txt         Specify the intensity->graph hint file (see below) " 
  "\n   -p N1 N2 N3         define cut plane at dimension N1, slice N2" 
  "\n                       with relative edge strength N3" 
  "\n   -o                  use optimization to refine partition weights" 
  "\nhint files: " 
  "\n   The hint file is used to convert an image into a graph. It specifies "
  "\n   the weights assigned to the vertices and edges in the graph based on"
  "\n   the intensities of the pixels that connect the edges. The vertex "
  "\n   weight of zero means that the pixels with a given intensity are not"
  "\n   included in the graph. The hint file contains the following types of"
  "\n   entries."
  "\n   "
  "\n     V int wgt"
  "\n     E int1 int2 wgt"
  "\n   "
  "\n   The first entry specifies vertex weights. The second specifies edge "
  "\n   weights. You can use asterisk (*) as a wildcard for any intensity value."
  "\n   Weights can be specified as an integer value, or as a randomly generated"
  "\n   value. For the latter, use the notation 'U n1 n2' for the uniform dist."
  "\n   and 'N n1 n2' for the normal dist. with mean n1 and s.d. n2. The order in"
  "\n   which the rules are specified matters, as the later rules replace the"
  "\n   earlier ones.";
  
  cout << usage << endl;
  return -1;
}

/* ***************************************************************************
 * GLOBAL TYPE DEFINITIONS
 * *************************************************************************** */
typedef itk::Image< short, 3 > ImageType;
typedef ImageToGraphFilter< ImageType > GraphFilter;
typedef GraphFilter::WeightFunctorType BaseWeightFunctor;
typedef vnl_vector<float> Vec;

/* ***************************************************************************
 * WEIGHT TABLE CODE
 * *************************************************************************** */
class MyWeightFunctor : public BaseWeightFunctor
{
public:

  /** Read table from file */
  bool ReadTable(const char *file) {};

  /** Check inclusion */
  virtual bool IsPixelAVertex(short i1)
    {
    return i1 != 0;
    }
  
  /** Compute edge weight */
  virtual int GetEdgeWeight(short i1, short i2)
    {
    return 1;
    }
      
  /** Compute vertex weight */
  virtual int GetVertexWeight(short i1)
    {
    return 1;
    }
};

/* ***************************************************************************
 * GRAPH VERIFICATION
 * *************************************************************************** */
template<class T, class S>
void VerifyGraph(int n, T *ai, T *a, S *wv, S *we)
{
  for(T i=0; i < n; i++)
    {
    T k = ai[i+1] - ai[i];
    for(T j=0;j < k;j++)
      {
      // Neighbor of i
      T m = a[ai[i] + j];
      S w = we[ai[i] + j];

      // Check the match
      T l = ai[m+1] - ai[m];
      bool match = false;
      for(T p=0;p < l;p++)
        if(a[ai[m]+p] == i && we[ai[m]+p] == w)
          { match = true; break; }

      if(!match)
        {
        cout << "Mismatch at node " << i << " edge to " << m << " weight " << w << endl;
        }
      if(w <= 0 || w > 1000 || wv[i] <= 0 || wv[i] > 1000)
        {
        cout << " bad weight" << i << endl;
        }
      }
    }
}

Vec OptimizeMETISPartition(GraphFilter *fltGraph, const Vec &xWeights)
{
  // Create a METIS problem based on the graph and weights
  MetisPartitionProblem::Pointer mp = MetisPartitionProblem::New();
  mp->SetProblem(fltGraph, xWeights.size() - 1);

  // Get the starting solution
  MetisPartitionProblem::ParametersType x( xWeights.size() - 1 );
  for(unsigned int j = 0; j < x.size(); j++)
    x[j] = xWeights[j];
  
  typedef OnePlusOneEvolutionaryOptimizer Optimizer;
  Optimizer::Pointer opt = Optimizer::New();

  typedef itk::Statistics::NormalVariateGenerator Generator;
  Generator::Pointer generator = Generator::New();

  opt->SetCostFunction(mp);
  opt->SetInitialPosition(x);
  opt->SetInitialRadius(0.005);
  opt->SetMaximumIteration(100);
  opt->SetNormalVariateGenerator(generator);
  opt->StartOptimization();

  Vec xResult(xWeights.size(), 0.0);
  xResult[xWeights.size() - 1] = 1.0;
  for(unsigned int i = 0; i < xWeights.size() - 1; i++)
    {
    xResult[i] = x[i];
    xResult[xWeights.size() - 1] -= x[i];
    }
  
  return xResult;
}


int main(int argc, char *argv[])
{
  // Check arguments
  if(argc < 4) return usage();

  // Variables to hold command line arguments
  string fnInput(argv[argc-3]), fnOutput(argv[argc-2]);
  int nParts = atoi(argv[argc-1]);
  vnl_vector<float> xWeights(nParts,1.0 / nParts);
  int iPlaneDim = -1, iPlaneSlice = -1, iPlaneStrength = 10;
  bool flagOptimize = false;
  
  // Read the options
  for(unsigned int iArg=1;iArg<argc-3;iArg++)
    {
    if(!strcmp(argv[iArg],"-w"))
      {
      unsigned int iPart = atoi(argv[++iArg]);
      double xWeightPart = atof(argv[++iArg]);
      if(iPart < nParts)
        {
        xWeights[iPart] = xWeightPart;
        }
      else
        {
        cerr << "Incorrent index in -w parameter" << endl;
        return usage();
        }
      }
    else if(!strcmp(argv[iArg],"-p"))
      {
      iPlaneDim = atoi(argv[++iArg]);
      iPlaneSlice = atoi(argv[++iArg]);
      iPlaneStrength = atoi(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg],"-o"))
      {
      flagOptimize = true;
      }
    else if(!strcmp(argv[iArg],"-seed"))
      {
      srand(atoi(argv[++iArg]));
      }
    else
      {
      cerr << "unknown option!" << endl;
      return usage();
      }
    }

  // Write partition information
  cout << "will generate " << nParts << " partitions" << endl;
  float xWeightSum = 0.0f;
  for(unsigned int iPart = 0;iPart < nParts;iPart++)
    {
    cout << "   part " << iPart << "\t weight " << xWeights[iPart] << endl;
    xWeightSum += xWeights[iPart];
    }
  cout << "   total of weights : " << xWeightSum << endl;
  if(iPlaneDim >= 0)
    {
    cout << "will insert cut plane at slice " << iPlaneSlice 
      << " in dimension " << iPlaneDim 
      << " with strength " << iPlaneStrength << endl;
    }
  cout << endl;

  // Read the input image image
  cout << "reading input image" << endl;

  typedef ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fnInput.c_str());
  fltReader->Update();
  ImageType::Pointer img = fltReader->GetOutput();

  cout << "   image has dimensions " << img->GetBufferedRegion().GetSize() 
    << ", nPixels = " << img->GetBufferedRegion().GetNumberOfPixels() << endl;

  // Begin converting the image to a graph
  cout << "converting image to graph structure" << endl;

  // Create the weight functor
  MyWeightFunctor fnWeight;
  
  // Create the graph filter
  typedef ImageToGraphFilter<ImageType,idxtype> GraphFilter;
  GraphFilter::Pointer fltGraph = GraphFilter::New();
  fltGraph->SetInput(img);
  fltGraph->SetWeightFunctor(&fnWeight);
  fltGraph->Update();

  // If asked to optimize, compute the best set of weights
  if(flagOptimize)
    {
    // Run the experimental optimization
    xWeights = OptimizeMETISPartition(fltGraph, xWeights);
    }

  // Run METIS once, using the specified weights
  cout << "Running METIS ... " << flush; 
  int *iPartition = new int[fltGraph->GetNumberOfVertices()];
  int xCut = RunMETISPartition<ImageType>(
    fltGraph, xWeights.size(), xWeights.data_block(), iPartition );
  cout << "done. Cut value = " << xCut << endl;

  // Create output image 
  ImageType::Pointer imgOut = ImageType::New();
  imgOut->SetRegions(img->GetBufferedRegion());
  imgOut->Allocate();
  imgOut->FillBuffer(0);
    
  // Apply partition to output image
  for(unsigned int iVertex = 0; iVertex < fltGraph->GetNumberOfVertices(); iVertex++)
    {
    GraphFilter::IndexType idx = fltGraph->GetVertexImageIndex(iVertex);
    imgOut->SetPixel(idx, iPartition[iVertex] + 1);
    }

  // Delete the partition
  delete iPartition;
    
  // Write the image
  cout << "writing output image" << endl;
  
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetInput(imgOut);
  fltWriter->SetFileName(fnOutput.c_str());
  fltWriter->Update();

  // Done!
  return 0;
}
