/**
 * Read metis output, create an image
 */
#include "ReadWriteImage.h"
#include <iostream>

// #include "metis.h"
typedef int idxtype;
extern "C" {
  void METIS_WPartGraphKway(
  	int *, idxtype *, idxtype *, idxtype *, idxtype *, 
  	int *, int *, int *, float *, int *, int *, idxtype *);

  void METIS_WPartGraphRecursive(
  	int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, 
    int *, int *, float *, int *, int *, idxtype *);
}

using namespace std;
using namespace itk;

int usage()
{
  cout << "usage: metisseg [options] input.img output.img num_part" << endl;
  cout << "   uses METIS to segment a binary image into num_part partitions" << endl;
  cout << "options: " << endl;
  cout << "   -w N XX.X           weight of partition N is XX.X" << endl;
  cout << "   -p N1 N2 N3         define cut plane at dimension N1, slice N2" << endl;
  cout << "                       with relative edge strength N3" << endl;
  cout << "   -o                  use optimization to refine partition weights" << endl;
  return -1;
}

#include <itkConstantBoundaryCondition.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkSingleValuedCostFunction.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkNormalVariateGenerator.h>

template<class TImage, class TVertex = int> 
class BinaryImageToGraphFilter : public itk::ProcessObject
{
public:
  typedef BinaryImageToGraphFilter Self;
  typedef itk::ProcessObject Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef TImage ImageType;
  typedef typename TImage::IndexType IndexType;
  typedef TVertex VertexType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);
  
  /** Constructor */
  BinaryImageToGraphFilter()
    {
    m_Adjacency = NULL;
    m_AdjacencyIndex = NULL;
    m_ImageIndex = NULL;
    m_NumberOfVertices = 0;
    m_NumberOfEdges = 0;
    m_SpareEdges = m_SpareVertices = 0;
    }

  /** Destructor */
  ~BinaryImageToGraphFilter()
    {
    if(m_AdjacencyIndex) 
      {
      delete m_AdjacencyIndex;
      delete m_Adjacency;
      delete m_ImageIndex;
      }
    }

  /** Set the input */
  void SetInput(TImage *image) { this->SetNthInput(0,image); }

  /** Request a number of 'empty' vertices to be
    allocated in the vertex array */
  itkSetMacro(SpareVertices, unsigned int);
  itkGetMacro(SpareVertices, unsigned int);

  /** Request a number of 'empty' directional edges to be allocated in 
    the edges array */
  itkSetMacro(SpareEdges, unsigned int);
  itkGetMacro(SpareEdges, unsigned int);

  /** Update method (why?) */
  void Update() { this->GenerateData(); }

  /** Get the adjacency index */
  VertexType *GetAdjacencyIndex() { return m_AdjacencyIndex; }

  /** Get the adjacency index */
  VertexType *GetAdjacencyArray() { return m_Adjacency; }

  /** Get the number of vertices */
  unsigned int GetNumberOfVertices() { return m_NumberOfVertices; }

  /** Get the number of directional edges (2x symmetric edges) */
  unsigned int GetNumberOfEdges() { return m_NumberOfEdges; }

  /** Get number of vertex's adjacencies */
  unsigned int GetVertexNumberOfNeighbors(unsigned int iVertex)
    {
    return m_AdjacencyIndex[iVertex+1] - m_AdjacencyIndex[iVertex];
    }
  
  /** Get the adjacency list for a vertex */
  VertexType *GetVertexNeighbors(unsigned int iVertex)
    {
    return m_Adjacency + m_AdjacencyIndex[iVertex];
    }

  /** Get the image index associated with a vertex */
  IndexType GetVertexImageIndex(unsigned int iVertex) 
    {
    return m_ImageIndex[iVertex];
    }

  /** Generate data */
  void GenerateData() 
    {
    // Get the image
    typename ImageType::ConstPointer image =
      reinterpret_cast<ImageType *>(this->GetInput(0));
    
    // Typedefs 
    typedef typename itk::ConstantBoundaryCondition<ImageType> ConditionType;
    typedef typename itk::ConstNeighborhoodIterator<ImageType, ConditionType> 
    	IteratorType;

    // Initialize the iterator
    itk::Size<ImageDimension> radius;
    radius.Fill(1);
    IteratorType it(radius, image, image->GetBufferedRegion());

    // The vector of neighboring pixels
    unsigned int iCenter = it.Size() >> 1;
    unsigned int iIndex[ImageDimension * 2];
    for(unsigned int iDim=0;iDim<ImageDimension;iDim++)
      {
      iIndex[iDim << 1] = iCenter - it.GetStride(iDim);
      iIndex[(iDim << 1) + 1] = iCenter + it.GetStride(iDim);
      }

    // Create an image that maps connected pixels to graph index
    typedef itk::Image<int, ImageDimension> IndexImage;
    typename IndexImage::Pointer imgIndex = IndexImage::New();
    imgIndex->SetRegions(image->GetBufferedRegion());
    imgIndex->Allocate();
    imgIndex->FillBuffer(-1);

    // We will count edges and vertices
    m_NumberOfVertices = 0;
    m_NumberOfEdges = 0;

    // Continue until iterator finishes
    while(!it.IsAtEnd())
      {
      if(it.GetCenterPixel() != 0) 
        {
        // See if the pixel has any non-zero vertices
        unsigned int nNeighbors = 0;
        for(unsigned int i=0;i<(ImageDimension << 1);i++)
          if(it.GetPixel(iIndex[i]) != 0) 
            ++nNeighbors;

        if(nNeighbors > 0)
          {
          m_NumberOfEdges += nNeighbors;
          imgIndex->SetPixel(it.GetIndex(iCenter), m_NumberOfVertices);
          m_NumberOfVertices++;
          }
        }
      ++it;
      }

    // Clear the arrays if needed
    if(m_AdjacencyIndex) 
      {
      delete m_AdjacencyIndex;
      delete m_Adjacency;
      delete m_ImageIndex;
      }
      
    // Allocate the arrays 
    m_AdjacencyIndex = new VertexType[m_NumberOfVertices + m_SpareVertices + 1];
    m_Adjacency = new VertexType[m_NumberOfEdges + m_SpareEdges];
    m_ImageIndex = new IndexType[m_NumberOfVertices];

    // Now, create the adjacency structure
    it.GoToBegin();
    unsigned int iVertex = 0, iEdge = 0;
    while(!it.IsAtEnd())
      {
      int offset = imgIndex->GetPixel(it.GetIndex(iCenter));
      if(offset >= 0)
        {
        // Pixel should be included
        m_ImageIndex[iVertex] = it.GetIndex(iCenter);
        m_AdjacencyIndex[iVertex] = iEdge;
        ++iVertex;

        // Add all the edges of the pixel 
        for(unsigned int i=0;i<(ImageDimension << 1);i++)
          {
          if(it.GetPixel(iIndex[i]) != 0)
            {
            // cout << it.GetIndex(iIndex[i]) << endl;
            m_Adjacency[iEdge++] = imgIndex->GetPixel(it.GetIndex(iIndex[i]));
        
            }
          }
        }
      ++it;
      }

    // Complete xAdjIndex
    m_AdjacencyIndex[iVertex] = iEdge;
    }

protected:

  /** Mapping from vertices to indices in the adjacency array */
  VertexType *m_AdjacencyIndex;

  /** List of adjacencies indexed by the above array */
  VertexType *m_Adjacency;

  /** Mapping from vertices to input image indices */
  IndexType *m_ImageIndex;

  /** The number of directed edges (2x undirected). These do not include 
   the spare vertices and edges */
  unsigned int m_NumberOfVertices, m_NumberOfEdges;

  /** Spare edges and vertices */
  unsigned int m_SpareVertices, m_SpareEdges;
};

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

#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_powell.h>

class MetisPartitionProblem : public SingleValuedCostFunction
{
public:
  typedef MetisPartitionProblem Self;
  typedef SingleValuedCostFunction Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  itkTypeMacro(MetisPartitionProblem, SingleValuedCostFunction);

  itkNewMacro(Self);
  
  typedef itk::Image<short,3> ImageType; 
  typedef BinaryImageToGraphFilter<ImageType> GraphFilter;

  typedef Superclass::MeasureType MeasureType;
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::DerivativeType DerivativeType;

  /** Set the problem parameters */
  void SetProblem(GraphFilter *fltGraph, int *wgtEdges, unsigned int nParts);

  /** Return the number of parameters */
  unsigned int GetNumberOfParameters() const
    { return m_NumberOfParameters; }
  
  /** Virtual method from parent class */
  MeasureType GetValue(const ParametersType &x) const;

  /** Not used, since there are no derivatives to evaluate */
  void GetDerivative(const ParametersType &, DerivativeType &) const {}

  /** Get the result of the last partition */
  idxtype *GetLastPartition() { return m_Partition; }
  
private:
  /** The stored graph information */
  GraphFilter *m_Graph;

  /** Edge weights */
  int *m_EdgeWeights;
  
  /** The partition array */
  idxtype *m_Partition;

  /** Problem size */
  unsigned int m_NumberOfParameters;
};

void
MetisPartitionProblem
::SetProblem(GraphFilter *fltGraph, int *wgtEdges, unsigned int nParts)
{
  m_Graph = fltGraph;
  m_EdgeWeights = wgtEdges;
  m_Partition = new idxtype[m_Graph->GetNumberOfVertices()];
  m_NumberOfParameters = nParts;
}

MetisPartitionProblem::MeasureType
MetisPartitionProblem
::GetValue(const ParametersType &x) const 
{
  // Convert the vector to floating array
  vnl_vector<float> wgt(x.size() + 1);
  float sum = 0;
  for(unsigned int i=0;i<x.size();i++)
    {
    wgt[i] = (float) x[i];
    sum += wgt[i];
    }
  wgt[x.size()] = 1.0f - sum;

  // Variables used to call METIS
  int nVertices = m_Graph->GetNumberOfVertices();
  int nParts = x.size() + 1;
  int wgtflag = 1;
  int numflag = 0;
  int options[] = {0,0,0,0,0};
  int edgecut = 0;

  // Apply METIS to the graph
  cout << "applying METIS with weights " << wgt << flush;

  METIS_WPartGraphRecursive(
    &nVertices,
    m_Graph->GetAdjacencyIndex(),
    m_Graph->GetAdjacencyArray(),
    NULL,
    m_EdgeWeights,
    &wgtflag,
    &numflag,
    &nParts,
    wgt.data_block(),
    options,
    &edgecut,
    m_Partition);

  // Report the edge cut
  cout << "  - done - edge cut " << edgecut << endl;

  return (MeasureType) edgecut;
}

int main(int argc, char *argv[])
{
  // Check arguments
  if(argc < 4) return usage();

  // Variables to hold command line arguments
  string fnInput(argv[argc-3]), fnOutput(argv[argc-2]);
  int nParts = atoi(argv[argc-1]);
  vnl_vector<double> xWeights(nParts,1.0 / nParts);
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

  typedef itk::Image<short,3> ImageType;
  ImageType::Pointer img = ImageType::New();
  ReadImage(img,fnInput.c_str());

  cout << "   image has dimensions " << img->GetBufferedRegion().GetSize() 
    << ", nPixels = " << img->GetBufferedRegion().GetNumberOfPixels() << endl;

  // Convert the image to a graph
  cout << "converting image to graph structure" << endl;

  typedef BinaryImageToGraphFilter<ImageType,idxtype> GraphFilter;
  GraphFilter::Pointer fltGraph = GraphFilter::New();
  fltGraph->SetInput(img);
  fltGraph->Update();
  
  // Get the number of vertices and edges
  int nVertices = fltGraph->GetNumberOfVertices();
  int nEdges = fltGraph->GetNumberOfEdges();
  int iVertex, iEdge = 0;
  
  // Compute weight arrays
  GraphFilter::VertexType *xEdgeWeight = new GraphFilter::VertexType [nEdges];
  for(iVertex = 0;iVertex < nVertices;iVertex++)
  	{
  	// Get the index of the vertex of interest
  	GraphFilter::IndexType idx = fltGraph->GetVertexImageIndex(iVertex);
  	
  	// Look at all the edges for this vertex
  	unsigned int nNbr = fltGraph->GetVertexNumberOfNeighbors(iVertex);
  	for(unsigned int iNbr=0; iNbr < nNbr; iNbr++)
  		{
  		// Get the neighbor vertex
  		GraphFilter::VertexType n = fltGraph->GetVertexNeighbors(iVertex)[iNbr];
  		GraphFilter::IndexType nidx = fltGraph->GetVertexImageIndex(n);	
  			
  		// Check if the vertex is on the edge
  		if(iPlaneDim >= 0 && (
  			(idx[iPlaneDim] == iPlaneSlice && nidx[iPlaneDim] == iPlaneSlice+1) ||
  			(nidx[iPlaneDim] == iPlaneSlice && idx[iPlaneDim] == iPlaneSlice+1)))
  			{
        // cout << "edge " << iEdge << " from " << idx << " to " << nidx << endl;
  			xEdgeWeight[iEdge++] = 1;  			  				
  			}
  		else
  		  {
  		  xEdgeWeight[iEdge++] = iPlaneStrength;  			  					
  		  }	
  		}
  	}

  // Create a METIS problem based on the graph and weights
  MetisPartitionProblem::Pointer mp = MetisPartitionProblem::New();
  mp->SetProblem(fltGraph,xEdgeWeight,nParts - 1);

  // Get the starting solution
  MetisPartitionProblem::ParametersType x(nParts - 1);
  x.set(xWeights.data_block());

  if(flagOptimize)
    {
    typedef OnePlusOneEvolutionaryOptimizer Optimizer;
    Optimizer::Pointer opt = Optimizer::New();

    typedef itk::Statistics::NormalVariateGenerator Generator;
    Generator::Pointer generator = Generator::New();

    opt->SetCostFunction(mp);
    opt->SetInitialPosition(x);
    opt->SetInitialRadius(0.001);
    opt->SetMaximumIteration(100);
    opt->SetNormalVariateGenerator(generator);
    opt->StartOptimization();

    x = opt->GetCurrentPosition();
    }

  // Evaluate the problem using the best value
  cout << "   edge cut value is " << mp->GetValue(x) << endl;
  idxtype *xPartition = mp->GetLastPartition();

  // Verify the edge cut ourselves
  /*
  iEdge = 0;
  int cut = 0;
  for(iVertex = 0;iVertex < nVertices;iVertex++)
  	{
  	// Look at all the edges for this vertex
  	unsigned int nNbr = fltGraph->GetVertexNumberOfNeighbors(iVertex);
  	for(unsigned int iNbr=0; iNbr < nNbr; iNbr++)
  		{
  		// Get the neighbor vertex
  		GraphFilter::VertexType n = fltGraph->GetVertexNeighbors(iVertex)[iNbr];
  
      if(xPartition[iVertex] != xPartition[n])
        cut += xEdgeWeight[fltGraph->GetAdjacencyIndex()[iVertex] + iNbr];
      }
  	}
    cout << "verified cut is " << (cut/2) << endl;
  */

  // Create output image 
  ImageType::Pointer imgOut = ImageType::New();
  imgOut->SetRegions(img->GetBufferedRegion());
  imgOut->Allocate();
  imgOut->FillBuffer(0);
    
  // Apply partition to output image
  for(unsigned int iVertex = 0; iVertex < nVertices; iVertex++)
    {
    GraphFilter::IndexType idx = fltGraph->GetVertexImageIndex(iVertex);
    imgOut->SetPixel(idx,xPartition[iVertex] + 1);
    }
    
  // Write the image
  cout << "writing output image" << endl;

  WriteImage(imgOut,fnOutput.c_str());

  return 0;
}
