/**
 * Read metis output, create an image
 */
#include "ReadWriteImage.h"
#include <iostream>

// #include "metis.h"
typedef int idxtype;
extern "C" {
  void METIS_WPartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, 
                       int *, int *, float *, int *, int *, idxtype *);

  void METIS_WPartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, 
                       int *, int *, float *, int *, int *, idxtype *);
}

using namespace std;

int usage()
{
  cout << "usage: metisseg [options] input.img output.img num_part" << endl;
  cout << "   uses METIS to segment a binary image into num_part partitions" << endl;
  cout << "options: " << endl;
  cout << "   -w N XX.X           weight of partition N is XX.X" << endl;
  cout << "   -s N                number of 'slack nodes'" << endl;
  cout << "   -p N1 N2            define cut plane at dimension N1, slice N2" << endl;
  return -1;
}

#include <itkConstantBoundaryCondition.h>
#include <itkConstNeighborhoodIterator.h>

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

  /** Update method (why?) */
  void Update() { this->GenerateData(); }

  /** Get the adjacency index */
  VertexType *GetAdjacencyIndex() { return m_AdjacencyIndex; }

  /** Get the adjacency index */
  VertexType *GetAdjacencyArray() { return m_Adjacency; }

  /** Get the number of vertices */
  unsigned int GetNumberOfVertices() { return m_NumberOfVertices; }

  /** Get the number of edges */
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
    typedef typename itk::ConstNeighborhoodIterator<ImageType, ConditionType> IteratorType;

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
    m_AdjacencyIndex = new VertexType[m_NumberOfVertices + 1];
    m_Adjacency = new VertexType[m_NumberOfEdges];
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

  /** The number of directed edges (2x undirected) */
  unsigned int m_NumberOfVertices, m_NumberOfEdges;
};

int main(int argc, char *argv[])
{
  // Check arguments
  if(argc < 4) return usage();

  // Variables to hold command line arguments
  string fnInput(argv[argc-3]), fnOutput(argv[argc-2]);
  int nParts = atoi(argv[argc-1]);
  vnl_vector<float> xWeights(nParts,1.0 / nParts);
  unsigned int nSlack = 0;
  int iPlaneDim = -1, iPlaneSlice = -1;
  
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
    else if(!strcmp(argv[iArg],"-s"))
      {
      nSlack = atoi(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg],"-p"))
      {
      iPlaneDim = atoi(argv[++iArg]);
      iPlaneSlice = atoi(argv[++iArg]);
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
    cout << "will insert cut plane at slice " << iPlaneSlice << " in dimension " << iPlaneDim << endl;
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

  typedef BinaryImageToGraphFilter<ImageType> GraphFilter;
  GraphFilter::Pointer fltGraph = GraphFilter::New();
  fltGraph->SetInput(img);
  fltGraph->Update();

  cout << "   graph has " << fltGraph->GetNumberOfVertices() << " vertices and " 
    << fltGraph->GetNumberOfEdges() << " edges" << endl;

  // Verify graph
  cout << "verifying graph" << endl;
  for(unsigned int iVertex=0;iVertex<fltGraph->GetNumberOfVertices();iVertex++)
    {
    unsigned int nAdj = fltGraph->GetVertexNumberOfNeighbors(iVertex);
    GraphFilter::VertexType *adj = fltGraph->GetVertexNeighbors(iVertex);
    GraphFilter::IndexType idxCenter = fltGraph->GetVertexImageIndex(iVertex);
    for(unsigned int iAdj=0;iAdj<nAdj;iAdj++)
      {
      GraphFilter::IndexType idxAdj = fltGraph->GetVertexImageIndex(adj[iAdj]);
      unsigned int diff = 0;
      for(unsigned int d=0;d<3;d++)
        diff += abs(idxCenter[d] - idxAdj[d]);
      if(diff != 1)
        cout << "   error at node " << iVertex << endl;
      }
    }
  
  // Compute weight arrays
  

  // Apply METIS to the graph
  int nVertices = fltGraph->GetNumberOfVertices();
  int wgtflag = 0;
  int numflag = 0;
  int options[] = {0,0,0,0,0};
  int edgecut = 0;
  idxtype *xPartition = new idxtype[nVertices];

  cout << "running METIS grapg partition algorithm" << endl;

  METIS_WPartGraphRecursive(
    &nVertices,
    (idxtype *) fltGraph->GetAdjacencyIndex(),
    (idxtype *) fltGraph->GetAdjacencyArray(),
    NULL,
    NULL,
    &wgtflag,
    &numflag,
    &nParts,
    xWeights.data_block(),
    options,
    &edgecut,
    (idxtype *)xPartition);

  cout << "   edge cut value is " << edgecut << endl;

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
