#ifndef __ImageToGraphFilter_h_
#define __ImageToGraphFilter_h_

#include <itkImage.h>
#include <itkConstantBoundaryCondition.h>
#include <itkConstNeighborhoodIterator.h>

/* ***************************************************************************
 * FUNCTOR DEFINITIONS
 * *************************************************************************** */

/** 
 * This is the base class for functors that can be used to tell 
 * the ImageToGraphFilter how to compute edge and vertex weights
 */
template < class TImage, class TWeight = int >
class AbstractGraphWeightFunctor
{
public:
  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  typedef typename TImage::PixelType TPixel;
  typedef itk::Point<float, ImageDimension> Point;

  /** Return true if the vertex should be included in the graph */
  virtual bool IsPixelAVertex( TPixel i, Point x )
    { return IsPixelAVertex(i); }

  /** Compute the edge weight between two pixels */
  virtual TWeight GetEdgeWeight(
    TPixel i1, const Point &x1, TPixel i2, const Point &x2)
    { return GetEdgeWeight( i1, i2 ); }
  
  /** Compute the vertex weight at a pixel */
  virtual TWeight GetVertexWeight(TPixel i, const Point &x)
    { return GetVertexWeight(i); }

protected:

  /** Return true if the vertex should be included in the graph */
  virtual bool IsPixelAVertex( TPixel i ) = 0;

  /** Compute the edge weight between two pixels */
  virtual TWeight GetEdgeWeight(TPixel i1, TPixel i2) = 0;
  
  /** Compute the vertex weight at a pixel */
  virtual TWeight GetVertexWeight(TPixel i) = 0;
};

/**
 * This is the standard weight functor that is applicable for binary
 * images. It considers pixels with non-zero intensities to be vertices
 * and assigns unit weights to edges and vertices 
 */
template <class TImage, class TWeight = int>
class BinaryGraphWeightFunctor 
: public AbstractGraphWeightFunctor<TImage, TWeight>
{
public:
  typedef AbstractGraphWeightFunctor<TImage, TWeight> Superclass;
  typedef typename Superclass::TPixel TPixel;

protected:
  /** Vertices with zero intensity are excluded */
  virtual bool IsPixelAVertex( TPixel i )
    { return (i != 0); }

  /** Compute the edge weight between two pixels */
  virtual TWeight GetEdgeWeight(TPixel i1, TPixel i2)
    { return 1; }
  
  /** Compute the vertex weight at a pixel */
  virtual TWeight GetVertexWeight(TPixel i)   
    { return 1; }
};


/* ***************************************************************************
 * MAIN CLASS DEFINITION
 * *************************************************************************** */

/**
 * \class ImageToGraphFilter
 * \brief Constructs a graph from a 3D image
 *
 * This class is intended to be used with the METIS and CLUTO graph 
 * partitioning libraries, but it may have other uses as well. Voxels
 * in a 3D image are converted to vertices in a graph, and pairs of
 * neighboring vertices are converted to edges.
 *
 * The user can specify a functor derived from AbstractGraphWeightFunctor
 * that will be used to compute the weights of the edges and vertices
 * in the graph.
 */
template<class TImage, class TVertex = int, class TWeight = int> 
class ImageToGraphFilter : public itk::ProcessObject
{
public:
  typedef ImageToGraphFilter Self;
  typedef itk::ProcessObject Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef TImage ImageType;
  typedef typename TImage::IndexType IndexType;
  typedef TVertex VertexType;
  typedef AbstractGraphWeightFunctor<TImage, TWeight> WeightFunctorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);
  
  /** Constructor */
  ImageToGraphFilter()
    {
    m_Adjacency = NULL;
    m_AdjacencyIndex = NULL;
    m_ImageIndex = NULL;
    m_VertexWeights = NULL;
    m_EdgeWeights = NULL;

    m_NumberOfVertices = 0;
    m_NumberOfEdges = 0;
    m_SpareEdges = m_SpareVertices = 0;
    m_WeightFunctor = &m_DefaultWeightFunctor;
    }

  /** Destructor */
  ~ImageToGraphFilter()
    {
    CleanUpGraphArrays();
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

  /** Get and set the weight table */
  // itkSetMacro(WeightFunctor, WeightFunctorType *);
  itkGetMacro(WeightFunctor, WeightFunctorType *);

  void SetWeightFunctor(WeightFunctorType *in_Functor)
    { m_WeightFunctor = in_Functor; }

  /** Update method (why?) */
  void Update() { this->GenerateData(); }

  /** Get the adjacency index */
  itkGetMacro( AdjacencyIndex, VertexType * );

  /** Get the adjacency index */
  itkGetMacro( Adjacency, VertexType * );

  /** Get the array of vertex weights */
  itkGetMacro( VertexWeights, TWeight * );

  /** Get the array of edge weights */
  itkGetMacro( EdgeWeights, TWeight * );

  /** Get the number of vertices */
  itkGetMacro( NumberOfVertices, unsigned int );

  /** Get the number of directional edges (2x symmetric edges) */
  itkGetMacro( NumberOfEdges, unsigned int );

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
    
    // Initialize the iterator
    itk::Size<ImageDimension> radius;
    radius.Fill(1);
    NeighborhoodIterator it(radius, image, image->GetBufferedRegion());

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
      // Only consider the vertex if it has non-zero weight
      if(IsPixelAVertex(it, iCenter)) 
        {
        // See if any of the pixel's neighbors qualify
        unsigned int nNeighbors = 0;
        for(unsigned int i=0;i<(ImageDimension << 1);i++)
          {
          // If the neighbor has a positive weight, increment the number of nbrs
          if(IsPixelAVertex(it, iIndex[i])) 
            ++nNeighbors;
          }

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
    CleanUpGraphArrays();
      
    // Allocate the arrays 
    m_AdjacencyIndex = new VertexType[m_NumberOfVertices + m_SpareVertices + 1];
    m_Adjacency = new VertexType[m_NumberOfEdges + m_SpareEdges];
    m_ImageIndex = new IndexType[m_NumberOfVertices];
    m_VertexWeights = new TWeight[m_NumberOfVertices + m_SpareVertices];
    m_EdgeWeights = new TWeight[m_NumberOfEdges + m_SpareEdges];

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

        // Compute the weight of the vertex 
        m_VertexWeights[iVertex++] = GetVertexWeight(it, iCenter);

        // Add all the edges of the pixel 
        for(unsigned int i=0;i<(ImageDimension << 1);i++)
          {
          if(it.GetPixel(iIndex[i]) != 0)
            {
            // Add the edge to the adjacency list
            m_Adjacency[iEdge] = imgIndex->GetPixel(it.GetIndex(iIndex[i]));

            // Compute the weight of the edge
            m_EdgeWeights[iEdge++] = GetEdgeWeight(it, iCenter, iIndex[i]);
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

  /** List of weights associated with each vertex */
  TWeight *m_VertexWeights;

  /** List of weights associated with each adjacency (directed edge) */
  TWeight *m_EdgeWeights;
  
  /** A table that determines the weights for edges and vertices */
  WeightFunctorType *m_WeightFunctor;

  /** The default weight table */
  BinaryGraphWeightFunctor<TImage, TWeight> m_DefaultWeightFunctor;

  /** The number of directed edges (2x undirected). These do not include 
   the spare vertices and edges */
  unsigned int m_NumberOfVertices, m_NumberOfEdges;

  /** Spare edges and vertices */
  unsigned int m_SpareVertices, m_SpareEdges;

  // Neighborhood iterator typedefs
  typedef typename itk::ConstantBoundaryCondition<ImageType> ConditionType;
  typedef typename itk::ConstNeighborhoodIterator<ImageType, ConditionType> 
    NeighborhoodIterator;

  /** Get the weight associated with a vertex (internal method) */
  int IsPixelAVertex( const NeighborhoodIterator &it, int index )
    {
    // Get the position of the vertex in image coordinates
    WeightFunctorType::Point xVertex;
    it.GetImagePointer()->TransformIndexToPhysicalPoint(
      it.GetIndex(index), xVertex);

    // Get the weight of the vertex from the weight table
    return m_WeightFunctor->IsPixelAVertex( it.GetPixel(index), xVertex );
    }

  /** Get the weight associated with a vertex (internal method) */
  int GetVertexWeight( const NeighborhoodIterator &it, int index )
    {
    // Get the position of the vertex in image coordinates
    WeightFunctorType::Point xVertex;
    it.GetImagePointer()->TransformIndexToPhysicalPoint(
      it.GetIndex(index), xVertex);

    // Get the weight of the vertex from the weight table
    return m_WeightFunctor->GetVertexWeight( it.GetPixel(index), xVertex );
    }

  /** Get the weight associated with an edge (internal method) */
  int GetEdgeWeight( const NeighborhoodIterator &it, int i1, int i2 )
    {
    // Get the position of the vertex in image coordinates
    WeightFunctorType::Point xVertex1, xVertex2;
    it.GetImagePointer()->TransformIndexToPhysicalPoint( it.GetIndex(i1), xVertex1);
    it.GetImagePointer()->TransformIndexToPhysicalPoint( it.GetIndex(i2), xVertex2);

    // Get the weight of the vertex from the weight table
    return m_WeightFunctor->GetEdgeWeight( 
      it.GetPixel(i1), xVertex1, it.GetPixel(i2), xVertex2 );
    }

  /** Memory cleanup for graph arrays */
  void CleanUpGraphArrays()
    {
    if(m_AdjacencyIndex) 
      {
      delete m_AdjacencyIndex;
      delete m_Adjacency;
      delete m_ImageIndex;
      delete m_EdgeWeights;
      delete m_VertexWeights;
      }
    }
};

#endif // __ImageToGraphFilter_h_
