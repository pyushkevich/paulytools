#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"

#include "Danielsson.h"
#include "DrawTriangles.h"

#include "matrix.h"
#include "EvolutionaryStrategy.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

using namespace itk;
using namespace std;
using namespace pauly;

typedef vnl_vector_fixed<double,3> Vec;
typedef pauly::Vector OptVec;

struct CurvePoint {
  // A point painted manually on the grey matter
  Vec xGrey;
    
  // A matching point generated on the white matter
  Vec xWhite;
};

typedef vector<CurvePoint> CurveType;
typedef Image<char,3> CharImageType;
typedef Image<short,3> ShortImageType;
typedef ImageFileReader<CharImageType> CharReaderType;
typedef ImageFileReader<ShortImageType> ShortReaderType;
typedef ImageFileWriter<ShortImageType> ShortWriterType;
typedef ImageFileWriter<CharImageType> CharWriterType;

double AreaOfTriangle(double a, double b, double c)
{
  double s = 0.5 * (a + b + c);
  double area = sqrt(s * (s-a) * (s-b) * (s-c));
  return area; 
}

// A problem associated with optimizing grey matter partition ribbons
class RibbonProblem : public Function
{
public:
  // Construct a problem
  RibbonProblem(const CurveType &curve, CharImageType *imgWhite,
    unsigned int iStart, unsigned int iLength)
    {
    m_Curve = curve;
    m_Start = iStart;
    m_Length = iLength;
    m_Image = imgWhite;
    }

  // Generate a vector corresponding to the current ribbon state
  OptVec GetVectorRepresentation() 
    {
    OptVec v(m_Length * 3);
    unsigned int j = 0;
    for(unsigned int i = 0;i < m_Length;i++)
      {
      v(j++) = m_Curve[i + m_Start].xWhite[0];
      v(j++) = m_Curve[i + m_Start].xWhite[1];
      v(j++) = m_Curve[i + m_Start].xWhite[2];
      }
    return v;
    }

  // Get upper and lower bounds
  void ComputeVectorBounds(OptVec &upper, OptVec &lower)
    {
    upper.setSize(m_Length * 3);
    lower.setSize(m_Length * 3);
    lower.setAll(1);

    unsigned int j = 0;
    for(unsigned int i = 0;i < m_Length;i++)
      {
      upper(j++) = m_Image->GetBufferedRegion().GetSize(0) - 2;
      upper(j++) = m_Image->GetBufferedRegion().GetSize(1) - 2;
      upper(j++) = m_Image->GetBufferedRegion().GetSize(2) - 2;
      }
    }

  // Update the current ribbon state with a vector
  CurveType &UpdateRibbon(const OptVec &v)
    {
    assert(v.size() == 3 * m_Length);
    unsigned int j = 0;
    for(unsigned int i = 0;i < m_Length;i++)
      {
      m_Curve[i + m_Start].xWhite[0] = v(j++);
      m_Curve[i + m_Start].xWhite[1] = v(j++);
      m_Curve[i + m_Start].xWhite[2] = v(j++);
      }
    return m_Curve;
    }

  void ComputeObjectives(double &area, double &greyRatio)
    {
    // Get a raw pointer to the data as well as the buffer access steps
    char *rawimg = m_Image->GetBufferPointer();
    unsigned int nLine = m_Image->GetBufferedRegion().GetSize(0); 
    unsigned int nSlice = m_Image->GetBufferedRegion().GetSize(1) * nLine;
    
    // Measure the costs associated with the ribbon
    double xArea = 0.0;
    double xWidth = 0.0;
    double xLenWhite = 0.0;
    double xGreyCount = 0.0;
    
    for(unsigned int i=0;i<m_Curve.size() - 1;i++)
      {
      // Get the two points
      CurvePoint &p1 = m_Curve[i];
      CurvePoint &p2 = m_Curve[i+1];

      // Compute the edges of the ribbon
      double w1 = (p1.xGrey - p1.xWhite).two_norm();
      double w2 = (p2.xGrey - p2.xWhite).two_norm();
      double lg = (p1.xGrey - p2.xGrey).two_norm();
      double lw = (p1.xWhite - p2.xWhite).two_norm();
        
      // Width, Length integration 
      xWidth += 0.5 * (w2 - w1); 
      xLenWhite += lw;

      // Area integration
      double h = (p2.xWhite-p1.xGrey).two_norm(); 
      xArea += AreaOfTriangle(w1,lw,h);
      xArea += AreaOfTriangle(w2,lg,h);

      // Generate a list of voxels that the line passes through
      for(double t=0;t<1;t+=0.05)
        {
        // Compute the corresponding image voxel
        int x = (int) round( p1.xWhite[0] * (1-t) + p2.xWhite[0] * t );
        int y = (int) round( p1.xWhite[1] * (1-t) + p2.xWhite[1] * t );
        int z = (int) round( p1.xWhite[2] * (1-t) + p2.xWhite[2] * t );

        // Get the offset value
        int offset = z * nSlice + y * nLine + x;
        double pixel = rawimg[offset] == 0 ? 1.0 : 0.0;
       
        // Index<3> index = {{x,y,z}};
        // assert(rawimg[offset] == m_Image->GetPixel(index));
        
        // Find the intensity at the voxel
        xGreyCount += 0.05 * lw * pixel;
        }
      }

    area = xArea;
    greyRatio = xGreyCount / xLenWhite;
    }
  
  double evaluate(const OptVec &x)
    {
    // Update the ribbon with the parameter vector
    UpdateRibbon(x);

    // Compute the components of the objective
    double xArea, xGreyRatio;
    ComputeObjectives(xArea,xGreyRatio);

    // Compute overall objective value based on the two components
    double obj = 10000 * xGreyRatio + xArea;

    // Increase the iterations counter
    setEvaluationCost(getEvaluationCost() + 1);

    // Return the objective
    return obj;
    }

  // Generate the initial search space given the current state
  void InitializeSolutionSpace(SolutionSpace **space) 
    {
    // Create a vector for the current state
    OptVec mean = GetVectorRepresentation();

    // Generate a standard deviation vector
    OptVec sd(mean.size(),2);
   
    // Return a Gaussian solution space
    *space = new GaussianSS(mean,sd);
    }

private:
  // The curve undergoing optimization
  CurveType m_Curve;

  // The range of optimization
  unsigned int m_Start, m_Length;

  // The white matter image
  CharImageType::Pointer m_Image;
};

int main(int argc, char *argv[])
{
  // This program projects curves painted on the surface of the 
  // grey matter into the white matter using distance transforms
  //
  // The parameters to this program are:
  //    Grey Matter Image
  //    Write Matter Image
  //    Curve File
  //
  // The output from this program is:
  //    Grey Matter Image with curves marked with a zero intensity
  //    and projecting through the grey matter

  // Images that we will be using
  CharImageType::Pointer imgWhite = NULL;
  ShortImageType::Pointer imgGrey = NULL;

  // Check parameters
  if(argc != 4)
    {
    cout << "Usage: bpaint greyimg whiteimg curvefile" << endl;
    return -1;
    }

  try 
    {
    // Read in the grey matter image
    ShortReaderType::Pointer readerGrey = ShortReaderType::New();
    readerGrey->SetFileName(argv[1]);
    readerGrey->Update();
    imgGrey = readerGrey->GetOutput();

    // Read in the white matter image
    CharReaderType::Pointer readerWhite = CharReaderType::New();
    readerWhite->SetFileName(argv[2]);
    readerWhite->Update();
    imgWhite = readerWhite->GetOutput();
    }
  catch(...)
    {
    cout << "Image read error" << endl;
    return -1;
    }

  // Allocate an array for curves
  typedef vector<CurvePoint> CurveType;
  vector<CurveType> vCurves;

  // Read in the curves on the surface
  FILE *fCurves = fopen(argv[3],"rt");
  char buffer[1024];
  while(!feof(fCurves)) 
    {
    // Read a line from the curve file
    fgets(buffer,1024,fCurves);
    string line(buffer);

    // Check the first 6 characters
    if(line.substr(0,6) == "NODES:")
      {
      CurveType c;
      vCurves.push_back(c);
      }
    else if (line.substr(0,6) == "      ")
      {
      CurvePoint p;
      istringstream iss(buffer);
      iss >> p.xGrey[0];
      iss >> p.xGrey[1];
      iss >> p.xGrey[2];

      vCurves.back().push_back(p);
      }
    }

  // Go through the points and trim away some

  // At this point we've loaded the input into memory. We now compute the 
  // vector distance transform of the white matter

  // Create a copy of the white image
  char *dataWhite = new char[imgWhite->GetBufferedRegion().GetNumberOfPixels()];
  memcpy(dataWhite,imgWhite->GetBufferPointer(),imgWhite->GetBufferedRegion().GetNumberOfPixels());
  
  short *xDistance[3];
  edt3ddan(dataWhite,
    (int) imgWhite->GetBufferedRegion().GetSize(0),
    (int) imgWhite->GetBufferedRegion().GetSize(1),
    (int) imgWhite->GetBufferedRegion().GetSize(2),
    6,xDistance,xDistance+1,xDistance+2);

  delete dataWhite;

  // Create the distance transform images
  typedef Image<short,3> OffsetImageType;
  OffsetImageType::Pointer imgOffset[3];
  for(unsigned int i=0;i<3;i++)
    {
    imgOffset[i] = OffsetImageType::New();
    imgOffset[i]->SetRegions(imgWhite->GetBufferedRegion());
    imgOffset[i]->Allocate();
    memcpy(imgOffset[i]->GetBufferPointer(),xDistance[i],
      imgOffset[i]->GetBufferedRegion().GetNumberOfPixels() * sizeof(short));
    free(xDistance[i]);
    }
  
  // Find the nearest white matter pixel for each of the curves
  unsigned int iCurve, iPoint;
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];
    for(iPoint=0;iPoint<curve.size();iPoint++)
      {
      CurvePoint &p = curve[iPoint];
      Index<3> idx = {{
        (unsigned int)p.xGrey[0],
        (unsigned int)p.xGrey[1],
        (unsigned int)p.xGrey[2]}};

      for(unsigned int j=0;j<3;j++)
        p.xWhite[j] = p.xGrey[j] + imgOffset[j]->GetPixel(idx);

      }
    }

  
  // Save the distance transform
  ImageFileWriter<OffsetImageType>::Pointer wtmp = ImageFileWriter<OffsetImageType>::New();
  wtmp->SetInput(imgOffset[0]);
  wtmp->SetFileName("whitedist.hdr");
  wtmp->Update();

  // This is the most interesting part of the program. We attempt to 
  // improve the ribbons by making their free side pass entirely through
  // the white matter
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    // Create a curve
    CurveType &curve = vCurves[iCurve];

    // Create an overall problem
    RibbonProblem pbTotal(curve,imgWhite,0,curve.size());

    // Compute the objective values
    double xArea, xGreyRatio;
    pbTotal.ComputeObjectives(xArea,xGreyRatio);
    cout << "Ribbon " << iCurve << " of length " << curve.size() << " has area " << xArea 
      << " and non-white ratio " << xGreyRatio << endl;

    // Create a list of starting points
    vector<unsigned int> vStarts;
    unsigned int nSpan = 4;
    for(unsigned int iStart=0;iStart < curve.size()-nSpan;iStart+=1)
      vStarts.push_back(iStart);
    
    // Iteratevely optimize over ribbon segments
    for(unsigned int iIteration = 0;iIteration < 40;iIteration++)
      {
      // State the iteration
      cout << "   Running iteration " << iIteration << endl;
      
      // Perform a random permutation of the starting points
      random_shuffle(vStarts.begin(),vStarts.end());

      // Peform optimizations
      for(unsigned int iStep=0;iStep<vStarts.size();iStep++)
        {
        // Create a Ribbon problem
        RibbonProblem pbChunk(curve,imgWhite,vStarts[iStep],nSpan);

        // Create a corresponding solution space
        SolutionSpace *ss;
        pbChunk.InitializeSolutionSpace(&ss);
        
        // Create an evolutionary strategy
        EvolutionaryStrategy es(pbChunk,*ss,2,4,SELECTION_MuPlusLambda);
        es.setNth(0,ss->getMean());
        
        // Run the optimization 
        while(pbChunk.getEvaluationCost() < 100)
          es.performIteration();
      
        // State the result
        cout << "      Chunk " << vStarts[iStep] << "\t Best Value : " 
          << es.getBestEverValue() << endl;

        // Update the curve with the curve from the solution
        curve = pbChunk.UpdateRibbon(es.getBestEverX());
        }
      }
    
    // Create an overall problem
    RibbonProblem pbFinal(curve,imgWhite,0,curve.size());

    // Compute the objective values
    pbFinal.ComputeObjectives(xArea,xGreyRatio);
    cout << "   Final curve has area " << xArea 
      << " and non-white ratio " << xGreyRatio << endl;
    }

  // At this point we have hopefully optimized the ribbons, and we proceed to
  // embed them into the grey matter image to form segmentation boundaries
  
  // Compute the number of triangles in the ribbons
  unsigned int nTriangles = 0, iTriangle = 0;
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    nTriangles += vCurves[iCurve].size() * 2 - 2;
    }

  // Allocate the ribbon array
  double **xTriangles = new double*[nTriangles * 3];
  
  // Allocate and form the triangles
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];
    for(iPoint=0;iPoint<curve.size() - 1;iPoint++)
      {
      CurvePoint &p1 = curve[iPoint];
      CurvePoint &p2 = curve[iPoint + 1];
      xTriangles[iTriangle++] = (double *)(p1.xGrey.data_block());
      xTriangles[iTriangle++] = (double *)(p1.xWhite.data_block());
      xTriangles[iTriangle++] = (double *)(p2.xGrey.data_block());
      xTriangles[iTriangle++] = (double *)(p1.xWhite.data_block());
      xTriangles[iTriangle++] = (double *)(p2.xWhite.data_block());
      xTriangles[iTriangle++] = (double *)(p2.xGrey.data_block());
      }
    }

  // Allocate an array for the ribbon mask
  typedef Image<unsigned char,3> RibbonImage;
  RibbonImage::Pointer imgRibbon = RibbonImage::New();
  imgRibbon->SetRegions(imgGrey->GetBufferedRegion());
  imgRibbon->Allocate();
  imgRibbon->FillBuffer(0);
  
  // Scan convert the triangles in the grey scale image
  int dim[3];
  dim[0] = (int) imgRibbon->GetBufferedRegion().GetSize(0);
  dim[1] = (int) imgRibbon->GetBufferedRegion().GetSize(1);
  dim[2] = (int) imgRibbon->GetBufferedRegion().GetSize(2);
  drawBinaryTrianglesSheetFilled(
    imgRibbon->GetBufferPointer(),dim,xTriangles,nTriangles);

  // Save the ribbon image
  ImageFileWriter<RibbonImage>::Pointer wrib = ImageFileWriter<RibbonImage>::New();
  wrib->SetInput(imgRibbon);
  wrib->SetFileName("ribbon.hdr");
  wrib->Update();

  // Merge the ribbon image with the grey image
  typedef ImageRegionConstIterator<RibbonImage> RibbonIterator;
  typedef ImageRegionIterator<ShortImageType> GreyIterator;
  typedef ImageRegionIterator<CharImageType> WhiteIterator;
  RibbonIterator itRibbon(imgRibbon,imgRibbon->GetBufferedRegion());
  GreyIterator itGrey(imgGrey,imgGrey->GetBufferedRegion());
  WhiteIterator itWhite(imgWhite,imgWhite->GetBufferedRegion());
  while(!itGrey.IsAtEnd())
    {
    
    if(itRibbon.Get() > 0) 
      {
      if(itGrey.Get() != 0)
        itGrey.Set(0xff);
      
      if(itWhite.Get() != 0)
        itWhite.Set(0x40);
      }

    ++itGrey;++itRibbon;++itWhite;
    }

  // Save the gree scale image
  ShortWriterType::Pointer writer = ShortWriterType::New();
  writer->SetFileName("outGrey.hdr");
  writer->SetInput(imgGrey);
  writer->Update();

  // Save the white scale image
  CharWriterType::Pointer writerWhite = CharWriterType::New();
  writerWhite->SetFileName("outWhite.hdr");
  writerWhite->SetInput(imgWhite);
  writerWhite->Update();
}
