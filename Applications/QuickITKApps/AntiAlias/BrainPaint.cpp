#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"

#include "Danielsson.h"
#include "DrawTriangles.h"

#include "matrix.h"
#include "EvolutionaryStrategy.h"
#include "ConjugateGradientMethod.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

using namespace itk;
using namespace std;
using namespace pauly;

#include <vtkProperty.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkSTLReader.h>
#include <vtkTriangleFilter.h>

#include "bspline.h"

typedef vnl_vector_fixed<double,3> Vec;
typedef pauly::Vector OptVec;

struct CurvePoint {
  // A point painted manually on the grey matter
  Vec xGrey;
    
  // A matching point generated on the white matter
  Vec xWhite;
};

struct Curve {
  // A vector of curve points
  vector <CurvePoint> points;

  // A b-spline representing the grey points
  SimpleBSplineCurve *splWhite, *splGrey;
};

typedef struct Curve CurveType;
typedef Image<char,3> CharImageType;
typedef Image<short,3> ShortImageType;
typedef Image<float,3> FloatImageType;
typedef ImageFileReader<CharImageType> CharReaderType;
typedef ImageFileReader<ShortImageType> ShortReaderType;

template <class TImageType> 
void ReadImage(SmartPointer<TImageType> &target, const char *file)
{
  ImageFileReader<TImageType>::Pointer reader = 
    ImageFileReader<TImageType>::New();
  reader->SetFileName(file);
  reader->Update();
  target = reader->GetOutput();
}

template <class TImageType> 
void WriteImage(SmartPointer<TImageType> image, const char *file)
{
  ImageFileWriter<TImageType>::Pointer writer = 
    ImageFileWriter<TImageType>::New();
  writer->SetFileName(file);
  writer->SetInput(image);
  writer->Update();
}


double AreaOfTriangle(double a, double b, double c)
{
  double s = 0.5 * (a + b + c);
  double area = sqrt(s * (s-a) * (s-b) * (s-c));
  return area; 
}

void drawRibbons(const vector<CurveType> &curves, unsigned int currentCurve)
{
  // Create a renderer
  vtkRenderer *ren = vtkRenderer::New();
  ren->SetBackground(0.8,0.8,0.8);

  // Add all ribbons to it
  for(unsigned int i=0;i<curves.size();i++) 
    {
    const CurveType &curve = curves[i];
    vtkPolyData *poly = vtkPolyData::New();
    vtkCellArray *cells = vtkCellArray::New();
    vtkPoints *points = vtkPoints::New();

    // Create the array of points
    for(unsigned int j=0;j<curve.splGrey->GetInterpolationGridSize();j++)
      {
      float xGrey[3], xWhite[3];
      curve.splGrey->InterpolateGrid(j,xGrey);
      curve.splWhite->InterpolateGrid(j,xWhite);
     
              
      points->InsertNextPoint(xGrey);
      points->InsertNextPoint(xWhite);
      }

    // Create an array of cells
    for(unsigned int k=0;k<curve.points.size()-1;k++)
      {
      vtkIdType id1[3] = {k,k+1,k+3};
      vtkIdType id2[3] = {k+1,k+3,k+2};
      cells->InsertNextCell(3,id1);
      cells->InsertNextCell(3,id2);
      }

    // Put the cells and points into PD
    poly->SetPoints(points);
    poly->SetPolys(cells);

    // Delete junk
    //points->Delete();
    //cells->Delete();

    // Create a property
    vtkProperty *prop = vtkProperty::New();
    if(i == currentCurve)
      prop->SetColor(1.0f,0.0f,0.0f);
    else
      prop->SetColor(0.0f,0.0f,1.0f);
    // prop->SetRepresentationToWireframe();
    prop->SetInterpolationToFlat();

    // Create a mapper and an actor
    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInput(poly);
    
    vtkLODActor *actor = vtkLODActor::New();
    actor->SetMapper(mapper);   
    actor->SetProperty(prop);
    ren->AddActor(actor);
    }

  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  renWin->SetSize(500,500);

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  iren->Initialize();

  renWin->Render();
  iren->Start();

  iren->Delete();
  renWin->Delete();
  ren->Delete();
}



// A problem associated with optimizing grey matter partition ribbons
class RibbonProblem : public Function
{
public:
  // Construct a problem
  RibbonProblem(const CurveType &curve, FloatImageType *imgWhite,
    unsigned int iStart, unsigned int iLength, double xGreyscaleCrossingFactor = 100)
    {
    m_Curve = curve;
    m_Image = imgWhite;
    m_Factor = xGreyscaleCrossingFactor;

    // Compute the fixed components of the objective function
    if(iStart > 0 || iLength < m_Curve.points.size()) 
      {
      double xFullArea, xFullWhiteLength, xFullWhiteDistance;
      double xSubArea, xSubWhiteLength, xSubWhiteDistance;

      // Compute the full length
      m_Start = 0; m_Length = m_Curve.points.size();
      ComputeObjectives(xFullArea,xFullWhiteDistance,xFullWhiteLength);

      // Compute the sub-range 
      m_Start = iStart; m_Length = iLength;
      ComputeObjectives(xSubArea,xSubWhiteDistance,xSubWhiteLength);

      // Store the differences
      m_FixedAreaObjective = xFullArea - xSubArea;
      m_FixedWhiteLengthObjective = xFullWhiteLength - xSubWhiteLength;
      m_FixedWhiteDistanceObjective = xFullWhiteDistance - xSubWhiteDistance;
      }
    else
      {
      m_FixedAreaObjective = 0.0;
      m_FixedWhiteLengthObjective = 0.0;
      m_FixedWhiteDistanceObjective = 0.0;
      m_Start = iStart; m_Length = iLength;
      }

    // Compute the upper and lower bounds
    ComputeVectorBounds(m_Upper,m_Lower);
    }

  // Generate a vector corresponding to the current ribbon state
  OptVec GetVectorRepresentation() 
    {
    OptVec v(m_Length * 3);
    unsigned int j = 0;
    for(unsigned int i = 0;i < m_Length;i++)
      {
      v(j++) = m_Curve.points[i + m_Start].xWhite[0];
      v(j++) = m_Curve.points[i + m_Start].xWhite[1];
      v(j++) = m_Curve.points[i + m_Start].xWhite[2];
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
      m_Curve.points[i + m_Start].xWhite[0] = v(j++);
      m_Curve.points[i + m_Start].xWhite[1] = v(j++);
      m_Curve.points[i + m_Start].xWhite[2] = v(j++);
      }
    return m_Curve;
    }

  // Computes three separate objectives based for the range m_Start to m_End
  void ComputeObjectives(double &area, double &distance, double &length)
    {
    // Get a raw pointer to the data as well as the buffer access steps
    float *rawimg = m_Image->GetBufferPointer();
    unsigned int nLine = m_Image->GetBufferedRegion().GetSize(0); 
    unsigned int nSlice = m_Image->GetBufferedRegion().GetSize(1) * nLine;
    
    // Measure the costs associated with the ribbon
    double xArea = 0.0;
    double xWidth = 0.0;
    double xLenWhite = 0.0;
    double xGreyCount = 0.0;

    unsigned int iStart = (m_Start == 0) ? 0 : m_Start - 1;
    unsigned int iEnd = (m_Start + m_Length == m_Curve.points.size()) 
      ? m_Curve.points.size() : m_Start + m_Length + 1;
    
    // for(unsigned int i=iStart;i<iEnd-1;i++)
    for(unsigned int i=iStart; i<iEnd ; i++)
      {
      // Get the two points
      CurvePoint &p1 = m_Curve.points[i];
      CurvePoint &p2 = m_Curve.points[i+1];

      // Compute the edges of the ribbon
      double w1 = (p1.xGrey - p1.xWhite).two_norm();
      double w2 = (p2.xGrey - p2.xWhite).two_norm();
      double lg = (p1.xGrey - p2.xGrey).two_norm();
      double lw = (p1.xWhite - p2.xWhite).two_norm();
        
      // Width, Length integration 
      xWidth += 0.5 * (w2 + w1); 
      xLenWhite += lw;

      // Area integration
      double h = (p2.xWhite-p1.xGrey).two_norm(); 
      xArea += AreaOfTriangle(w1,lw,h);
      xArea += AreaOfTriangle(w2,lg,h);

      // Generate a list of voxels that the line passes through
      for(double t=0;t<1;t+=0.05)
        {
        // Compute the current position in the image
        double x = p1.xWhite[0] * (1-t) + p2.xWhite[0] * t;
        double y = p1.xWhite[1] * (1-t) + p2.xWhite[1] * t;
        double z = p1.xWhite[2] * (1-t) + p2.xWhite[2] * t;

        // Compute the interpolation offsets
        int ix = (int) x, iy = (int) y, iz = (int) z;
        double ux = x - ix, vx = 1.0 - ux;
        double uy = y - iy, vy = 1.0 - uy;
        double uz = z - iz, vz = 1.0 - uz;

        // Interpolate the distance image at this position
        int offset = iz * nSlice + iy * nLine + ix;
        double x000 = rawimg[offset];
        double x001 = rawimg[offset+1];
        double x010 = rawimg[offset+nLine];
        double x011 = rawimg[offset+nLine+1];
        double x100 = rawimg[offset+nSlice];
        double x101 = rawimg[offset+nSlice+1];
        double x110 = rawimg[offset+nSlice+nLine];
        double x111 = rawimg[offset+nSlice+nLine+1];

        // Compute the interpolated value - brute force way
        double val = 
          x000 * vz * vy * vx + 
          x001 * vz * vy * ux + 
          x010 * vz * uy * vx + 
          x011 * vz * uy * ux + 
          x100 * uz * vy * vx + 
          x101 * uz * vy * ux + 
          x110 * uz * uy * vx + 
          x111 * uz * uy * ux;
               
        // Find the intensity at the voxel
        xGreyCount += 0.05 * lw * val;
        }
      }

    area = xWidth * xLenWhite; //xArea;
    distance = xGreyCount;
    length = xLenWhite;
    }
  
  double evaluate(const OptVec &x)
    {
    // Check the parameter range
    double penalty = 0.0;
    for(unsigned int i=0;i<x.size();i++)
      {
      if(x(i) <= m_Lower(i))
        penalty += 1.0e10 * (1 + (m_Lower(i) - x(i)));
      else if(x(i) >= m_Upper(i))
        penalty += 1.0e10 * (1 - (m_Upper(i) - x(i)));
      }
    if(penalty > 0.0)
      return penalty;

    // Update the ribbon with the parameter vector
    UpdateRibbon(x);

    // Compute the components of the objective
    double xArea, xWhiteDistance, xWhiteLength;
    ComputeObjectives(xArea,xWhiteDistance,xWhiteLength);

    // Convert the components to a full objective function
    xArea += m_FixedAreaObjective;
    xWhiteDistance += m_FixedWhiteDistanceObjective;
    xWhiteLength += m_FixedWhiteLengthObjective;

    // Compute overall objective value based on the two components
    double obj = m_Factor * (xWhiteDistance / xWhiteLength) + xArea;

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
    OptVec sd(mean.size());
    sd.setAll(0.02);
   
    // Return a Gaussian solution space
    *space = new GaussianSS(mean,sd);
    }

private:
  // The curve undergoing optimization
  CurveType m_Curve;

  // The range of optimization
  unsigned int m_Start, m_Length;

  // Weight assigned to the percentage of grey matter crossed
  double m_Factor;

  // The components of the objective function that are not dependent on 
  // the range of the problem, and don't change between iterations
  double m_FixedAreaObjective, m_FixedWhiteLengthObjective, m_FixedWhiteDistanceObjective;

  // Upper and lower bounds for disallowing bad input
  OptVec m_Lower, m_Upper;

  // The white matter image
  FloatImageType::Pointer m_Image;
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
  vector<CurveType> vCurves;

  // Read in the curves on the surface
  FILE *fCurves = fopen(argv[3],"rt");
  char buffer[1024];
  while(!feof(fCurves)) 
    {
    // Read a line from the curve file
    fgets(buffer,1024,fCurves);
    string line(buffer);

    // A vector that holds the last read point
    Vec xLast;

    // Check the first 6 characters
    if(line.substr(0,6) == "NODES:")
      {
      CurveType c;
      vCurves.push_back(c);
      xLast.fill(-1);
      }
    else if (line.substr(0,6) == "      ")
      {
      CurvePoint p;
      istringstream iss(buffer);
      iss >> p.xGrey[0];
      iss >> p.xGrey[1];
      iss >> p.xGrey[2];

      if(p.xGrey != xLast)
        vCurves.back().points.push_back(p);
      }
    }

  unsigned int iCurve, iPoint;
  /*
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];
    CurveType curveTrim;
    for(iPoint=0;iPoint<curve.size()-1;iPoint+=3)
      curveTrim.push_back(curve[iPoint]);
    curveTrim.push_back(curve.back());
    curve = curveTrim;
    }
  */

  // At this point we've loaded the input into memory. We now compute the 
  // vector distance transform of the white matter

  // Create a copy of the white image
  char *dataWhite = new char[imgWhite->GetBufferedRegion().GetNumberOfPixels()];
  memcpy(dataWhite,imgWhite->GetBufferPointer(),imgWhite->GetBufferedRegion().GetNumberOfPixels());

  // Compute DT
  cout << "Computing the distance transform..." << endl;
  short *xDistance[3];
  edt3ddan(dataWhite,
    (int) imgWhite->GetBufferedRegion().GetSize(0),
    (int) imgWhite->GetBufferedRegion().GetSize(1),
    (int) imgWhite->GetBufferedRegion().GetSize(2),
    6,xDistance,xDistance+1,xDistance+2);

  delete dataWhite;

  // Create a signed squared distance image
  FloatImageType::Pointer imgDistance = FloatImageType::New();
  imgDistance->SetRegions(imgWhite->GetBufferedRegion());
  imgDistance->Allocate();
  float *xDistanceBuffer = imgDistance->GetBufferPointer();
  for(unsigned int iPixel=0; iPixel < 
      imgDistance->GetBufferedRegion().GetNumberOfPixels(); iPixel++)
    {
    xDistanceBuffer[iPixel] = 
      xDistance[0][iPixel] * xDistance[0][iPixel] + 
      xDistance[1][iPixel] * xDistance[1][iPixel] + 
      xDistance[2][iPixel] * xDistance[2][iPixel];
    }

  // Blur the distance image
  cout << "Blurring the distance transform..." << endl;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType,FloatImageType> BlurFilter;
  BlurFilter::Pointer fltBlur = BlurFilter::New();
  double var[] = {1.0,1.0,1.0};
  fltBlur->SetInput(imgDistance);
  fltBlur->SetVariance(var);
  fltBlur->Update();
  imgDistance = fltBlur->GetOutput();

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
  
  // Save the distance transform
  WriteImage(imgDistance,"white_distance.hdr");

  // Find the nearest white matter pixel for each of the curves
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];
    for(iPoint=0;iPoint<curve.points.size();iPoint++)
      {
      CurvePoint &p = curve.points[iPoint];
      Index<3> idx = {{
        (unsigned int)p.xGrey[0],
        (unsigned int)p.xGrey[1],
        (unsigned int)p.xGrey[2]}};

      for(unsigned int j=0;j<3;j++)
        p.xWhite[j] = p.xGrey[j] + imgOffset[j]->GetPixel(idx);
      }
    }

  // Fit b-splines to the grey and white points
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];

    // Construct a matrix Q
    BSpline1D::MatrixType Q1(curve.points.size(),3);
    BSpline1D::MatrixType Q2(curve.points.size(),3);

    for(iPoint=0;iPoint<curve.points.size();iPoint++)
      {
      Q1(iPoint,0) = curve.points[iPoint].xGrey[0];
      Q1(iPoint,1) = curve.points[iPoint].xGrey[1];
      Q1(iPoint,2) = curve.points[iPoint].xGrey[2];
      Q2(iPoint,0) = curve.points[iPoint].xWhite[0];
      Q2(iPoint,1) = curve.points[iPoint].xWhite[1];
      Q2(iPoint,2) = curve.points[iPoint].xWhite[2];
      }

    // Create the b-splines
    curve.splGrey = new SimpleBSplineCurve(10,3);
    curve.splGrey->FitToPoints(Q1);
    curve.splGrey->SetUniformInterpolationGrid(200);
    
    curve.splWhite = new SimpleBSplineCurve(10,3);
    curve.splWhite->FitToPoints(Q2);
    curve.splWhite->SetUniformInterpolationGrid(200);
    }
  
  // This is the most interesting part of the program. We attempt to 
  // improve the ribbons by making their free side pass entirely through
  // the white matter
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    // Create a curve
    CurveType &curve = vCurves[iCurve];

    // Create an overall problem
    RibbonProblem pbTotal(curve,imgDistance,0,curve.points.size());

    // Compute the objective values
    double xArea, xWhiteLength, xWhiteDistance;
    pbTotal.ComputeObjectives(xArea,xWhiteDistance,xWhiteLength);
    cout << "Ribbon " << iCurve << " of length " << curve.points.size() << " has area " << xArea 
      << " and WM distance " << (xWhiteDistance / xWhiteLength) << endl;


    // Create a list of starting points. We are not allowed to optimize the 
    // first and last points
    vector<unsigned int> vStarts;
    unsigned int nSpan = 2;
    for(unsigned int iStart=1;iStart < curve.points.size()-(nSpan+2);iStart+=1)
      vStarts.push_back(iStart);
    
    // Iteratevely optimize over ribbon segments
    for(unsigned int iIteration = 0;iIteration < 200;iIteration++)
      {
      // State the iteration
      cout << "   Running iteration " << iIteration << flush;
      drawRibbons(vCurves,iCurve);
      
      // Perform a random permutation of the starting points
      random_shuffle(vStarts.begin(),vStarts.end());

      // Peform optimizations
      for(unsigned int iStep=0;iStep<vStarts.size();iStep++)
        {
        // Create a Ribbon problem
        RibbonProblem pbChunk(curve,imgDistance,vStarts[iStep],nSpan,100);

        // Create a corresponding solution space
        SolutionSpace *ss;
        pbChunk.InitializeSolutionSpace(&ss);
        
        // Initialize the problem
        NumericalMethod *method = NULL;
        NumericalFunction nf(pbChunk);

        if(iIteration < 000)
          {
          // Use conjugate descent
          ConjugateGradientMethod *cgm = new ConjugateGradientMethod(nf,ss->getMean());
          method = cgm;
          }
        else
          {
          // Create an evolutionary strategy
          EvolutionaryStrategy *es = 
            new EvolutionaryStrategy(pbChunk,*ss,2,4,SELECTION_MuPlusLambda);
          es->setNth(0,ss->getMean());

          // Set sigma evolution parameters
          OptVec ds(ss->getMean().size());
          ds.setAll(0.6);
          es->setDeltaSigma(ds);

          // Compute upper and lower bounds
          OptVec upper,lower;
          pbChunk.ComputeVectorBounds(upper,lower);
          es->setXBounds(lower,upper);
          method = es;
          }
        
        // Run the optimization 
        while(pbChunk.getEvaluationCost() < 200 && !method->isFinished())
        //while(!method->isFinished())
          method->performIteration();
      
        // State the result
        //cout << "      Chunk " << vStarts[iStep] << "\t Best Value : " 
        //  << es.getBestEverValue() << endl;

        // Update the curve with the curve from the solution
        curve = pbChunk.UpdateRibbon(method->getBestEverX());

        delete method;
        }

      // Create an overall problem
      RibbonProblem pbFinal(curve,imgDistance,0,curve.points.size(),100);

      // Compute the objective values
      pbFinal.ComputeObjectives(xArea,xWhiteDistance,xWhiteLength);
      cout << "\t area: " << xArea 
        << "\t WM distance: " << (xWhiteDistance / xWhiteLength) 
        << "\t obj : " << pbFinal.evaluate(pbFinal.GetVectorRepresentation()) << endl;
      }

    }
    
  // At this point we have hopefully optimized the ribbons, and we proceed to
  // embed them into the grey matter image to form segmentation boundaries
  
  // Compute the number of triangles in the ribbons
  unsigned int nTriangles = 0, iTriangle = 0;
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    nTriangles += vCurves[iCurve].points.size() * 2 - 2;
    }

  // Allocate the ribbon array
  double **xTriangles = new double*[nTriangles * 3];
  
  // Allocate and form the triangles
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];
    for(iPoint=0;iPoint<curve.points.size() - 1;iPoint++)
      {
      CurvePoint &p1 = curve.points[iPoint];
      CurvePoint &p2 = curve.points[iPoint + 1];
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
  WriteImage(imgRibbon,"ribbon.hdr");

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
  WriteImage(imgGrey,"grey_with_ribbon.hdr");

  // Save the white scale image
  WriteImage(imgWhite,"white_with_ribbon.hdr");
}
