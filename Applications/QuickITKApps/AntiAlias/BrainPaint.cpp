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

// Include my BSpline library
#include "BSplineCurve.h"

// Typedef of the spline
typedef BSplineCurve<3> SplineType;
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
  SplineType splWhite, splGrey;

  // Interpolation grids for the splines (speedier)
  SplineType::EvaluationGrid gridWhite, gridGrey;
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
  typename ImageFileReader<TImageType>::Pointer reader = 
    ImageFileReader<TImageType>::New();
  reader->SetFileName(file);
  reader->Update();
  target = reader->GetOutput();
}

template <class TImageType> 
void WriteImage(SmartPointer<TImageType> image, const char *file)
{
  typename ImageFileWriter<TImageType>::Pointer writer = 
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
    for(unsigned int j=0;j<curve.gridGrey.size();j++)
      {
      SplineType::Point xGrey = curve.splGrey.EvaluateGridPoint(curve.gridGrey, j, 0);
      SplineType::Point xWhite = curve.splWhite.EvaluateGridPoint(curve.gridWhite, j, 0);
      
      points->InsertNextPoint(xGrey.data_block());
      points->InsertNextPoint(xWhite.data_block());
      }

    // Create an array of cells
    for(unsigned int k=0;k<curve.gridGrey.size()-1;k++)
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
    m_Start = iStart; m_Length = iLength;

    // Compute the bounds on the image
    m_ImageMax[0] = (double) m_Image->GetBufferedRegion().GetSize(0) - 1;
    m_ImageMax[1] = (double) m_Image->GetBufferedRegion().GetSize(1) - 1;
    m_ImageMax[2] = (double) m_Image->GetBufferedRegion().GetSize(2) - 1;

    // Compute the upper and lower bounds
    ComputeVectorBounds(m_Upper,m_Lower);
    }

  // Generate a vector corresponding to the current ribbon state
  OptVec GetVectorRepresentation() 
    {
    OptVec v(m_Length * 3);
    
    unsigned int j = 0;
    for(unsigned int i = 0;i < m_Length; i++)
      {
      SplineType::Point point = m_Curve.splWhite.GetControlPoint(i + m_Start);
      v(j++) = point(0); v(j++) = point(1); v(j++) = point(2);
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
      SplineType::Point point;
      point(0) = v(j++); point(1) = v(j++); point(2) = v(j++);
      m_Curve.splWhite.SetControlPoint(i + m_Start, point);
      }
    return m_Curve;
    }

  // Different objectives computed in the optimization
  struct Objectives {
    double xArea;
    double xWidth;
    double xWhiteDistance;
    double xWhiteLength;
    double xWhiteIrregularity;
  };

  // Computes three separate objectives based for the range m_Start to m_End
  Objectives ComputeObjectives()
    {
    // Get a raw pointer to the data as well as the buffer access steps
    float *rawimg = m_Image->GetBufferPointer();
    unsigned int nLine = m_Image->GetBufferedRegion().GetSize(0); 
    unsigned int nSlice = m_Image->GetBufferedRegion().GetSize(1) * nLine;
    
    // Measure the costs associated with the ribbon
    Objectives obj;
    obj.xArea = 0.0;
    obj.xWidth = 0.0;
    obj.xWhiteLength = 0.0;
    obj.xWhiteDistance = 0.0;
    obj.xWhiteIrregularity = 0.0;

    // Integrate the function over the ribbon
    unsigned int nPoints = m_Curve.gridWhite.size();
    vector< SplineType::Point > xWhite(nPoints); 
    vector< SplineType::Point > xGrey(nPoints);
    vector<double> xImageVal(nPoints); 

    for(unsigned int j=0;j<nPoints;j++)
      {
      xWhite[j] = m_Curve.splWhite.EvaluateGridPoint(m_Curve.gridWhite,j,0);
      xGrey[j] = m_Curve.splGrey.EvaluateGridPoint(m_Curve.gridGrey,j,0);
      
      // Compute the interpolation offsets
      double x = xWhite[j](0), y = xWhite[j](1), z = xWhite[j](2);
      int ix = (int) x, iy = (int) y, iz = (int) z;
      double ux = x - ix, vx = 1.0 - ux;
      double uy = y - iy, vy = 1.0 - uy;
      double uz = z - iz, vz = 1.0 - uz;

      // Check the bounds on the coordinates
      double xPenalty = 0.0, xPenaltyBase = 1e6;
      if(x < 0.0) xPenalty += - ux;
      if(y < 0.0) xPenalty += - uy;
      if(z < 0.0) xPenalty += - uz;
      if(x >= m_ImageMax[0]) xPenalty += (ux - m_ImageMax[0]);
      if(y >= m_ImageMax[1]) xPenalty += (uy - m_ImageMax[1]);
      if(z >= m_ImageMax[2]) xPenalty += (uz - m_ImageMax[2]);
      
      // Penalty precludes image computatino
      if(xPenalty > 0.0)
        {
        xImageVal[j] = xPenaltyBase * (1 + xPenalty);
        }        
      else
        {
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
        xImageVal[j] = 
          x000 * vz * vy * vx + 
          x001 * vz * vy * ux + 
          x010 * vz * uy * vx + 
          x011 * vz * uy * ux + 
          x100 * uz * vy * vx + 
          x101 * uz * vy * ux + 
          x110 * uz * uy * vx + 
          x111 * uz * uy * ux;
        }
      }

    // Keep track of the sum of distances and the sum of squares for variance computation
    double xWhiteLenSqr = 0.0;
    
    // Compute the shape and white matter membership over the ribbon
    for(unsigned int i=0; i < xWhite.size() - 1 ; i++)
      {
      
      // Compute the edges of the ribbon square
      double w1 = (xGrey[i] - xWhite[i]).two_norm();
      double w2 = (xGrey[i+1] - xWhite[i+1]).two_norm();
      double lg = (xGrey[i] - xGrey[i+1]).two_norm();
      double lw = (xWhite[i] - xWhite[i+1]).two_norm();
        
      // Width, Length integration 
      obj.xWidth += 0.5 * (w2 + w1); 
      obj.xWhiteLength += lw;
      xWhiteLenSqr += lw * lw;
      obj.xWhiteDistance += 0.5 * (xImageVal[i] + xImageVal[i+1]) * lw;

      // Area integration
      double h = (xWhite[i+1]-xGrey[i]).two_norm(); 
      obj.xArea += AreaOfTriangle(w1,lw,h);
      obj.xArea += AreaOfTriangle(w2,lg,h);
      }

    // Compute the variance in the white length
    obj.xWhiteIrregularity = (xWhiteLenSqr - (obj.xWhiteLength * obj.xWhiteLength) / nPoints) / (nPoints - 1);
    
    return obj;
    }
  
  double evaluate(const OptVec &x)
    {
    // Update the ribbon with the parameter vector
    UpdateRibbon(x);

    // Compute the components of the objective
    Objectives obj = ComputeObjectives();

    // Compute overall objective value based on the two components
    double result = m_Factor * (obj.xWhiteDistance / obj.xWhiteLength) 
       + obj.xArea + 10000 * obj.xWhiteIrregularity;

    // Increase the iterations counter
    setEvaluationCost(getEvaluationCost() + 1);

    // Return the objective
    return result;
    }

  // Generate the initial search space given the current state
  void InitializeSolutionSpace(SolutionSpace **space) 
    {
    // Create a vector for the current state
    OptVec mean = GetVectorRepresentation();

    // Generate a standard deviation vector
    OptVec sd(mean.size());
    sd.setAll(0.01);
   
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

  // Image bounds
  double m_ImageMax[3];
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
  /*
  cout << "Blurring the distance transform..." << endl;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType,FloatImageType> BlurFilter;
  BlurFilter::Pointer fltBlur = BlurFilter::New();
  double var[] = {1.0,1.0,1.0};
  fltBlur->SetInput(imgDistance);
  fltBlur->SetVariance(var);
  fltBlur->Update();
  imgDistance = fltBlur->GetOutput();
  */

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
    unsigned int nPoints = curve.points.size();
    
    // Construct a matrix Q
    SplineType::Point *Q1 = new SplineType::Point[nPoints];
    SplineType::Point *Q2 = new SplineType::Point[nPoints];

    // Fill the matrix    
    for(iPoint=0;iPoint<nPoints;iPoint++)
      {
      Q1[iPoint] = curve.points[iPoint].xGrey;
      Q2[iPoint] = curve.points[iPoint].xWhite;      
      }

    // Create the b-splines
    curve.splGrey.SetNumberOfControlPoints(16);
    curve.splGrey.FitToPoints(nPoints,Q1);
    curve.splGrey.CreateUniformEvaluationGrid(200, curve.gridGrey);
    
    curve.splWhite.SetNumberOfControlPoints(16);
    curve.splWhite.FitToPoints(nPoints,Q2);
    curve.splWhite.CreateUniformEvaluationGrid(200, curve.gridWhite);
    }
  
  // This is the most interesting part of the program. We attempt to 
  // improve the ribbons by making their free side pass entirely through
  // the white matter
  for(iCurve=0; iCurve < vCurves.size() ;iCurve++)
    {
    // Create a curve
    CurveType &curve = vCurves[iCurve];

    // Create an overall problem
    double factor = 20000;
    unsigned int nControls = curve.splWhite.GetNumberOfControlPoints();
    RibbonProblem pbTotal(curve,imgDistance,0,nControls,factor);

    // Compute the objective values
    RibbonProblem::Objectives obj = pbTotal.ComputeObjectives();
    cout << "Ribbon " << iCurve << " of length " << curve.points.size() << " has area " << obj.xArea 
       << " irreg " << obj.xWhiteIrregularity  
      << " and WM distance " << (obj.xWhiteDistance / obj.xWhiteLength) << endl;

    // Create a list of starting points. We are not allowed to optimize the 
    // first and last points
    vector<unsigned int> vStarts;
    unsigned int nSpan = 4;
    for(unsigned int iStart=1;iStart < curve.splWhite.GetNumberOfControlPoints() - nSpan;iStart++)
      {
      vStarts.push_back(iStart);
      }
    
    // Iteratevely optimize over ribbon segments
    for(unsigned int iIteration = 0;iIteration < 200;iIteration++)
      {
      // drawRibbons(vCurves,iCurve);
       
      // State the iteration
      cout << "   Running iteration " << iIteration << flush;
            
      // Perform a random permutation of the starting points
      random_shuffle(vStarts.begin(),vStarts.end());

      // Peform optimizations
      for(unsigned int iStep=0;iStep<vStarts.size();iStep++)
        {
        // Create a Ribbon problem
        RibbonProblem pbChunk(curve,imgDistance,vStarts[iStep],nSpan,factor);

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
      RibbonProblem pbFinal(curve,imgDistance,0,nControls,factor);

      // Compute the objective values
      RibbonProblem::Objectives obj = pbFinal.ComputeObjectives();
      cout << "\t area: " << obj.xArea 
        << "\t WM distance: " << (obj.xWhiteDistance / obj.xWhiteLength) 
        << "\t irreg: " << obj.xWhiteIrregularity 
        << "\t obj : " << pbFinal.evaluate(pbFinal.GetVectorRepresentation()) << endl;
      }

      // drawRibbons(vCurves,iCurve);

    }
    
  // At this point we have hopefully optimized the ribbons, and we proceed to
  // embed them into the grey matter image to form segmentation boundaries
  
  // Compute the number of triangles in the ribbons
  unsigned int nTriangles = 0, iTriangle = 0;
  unsigned int nTotalPoints = 0, iTotalPoints = 0;
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    unsigned int nPoints = vCurves[iCurve].gridWhite.size();
    nTriangles += nPoints * 2 - 2;
    nTotalPoints += nPoints;
    }

  // Allocate the ribbon array
  SplineType::Point *xPoints = new SplineType::Point[nTotalPoints * 2];
  double **xTriangles = new double*[nTriangles * 3];
  
  // Allocate and form the triangles
  for(iCurve=0;iCurve<vCurves.size();iCurve++)
    {
    CurveType &curve = vCurves[iCurve];
    unsigned int nPoints = curve.gridWhite.size();

    xPoints[iTotalPoints++] = curve.splGrey.EvaluateGridPoint(curve.gridGrey,0,0);
    xPoints[iTotalPoints++] = curve.splWhite.EvaluateGridPoint(curve.gridWhite,0,0);
    for(iPoint=1;iPoint<nPoints;iPoint++)
      {
      SplineType::Point &g1 = xPoints[iTotalPoints-2];
      SplineType::Point &w1 = xPoints[iTotalPoints-1];
      SplineType::Point &g2 = xPoints[iTotalPoints++] 
        = curve.splGrey.EvaluateGridPoint(curve.gridGrey,iPoint,0);
      SplineType::Point &w2 =xPoints[iTotalPoints++] 
        = curve.splWhite.EvaluateGridPoint(curve.gridWhite,iPoint,0);

      xTriangles[iTriangle++] = g1.data_block();
      xTriangles[iTriangle++] = w1.data_block();
      xTriangles[iTriangle++] = g2.data_block();
      xTriangles[iTriangle++] = w1.data_block();
      xTriangles[iTriangle++] = w2.data_block();
      xTriangles[iTriangle++] = g2.data_block();
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
