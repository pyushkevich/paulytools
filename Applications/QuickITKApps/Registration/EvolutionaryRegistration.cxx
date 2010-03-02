#include "itkOrientedRASImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInterpolateImageFunction.h"
#include "itkImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include <iostream>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkNormalVariateGenerator.h"

using namespace std;
using namespace itk;


//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() { m_LastMetricValue = 0.0; };
public:
  typedef itk::OnePlusOneEvolutionaryOptimizer     OptimizerType;
  typedef   const OptimizerType *                  OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer = 
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      double currentValue = optimizer->GetValue();
      // Only print out when the Metric value changes
      if( vcl_fabs( m_LastMetricValue - currentValue ) > 1e-7 )
        { 
        std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << currentValue << "   ";
        std::cout << optimizer->GetCurrentPosition() << std::endl;
        m_LastMetricValue = currentValue;
        }
    }
private:
  double m_LastMetricValue;
};

void PrintMatrix(vnl_matrix<double> mat, const char *fmt, const char *prefix)
{
  // Print each row and column of the matrix
  char buffer[256];
  for(size_t i = 0; i < mat.rows(); i++)
    {
    cout << prefix;
    for(size_t j = 0; j < mat.columns(); j++)
      {
      sprintf(buffer, fmt, mat(i,j));
      cout << buffer;
      }
    cout << endl;
    }
}

// Check if input matrix is orthogonal to within tolerance
template<class TScalarType, unsigned int VDim>
bool
MatrixIsOrthogonal(
 const itk::Matrix<double,VDim,VDim> &matrix,
 double tolerance )
{
  typename itk::Matrix<double,VDim,VDim>::InternalMatrixType test =
    matrix.GetVnlMatrix() * matrix.GetTranspose();

  if( !test.is_identity( tolerance ) )
    {
    return false;
    }

  return true;
}


template<unsigned int VDim>
void ReadMatrix(const char *fname, itk::Matrix<double,VDim+1,VDim+1> &mat)
  {
  ifstream fin(fname);
  for(size_t i = 0; i < VDim+1; i++)
    for(size_t j = 0; j < VDim+1; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        std::cerr << "Unable to read matrix " <<  fname << std::endl;
        }
  fin.close();
  }


template <typename TAffine, unsigned int VDim>
void 
ReadTransformFile(typename TAffine::Pointer atran, char * fn, bool isrigid)
{
  // Read transform based on type
  string tran_fn(fn);
  string format; 

  if (tran_fn.compare(tran_fn.size()-4,4,".mat") == 0)
    format="matrix";
  else
    format="itk";
    
    
  if(format=="itk")
    {
    cout << "ITK transform file " << fn << endl;
    typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> MOTBType;
    typedef itk::AffineTransform<double, VDim> AffTran;
    typedef itk::VersorRigid3DTransform<double> RigTran;
    itk::TransformFactory<MOTBType>::RegisterTransform();
    itk::TransformFactory<AffTran>::RegisterTransform();
    itk::TransformFactory<RigTran>::RegisterTransform();

    itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
    fltReader->SetFileName(fn);
    fltReader->Update();

    itk::TransformBase *base = fltReader->GetTransformList()->front();
    MOTBType *motb = dynamic_cast<MOTBType *>(base);

    if(motb)
      {
      if (isrigid)
        {

        // get rotation matrix
        vnl_matrix_fixed<double, 4, 4> mat;
        {
        for(size_t r = 0; r < 3; r++)
          {
          for(size_t c = 0; c < 3; c++)
            {
            mat(r,c) = motb->GetMatrix()(r,c);
          }
          mat(r,3) = motb->GetOffset()[r];
          }
        mat(2,0) *= -1; mat(2,1) *= -1;
        mat(0,2) *= -1; mat(1,2) *= -1;
        mat(0,3) *= -1; mat(1,3) *= -1;
        }


        itk::Matrix<double,VDim,VDim> amat;
        amat.GetVnlMatrix().update( mat.extract(VDim, VDim));
        // Force orthogonality
        //PrintMatrix( amat.GetVnlMatrix(), "%12.12f ", "    ");
        //double det = vnl_determinant( amat.GetVnlMatrix() );
        //amat = amat / det;
        PrintMatrix( amat.GetVnlMatrix(), "%12.5f ", "    ");
        //if (MatrixIsOrthogonal<double,VDim>( amat, 1e-5 ))
        //  cout << "orthogonal" << endl;
        //else
        //  cout << "nooooooooooooooo" << endl;

        atran->SetMatrix(amat);
        }
        else
          atran->SetMatrix(motb->GetMatrix());
        atran->SetOffset(motb->GetOffset());
        }
      
    }
  else if(format=="matrix")
    {
    cout << "Matrix transform file " << fn << endl;
    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> matrix;
    itk::Matrix<double,VDim,VDim> amat, amatorig;
    itk::Vector<double, VDim> aoff, aofforig;

    ReadMatrix<VDim>(fn, matrix);
    amat.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(VDim, VDim));
    aoff.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(VDim).extract(VDim));

    amatorig.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(VDim, VDim));
    aofforig.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(VDim).extract(VDim));
    

    // Extrernal matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
    vnl_vector<double> v_lps_to_ras(VDim, 1.0);
    v_lps_to_ras[0] = v_lps_to_ras[1] = -1.0;
    vnl_diag_matrix<double> m_lps_to_ras(v_lps_to_ras);
    vnl_matrix<double> mold = amat.GetVnlMatrix();
    amat.GetVnlMatrix().update(m_lps_to_ras * mold * m_lps_to_ras);
    aoff.GetVnlVector().update(m_lps_to_ras * aoff.GetVnlVector());

    // Put the values in the transform
    if (isrigid)
      {
      PrintMatrix( amatorig.GetVnlMatrix(), "%12.5f ", "    ");
      //if (MatrixIsOrthogonal<double,VDim>( amatorig, 1e-5 ))
      //  cout << "orthogonal" << endl;
      //else
      //  cout << "nooooooooooooooo" << endl;
      atran->SetMatrix(amatorig);
      }
    else
      atran->SetMatrix(amat);
    atran->SetOffset(aoff);
    }


}

template <class TPixel, unsigned int VDim>
int
EvolutionaryRegistration(int argc, char * argv[])
{
  typedef  TPixel  PixelType;

  typedef itk::OrientedRASImage<TPixel, VDim> ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef vnl_vector_fixed<double, VDim> RealVector;
  typedef vnl_vector_fixed<int, VDim> IntegerVector;

  typedef ImageType FixedImageType;
  typedef ImageType MovingImageType;
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  typedef itk::AffineTransform< double, VDim > AffineTransformType;
  // TODO Generalize to 2D rigid
  typedef itk::VersorRigid3DTransform< double > VersorRigid3DTransformType;
    typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> MOTBType;
  typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType,MovingImageType> MattesMutualInformationImageToImageMetricType;
  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric< FixedImageType,MovingImageType> NormalizedMutualInformationHistogramImageToImageMetricType;
  typedef itk::MeanSquaresImageToImageMetric< FixedImageType,MovingImageType> MeanSquaresImageToImageMetricType;
  typedef itk::NormalizedCorrelationImageToImageMetric< FixedImageType,MovingImageType> NormalizedCorrelationImageToImageMetricType;
  typedef itk::OnePlusOneEvolutionaryOptimizer           OptimizerType;
  typedef itk::LinearInterpolateImageFunction<
                                    MovingImageType,
                                    double             > InterpolatorType;
  typedef itk::ImageRegistrationMethod<
                                    FixedImageType,
                                    MovingImageType    > RegistrationType;
  // Create registration
  typename RegistrationType::Pointer   registration  = RegistrationType::New();

  // Create optimizer
  typename OptimizerType::Pointer      optimizer     = OptimizerType::New();

  // Create readers
  typename FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  typename MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  fixedImageReader->Update();

  registration->SetFixedImageRegion( 
       fixedImageReader->GetOutput()->GetBufferedRegion() );

  typename FixedImageType::SizeType fixedsize = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize(); 

  // Create the appropriate metric
  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType> MetricType;
  typename MetricType::Pointer metric;


  std::string metric_name( argv[6] );

  if(!strcmp(metric_name.c_str(),"mmi"))
    {
    //MattesMutualInformationImageToImageMetricType *metric = dynamic_cast<MattesMutualInformationImageToImageMetricType *>( basemetric);
    metric = MattesMutualInformationImageToImageMetricType::New();
    //metric->SetNumberOfHistogramBins( 40 );
    //metric->SetNumberOfSpatialSamples( 10000 );
    registration->SetMetric( metric  );
    optimizer->MaximizeOff();
    }
  else if(!strcmp(metric_name.c_str(),"nmi"))
    {
    metric = NormalizedMutualInformationHistogramImageToImageMetricType::New();
    registration->SetMetric( metric  );
    optimizer->MaximizeOn();
    }
  else if(!strcmp(metric_name.c_str(),"msq"))
    {
    metric = MeanSquaresImageToImageMetricType::New();
    registration->SetMetric( metric  );
    optimizer->MaximizeOff();
    }
  else if(!strcmp(metric_name.c_str(),"ncc"))
    {
    metric = NormalizedCorrelationImageToImageMetricType::New();
    registration->SetMetric( metric  );
    optimizer->MaximizeOff();
    }
  else 
    {
    cerr << "Unknown metric type " << metric_name << endl;
    return EXIT_FAILURE;
    };


    cout << "Set metric " << metric_name << endl;
  

  // Create the appropriate transform and initialize transform with appropriate parameters
  typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> MOTBType;
  typename MOTBType::Pointer transform, finaltransform;
  
  std::string fn( argv[4] );
  bool isinit = (fn.compare(fn.size()-4,4,"none") != 0);
  std::string transform_name( argv[5] );

  int nParameters;
  vnl_matrix_fixed<double, VDim+1, VDim+1> amat(0.0);
  vnl_vector_fixed<double, VDim+1> atmp(1.0);

  typename RegistrationType::ParametersType initialParameters;
  typename OptimizerType::ScalesType paramScales;

  if(!strcmp(transform_name.c_str(),"rig"))
    {
    transform = VersorRigid3DTransformType::New();
    transform->SetIdentity();
    if (isinit)
      {
      ReadTransformFile<MOTBType, VDim>( transform, argv[4], true);
      cout << "Transform file " << fn << " read" << endl;
      }
    else
      cout << "Initial Transform set to identity" << endl;
    registration->SetTransform(     transform     );
    nParameters = transform->GetNumberOfParameters();
    initialParameters.SetSize( nParameters );
    paramScales.SetSize( nParameters );
    amat.update(transform->GetMatrix().GetVnlMatrix(), 0, 0);
    atmp.update(transform->GetOffset().GetVnlVector(), 0);
    amat.set_column(VDim, atmp);
    initialParameters = transform->GetParameters();
    for (size_t i=0; i < VDim; i++)
      paramScales[i] = 1.0;
    // 0.01 radian of rotation is set equivalent to 0.1*voxel size
    for (size_t i=0; i < VDim ; i++)
      paramScales[VDim + i] = fixedsize[i] * 0.1/0.01;

    }
  else if(!strcmp(transform_name.c_str(),"aff"))
    {
    transform = AffineTransformType::New();
    transform->SetIdentity();
    if (isinit)
      {
      ReadTransformFile<MOTBType, VDim>( transform, argv[4], false);
      cout << "Transform file " << fn << " read" << endl;
      }
    else
      cout << "Initial Transform set to identity" << endl;
    registration->SetTransform(     transform     );
    nParameters = transform->GetNumberOfParameters();
    initialParameters.SetSize( nParameters );
    paramScales.SetSize( nParameters );
    amat.update(transform->GetMatrix().GetVnlMatrix(), 0, 0);
    atmp.update(transform->GetOffset().GetVnlVector(), 0);
    amat.set_column(VDim, atmp);
    initialParameters = transform->GetParameters();
    paramScales.fill(1.0);
    for (size_t i=0; i < VDim; i++)
      paramScales[i*VDim + i] = 1.0;
    // 0.01 radian of rotation is set equivalent to 0.1*voxel size
    for (size_t i=0; i < VDim ; i++)
      paramScales[VDim*VDim + i] = fixedsize[i] * 0.1/0.01;

    }
  else 
    {
    cerr << "Unknown transform type " << transform_name << endl;
    return EXIT_FAILURE;
    }; 


  typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();


  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  
  optimizer->SetScales( paramScales );
  registration->SetInitialTransformParameters( initialParameters );

  
  
  // Output some information about initial transform and parameters
  PrintMatrix( amat, "%12.5f ", "    ");
  std::cout << "N = " << nParameters << std::endl << " params: " << registration->GetInitialTransformParameters() << std::endl;
  std::cout << "param scales are: " << optimizer->GetScales() << std::endl;
//std::cout << "Fixed paramaters are " ;
//if (ttype.compare(ttype.size()-3,3,"rig") == 0)  cout << rigidtransform->GetFixedParameters() ; else cout << affinetransform->GetFixedParameters(); std::endl;

  typedef itk::Statistics::NormalVariateGenerator  GeneratorType;

  typename GeneratorType::Pointer generator = GeneratorType::New();


  generator->Initialize(12345);

  optimizer->SetNormalVariateGenerator( generator );
  optimizer->Initialize( atof(argv[7]), atof(argv[8]), atof(argv[9]) );
  optimizer->SetEpsilon( atof(argv[10]) );
  optimizer->SetMaximumIteration( atoi(argv[11]) );


  // Create the Command observer and register it with the optimizer.
  //
  typename CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
 
  optimizer->DebugOn();
  //cout << "Debug state of optimizer is " << optimizer->GetDebug() << endl;
  try 
    { 
    std::cout << "Registration starts!" << std::endl;
    registration->StartRegistration(); 
    std::cout << "Registration completed!" << std::endl;
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return EXIT_FAILURE;
    } 

  typename RegistrationType::ParametersType finalParameters = registration->GetLastTransformParameters();
  
  
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  
  double bestValue = optimizer->GetValue();


  // Print out results
  //
  std::cout << "Result = " << std::endl;
  std::cout << " Final paramaters = " << finalParameters  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;



  typedef itk::ResampleImageFilter< 
                            MovingImageType, 
                            FixedImageType >    ResampleFilterType;

  typename AffineTransformType::Pointer finalaffineTransform = AffineTransformType::New();
  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();

  if(!strcmp(transform_name.c_str(),"rig"))
    {
    finaltransform = VersorRigid3DTransformType::New();
    finaltransform->SetParameters( finalParameters );
    resample->SetTransform( finaltransform );
  
    }
  else
    {
    finaltransform = AffineTransformType::New();
    finaltransform->SetParameters( finalParameters );
    resample->SetTransform( finaltransform );
    
    }

  resample->SetInput( movingImageReader->GetOutput() );

  typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( 0 );


  typedef ImageType OutputImageType;
  typedef itk::CastImageFilter< 
                        FixedImageType,
                        OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  typename WriterType::Pointer      writer =  WriterType::New();
  typename CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );

  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  return EXIT_SUCCESS;
};

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << "outputImagefile ";
    std::cerr << "InitialTransformFile(none for no initial) ";
    std::cerr << "RegType(rig or aff) ";
    std::cerr << "Metric ";
    std::cerr << "InitialRadius Grow Shrink Epsilon MaxIter ";
    return EXIT_FAILURE;
    }

    return EvolutionaryRegistration<double,3>(argc, argv);
  
}
