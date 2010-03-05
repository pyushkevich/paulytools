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
#include "itkSymmetricImageToImageMetric.h"
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

void ras_write(vnl_matrix<double> mat, const char *fname)
{
  

  ofstream fout(fname);
  for(size_t i = 0; i < mat.rows(); i++)
    for(size_t j = 0; j < mat.columns(); j++)
      fout << mat[i][j] << (j < mat.columns()-1 ? " " : "\n");

  fout.close();
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


template<unsigned int VDim>
void Flip_LPS_to_RAS(itk::Matrix<double,VDim+1,VDim+1> &matrix, itk::Matrix<double,VDim,VDim> &amat, itk::Vector<double, VDim> &aoff)
{   

    // Convert lps to ras
    vnl_vector<double> v_ras_to_lps(VDim, 1.0);
    v_ras_to_lps[0] = v_ras_to_lps[1] = -1.0;
    vnl_diag_matrix<double> m_ras_to_lps(v_ras_to_lps);

    vnl_matrix<double> amatvnl = amat.GetVnlMatrix();
    amatvnl = m_ras_to_lps * amatvnl * m_ras_to_lps;
    vnl_vector_fixed<double, VDim > aoffs ;
    vnl_vector_fixed<double, VDim + 1> aoffl ;
    aoffs = m_ras_to_lps * aoff.GetVnlVector();
    aoffl.fill(1.0);
    for (size_t i=0; i<VDim; i++)
      aoffl(i) = aoffs(i);    

    matrix.GetVnlMatrix().set_identity();
    matrix.GetVnlMatrix().update( amatvnl, 0, 0);
    matrix.GetVnlMatrix().set_column(VDim, aoffl);


}


template<unsigned int VDim>
void Flip_RAS_to_LPS(itk::Matrix<double,VDim+1,VDim+1> &matrix, itk::Matrix<double,VDim,VDim> &amat, itk::Vector<double, VDim> &aoff)
{

    amat.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(VDim, VDim));
    aoff.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(VDim).extract(VDim));


    // Extrernal matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
    vnl_vector<double> v_lps_to_ras(VDim, 1.0);
    v_lps_to_ras[0] = v_lps_to_ras[1] = -1.0;
    vnl_diag_matrix<double> m_lps_to_ras(v_lps_to_ras);
    vnl_matrix<double> mold = amat.GetVnlMatrix();
    amat.GetVnlMatrix().update(m_lps_to_ras * mold * m_lps_to_ras);
    aoff.GetVnlVector().update(m_lps_to_ras * aoff.GetVnlVector());


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
/*      if (isrigid)
        {

        // get rotation matrix
        vnl_matrix_fixed<double, VDim + 1, VDim + 1> mat;
        {
        for(size_t r = 0; r < VDim; r++)
          {
          for(size_t c = 0; c < VDim; c++)
            {
            mat(r,c) = motb->GetMatrix()(r,c);
          }
          mat(r,VDim) = motb->GetOffset()[r];
          }
        // TODO Generalize with VDim
        mat(2,0) *= -1; mat(2,1) *= -1;
        mat(0,2) *= -1; mat(1,2) *= -1;
        mat(0,3) *= -1; mat(1,3) *= -1;
        }


        itk::Matrix<double,VDim,VDim> amat;
        amat.GetVnlMatrix().update( mat.extract(VDim, VDim));
        // Force orthogonality
        //PrintMatrix( amat.GetVnlMatrix(), "%12.12f ", "    ");
        //double det = vnl_determinant( amat.GetVnlMatrix() );
        //double normfactor = vcl_pow( 1.0/det, 1/(double) VDim ); 
        //amat = normfactor * amat.GetVnlMatrix();
        //PrintMatrix( amat.GetVnlMatrix(), "%12.12f ", "    ");
        //if (MatrixIsOrthogonal<double,VDim>( amat, 1e-10 ))
        //  cout << "orthogonal" << endl;
        //else
        //  cout << "nooooooooooooooo" << endl;

        atran->SetMatrix(amat);
        }
        else */
          atran->SetMatrix(motb->GetMatrix());
        atran->SetOffset(motb->GetOffset());
        }
      
    }
  else if(format=="matrix")
    {
    cout << "Matrix transform file " << fn << endl;
    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> matrix;
    itk::Matrix<double,VDim,VDim> amat ;
    itk::Vector<double, VDim> aoff ;

    ReadMatrix<VDim>(fn, matrix);
    //cout << "matrix after reading.. " << endl;
    //  PrintMatrix( matrix.GetVnlMatrix(), "%12.5f ", "    ");

    Flip_RAS_to_LPS(matrix, amat, aoff);
    // Put the values in the transform
/*
    if (isrigid)
      {
      //PrintMatrix( amatorig.GetVnlMatrix(), "%12.5f ", "    ");
      //if (MatrixIsOrthogonal<double,VDim>( amatorig, 1e-5 ))
      //  cout << "orthogonal" << endl;
      //else
      //  cout << "nooooooooooooooo" << endl;
      atran->SetMatrix(amatorig);
      }
    else
*/
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
  typedef itk::SymmetricImageToImageMetric< FixedImageType,MovingImageType> SymmetricImageToImageMetricType;
  typename SymmetricImageToImageMetricType::Pointer   symmmetric  = SymmetricImageToImageMetricType::New();

  // Create optimizer
  typename OptimizerType::Pointer      optimizer     = OptimizerType::New();

  // Create readers
  typename FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  typename MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  symmmetric->SetFixedImage(    fixedImageReader->GetOutput()    );
  symmmetric->SetMovingImage(   movingImageReader->GetOutput()   );

  fixedImageReader->Update();
  movingImageReader->Update();

  symmmetric->SetFixedImageRegion( 
       fixedImageReader->GetOutput()->GetBufferedRegion() );

  typename FixedImageType::SizeType fixedsize = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize(); 

  // Create the appropriate metric
  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType> MetricType;
  typename MetricType::Pointer metric;


  std::string reg_symm( argv[7] );
  if(!strcmp(reg_symm.c_str(),"symm"))
    symmmetric->UseSymmetricOn();
  else if (!strcmp(reg_symm.c_str(),"asym"))
    symmmetric->UseSymmetricOff();
 else
    {
    cerr << "Unknown symmetry type " << reg_symm << endl;
    return EXIT_FAILURE;
    };

    


  std::string metric_name( argv[6] );

  if(!strcmp(metric_name.c_str(),"mmi"))
    {
    //MattesMutualInformationImageToImageMetricType *metric = dynamic_cast<MattesMutualInformationImageToImageMetricType *>( basemetric);
    metric = MattesMutualInformationImageToImageMetricType::New();
    //metric->SetNumberOfHistogramBins( 40 );
    //metric->SetNumberOfSpatialSamples( 10000 );
    symmmetric->SetAsymMetric( metric  );
    optimizer->MaximizeOff();
    }
  else if(!strcmp(metric_name.c_str(),"nmi"))
    {
    metric = NormalizedMutualInformationHistogramImageToImageMetricType::New();
    symmmetric->SetAsymMetric( metric  );
    optimizer->MaximizeOn();
    }
  else if(!strcmp(metric_name.c_str(),"msq"))
    {
    metric = MeanSquaresImageToImageMetricType::New();
    symmmetric->SetAsymMetric( metric  );
    optimizer->MaximizeOff();
    }
  else if(!strcmp(metric_name.c_str(),"ncc"))
    {
    metric = NormalizedCorrelationImageToImageMetricType::New();
    symmmetric->SetAsymMetric( metric  );
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
  itk::Matrix<double, VDim+1, VDim+1> amat;

  typename RegistrationType::ParametersType initialParameters;
  typename OptimizerType::ScalesType paramScales;

  if(!strcmp(transform_name.c_str(),"rig"))
    {
    transform = VersorRigid3DTransformType::New();
    transform->SetIdentity();
    if (isinit)
      {
      ReadTransformFile<MOTBType, VDim>( transform, argv[4], true);
      //cout << "Transform file " << fn << " read" << endl;
      //std::cout << "Just after reading params: " << transform->GetParameters() << std::endl;
      }
    else
      cout << "Initial Transform set to identity" << endl;
    symmmetric->SetTransform(     transform     );
    nParameters = transform->GetNumberOfParameters();
    initialParameters.SetSize( nParameters );
    paramScales.SetSize( nParameters );
    initialParameters = transform->GetParameters();
    itk::Matrix<double, VDim, VDim> tmat = transform->GetMatrix();
    itk::Vector<double, VDim> toff = transform->GetOffset();
    Flip_LPS_to_RAS( amat, tmat, toff);
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
      //cout << "Transform file " << fn << " read" << endl;
      }
    else
      cout << "Initial Transform set to identity" << endl;
    symmmetric->SetTransform(     transform     );
    nParameters = transform->GetNumberOfParameters();
    initialParameters.SetSize( nParameters );
    paramScales.SetSize( nParameters );
    initialParameters = transform->GetParameters();
    itk::Matrix<double, VDim, VDim> tmat = transform->GetMatrix();
    itk::Vector<double, VDim> toff = transform->GetOffset();
    Flip_LPS_to_RAS( amat, tmat, toff);
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


  symmmetric->SetInterpolator(  interpolator  );
  
  optimizer->SetScales( paramScales );
  symmmetric->SetTransformParameters( initialParameters );

  
  
  // Output some information about initial transform and parameters
  std::cout << "Initial Transform: " << std::endl;
  PrintMatrix( amat.GetVnlMatrix(), "%12.5f ", "   ");
  std::cout << "N = " << nParameters << std::endl << "params: " << symmmetric->GetTransformParameters() << std::endl;
  std::cout << "param scales are: " << optimizer->GetScales() << std::endl;

  optimizer->SetCostFunction( symmmetric );
  symmmetric->Initialize();
  optimizer->SetInitialPosition( symmmetric->GetTransformParameters() );

  typedef itk::Statistics::NormalVariateGenerator  GeneratorType;

  typename GeneratorType::Pointer generator = GeneratorType::New();


  generator->Initialize(12345);

  optimizer->SetNormalVariateGenerator( generator );
  optimizer->Initialize( atof(argv[8]), atof(argv[9]), atof(argv[10]) );
  optimizer->SetEpsilon( atof(argv[11]) );
  optimizer->SetMaximumIteration( atoi(argv[12]) );


  // Create the Command observer and register it with the optimizer.
  //
  typename CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
 
  optimizer->DebugOn();
  try 
    { 
    std::cout << "Registration starts!" << std::endl;
    optimizer->StartOptimization(); 
    std::cout << "Registration completed!" << std::endl;
    std::cout << "Optimizer stop condition: "
              << optimizer->GetStopConditionDescription()
              << std::endl;
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return EXIT_FAILURE;
    } 

  typename OptimizerType::ParametersType finalParameters = optimizer->GetCurrentPosition();
  
  
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

  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();

  itk::Matrix<double, VDim+1, VDim+1> finalmat;
  if(!strcmp(transform_name.c_str(),"rig"))
    {
    finaltransform = VersorRigid3DTransformType::New();
    finaltransform->SetParameters( finalParameters );
    itk::Matrix<double, VDim, VDim> tmat = finaltransform->GetMatrix();
    itk::Vector<double, VDim> toff = finaltransform->GetOffset();
    Flip_LPS_to_RAS( finalmat, tmat, toff);
    resample->SetTransform( finaltransform );
  
    }
  else
    {
    finaltransform = AffineTransformType::New();
    finaltransform->SetParameters( finalParameters );
    itk::Matrix<double, VDim, VDim> tmat = finaltransform->GetMatrix();
    itk::Vector<double, VDim> toff = finaltransform->GetOffset();
    Flip_LPS_to_RAS( finalmat, tmat, toff);

    resample->SetTransform( finaltransform );
    
    }

  std::cout << "Final Transform:" << std::endl;
  PrintMatrix( finalmat.GetVnlMatrix(), "%12.5f ", "   ");
  std::string outname( argv[3] );
  std::string outmatname = outname + ".mat";
  std::string outimname = outname + ".nii.gz";
  std::string outhwname = outname + "hwdef.nii.gz";

  ras_write( finalmat.GetVnlMatrix(), outmatname.c_str() );



  typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();


  resample->SetInput( symmmetric->GetMovingImage() );
  resample->SetInterpolator( interpolator );
  resample->SetDefaultPixelValue( 0.0 );
  resample->UseReferenceImageOn();
  resample->SetReferenceImage(symmmetric->GetFixedImage());


  // Write output image
  typedef itk::ImageFileWriter< FixedImageType > WriterType;
  typename WriterType::Pointer      writer =  WriterType::New();

  writer->SetFileName( outimname.c_str() );
  writer->SetInput( resample->GetOutput() );
  writer->Update();

  if (symmmetric->GetHalfwayImage())
    {
    // Write halfway image
    typename WriterType::Pointer      hwwriter =  WriterType::New();
    hwwriter->SetInput( symmmetric->GetHalfwayImage() );
    hwwriter->SetFileName( outhwname.c_str() );
    hwwriter->Update();
    }
  

  return EXIT_SUCCESS;
};

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << "outputNaming ";
    std::cerr << "InitialTransformFile(none for no initial) ";
    std::cerr << "RegType(rig or aff) ";
    std::cerr << "Metric ";
    std::cerr << "Symmetry ";
    std::cerr << "InitialRadius Grow Shrink Epsilon MaxIter "<< std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "fixedImageFile   : " << argv[1] << std::endl;
  std::cout << "movingImageFile  : " << argv[2] << std::endl;
  std::cout << "outputNaming     : " << argv[3] << std::endl;
  std::cout << "Initial Transform: " << argv[4] << std::endl;
  std::cout << "Registration Type: " << argv[5] << std::endl;
  std::cout << "Metric Type      : " << argv[6] << std::endl;
  std::cout << "Symmetry Type    : " << argv[7] << std::endl;
  std::cout << "Initial Radius   : " << argv[8] << std::endl;
  std::cout << "Growth Factor    : " << argv[9] << std::endl;
  std::cout << "Shrink Factor    : " << argv[10] << std::endl;
  std::cout << "Convergence Eps  : " << argv[11] << std::endl;
  std::cout << "Max Iteration    : " << argv[12] << std::endl;

    return EvolutionaryRegistration<double,3>(argc, argv);
  
}
