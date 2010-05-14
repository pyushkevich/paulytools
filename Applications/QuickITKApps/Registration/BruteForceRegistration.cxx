#include "itkOrientedRASImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInterpolateImageFunction.h"
#include "itkImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkExhaustiveOptimizer.h"
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

using namespace std;
using namespace itk;

enum InitType { UNCHANGED, FULLRAND, IDENTITY, ROTRAND };

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
  CommandIterationUpdate() { m_BestMetricValue = 0.0; m_CurrentIteration = 0; };
public:
  typedef itk::ExhaustiveOptimizer     OptimizerType;
  typedef   const OptimizerType *                  OptimizerPointer;
  int    m_CurrentIteration;

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
      double currentValue = optimizer->GetCurrentValue();
      m_CurrentIteration++;
      std::cout << m_CurrentIteration << "   ";
      std::cout << currentValue << "   ";
      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
private:
  double m_BestMetricValue;
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
  fout.precision(12);
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
ReadTransformFile(typename TAffine::Pointer atran, char * fn, InitType init, bool isrigid)
{
  // Read transform based on type
  if ( init == UNCHANGED || init == ROTRAND )
  {
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
      itk::Vector<double, VDim> aoff, rtran ;
  
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
      if (atran->GetDebug())
        {
        //rtran[0]=1.5;rtran[1]=2.5;rtran[2]=3.5;
        //cout << "calling transform print after setting random tranalation"<< endl;
        //atran->SetTranslation( rtran );
        //atran->Print( std::cout );
        typename TAffine::InputPointType rcenter;
        //rcenter[0]=11.5;rcenter[1]=0.5;rcenter[2]=7.5;
        //atran->SetCenter( rcenter );
        //cout << "calling transform print after setting random center"<< endl;
        //atran->Print( std::cout );
        }
        atran->SetMatrix(amat);
      if (atran->GetDebug())
        {
        cout << "calling transform print after setting matrix"<< endl;
        atran->Print( std::cout );
        }
      atran->SetOffset(aoff);
      if (atran->GetDebug())
        {
        cout << "calling transform print after setting offset"<< endl;
        atran->Print( std::cout );
        }
      }
  }
  if ( init == IDENTITY )
    {
    atran->SetIdentity();
    return;
    }
  if ( init == UNCHANGED )
    return;
  if ( init == ROTRAND )
    {
    if ( isrigid )
      {
      typename TAffine::ParametersType params = atran->GetParameters();
      typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
      GeneratorType::Pointer      generator = GeneratorType::New();
      generator->Initialize();
      // Max rotation -alpha to alpha degrees
      double MAXANGLE = 20.0*( vnl_math::pi )/180.0;
      double randomangle = generator->GetVariateWithClosedRange ( MAXANGLE ) - MAXANGLE/2.0;
      double randomXaxis = generator->GetVariateWithClosedRange ( 1.0 );
      double randomYaxis = generator->GetVariateWithClosedRange ( 1.0 );
      double randomZaxis = generator->GetVariateWithClosedRange ( 1.0 );
      double mag = vcl_sqrt( randomXaxis*randomXaxis + randomYaxis*randomYaxis + randomZaxis*randomZaxis );
      randomXaxis = randomXaxis/mag;
      randomYaxis = randomYaxis/mag;
      randomZaxis = randomZaxis/mag;
      double q0 = randomXaxis * vcl_sin( randomangle/2 );
      double q1 = randomYaxis * vcl_sin( randomangle/2 );
      double q2 = randomZaxis * vcl_sin( randomangle/2 );
      //double q3 = vcl_cos( randomangle/2 );

      params[0] = q0;
      params[1] = q1;
      params[2] = q2;

      std::cout << "Random angle: " << randomangle*180.0/vnl_math::pi << " degrees" << std::endl;
      std::cout << "Random axis: [ " << randomXaxis << " , " << randomYaxis << " , " << randomZaxis << " ]" << std::endl;
      std::cout << "Random rotation parameters: " << params << std::endl;
      atran->SetParameters( params );
      
      }
    else
      {
      // TODO randomly initialize the rotation part of affine parameters
      }
    return;
    }

  if ( init == FULLRAND )
    {  
    if ( isrigid )
      {
      typename TAffine::ParametersType params = atran->GetParameters();
      typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
      GeneratorType::Pointer      generator = GeneratorType::New();
      generator->Initialize();
      // Max rotation -5 to 5 degrees
      double MAXANGLE = 10.0*( vnl_math::pi )/180.0;
      double randomangle = generator->GetVariateWithClosedRange ( MAXANGLE ) - MAXANGLE/2.0;
      double randomXaxis = generator->GetVariateWithClosedRange ( 1.0 );
      double randomYaxis = generator->GetVariateWithClosedRange ( 1.0 );
      double randomZaxis = generator->GetVariateWithClosedRange ( 1.0 );
      double mag = vcl_sqrt( randomXaxis*randomXaxis + randomYaxis*randomYaxis + randomZaxis*randomZaxis );
      randomXaxis = randomXaxis/mag;
      randomYaxis = randomYaxis/mag;
      randomZaxis = randomZaxis/mag;
      double q0 = randomXaxis * vcl_sin( randomangle/2 );
      double q1 = randomYaxis * vcl_sin( randomangle/2 );
      double q2 = randomZaxis * vcl_sin( randomangle/2 );
      //double q3 = vcl_cos( randomangle/2 );

      // Max translation -20 to 20 mm
      double MAXDISP = 40.0;
      double randomXtrans = generator->GetVariateWithClosedRange ( MAXDISP ) - MAXDISP/2.0;
      double randomYtrans = generator->GetVariateWithClosedRange ( MAXDISP ) - MAXDISP/2.0;
      double randomZtrans = generator->GetVariateWithClosedRange ( MAXDISP ) - MAXDISP/2.0;

      params[0] = q0;
      params[1] = q1;
      params[2] = q2;
      params[3] = randomXtrans;
      params[4] = randomYtrans;
      params[5] = randomZtrans;

      atran->SetParameters( params );

      }   
    else
      {
      // TODO randomly initialize 12 affine parameters
      }
    return;
    }

  

}

template <class TPixel, unsigned int VDim>
int
BruteForceRegistration(int argc, char * argv[])
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
  typedef itk::ExhaustiveOptimizer           OptimizerType;
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
  typename FixedImageReaderType::Pointer  fixedMaskImageReader  = FixedImageReaderType::New();
  typename MovingImageReaderType::Pointer movingMaskImageReader = MovingImageReaderType::New();

  // Set up fixed and moving images
  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  //symmmetric->DebugOn();
  symmmetric->SetFixedImage(    fixedImageReader->GetOutput()    );
  symmmetric->SetMovingImage(   movingImageReader->GetOutput()   );
  
  fixedImageReader->Update();
  movingImageReader->Update();

  symmmetric->SetFixedImageRegion( 
       fixedImageReader->GetOutput()->GetBufferedRegion() );

  // Set up mask images
  std::string fixedmaskfn( argv[3] );
  std::string movingmaskfn( argv[4] );
  if (fixedmaskfn.compare(fixedmaskfn.size()-4,4,"none") != 0)
    {
    fixedMaskImageReader->SetFileName(  argv[3] );
    symmmetric->SetInternalFixedImageMask(    fixedMaskImageReader->GetOutput()    );
    fixedMaskImageReader->Update();
    }
  if (movingmaskfn.compare(movingmaskfn.size()-4,4,"none") != 0)
    {
    movingMaskImageReader->SetFileName(  argv[4] );
    symmmetric->SetInternalMovingImageMask(    movingMaskImageReader->GetOutput()    );
    movingMaskImageReader->Update();
    }
    


  // Get spacing and size for parameter scaling later
  typename FixedImageType::SizeType fixedsize = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize(); 
  typename FixedImageType::SpacingType fixedspacing = fixedImageReader->GetOutput()->GetSpacing();
  double minspacing = fixedspacing[0];
  for (size_t i=0; i< VDim; i++)
    if ( fixedspacing[i] < minspacing )
      minspacing = fixedspacing[i];


  // Create the appropriate metric
  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType> MetricType;
  typename MetricType::Pointer metric;


  std::string reg_symm( argv[10] );
  if(!strcmp(reg_symm.c_str(),"symm"))
    symmmetric->UseSymmetricOn();
  else if (!strcmp(reg_symm.c_str(),"asym"))
    symmmetric->UseSymmetricOff();
 else
    {
    cerr << "Unknown symmetry type " << reg_symm << endl;
    return EXIT_FAILURE;
    };

  // Use slave metric or not 
  std::string masterslave( argv[11] );
  if(!strcmp(masterslave.c_str(),"master"))
    symmmetric->UseSlaveMetricOff();
  else if (!strcmp(masterslave.c_str(),"slave"))
    symmmetric->UseSlaveMetricOn();
 else
    {
    cerr << "Unknown master/slave type " << masterslave << endl;
    return EXIT_FAILURE;
    };

    

  std::string metric_name( argv[9] );
  double isMaximize = true;

  if(!strcmp(metric_name.c_str(),"mutualinfo"))
    {
    typename MattesMutualInformationImageToImageMetricType::Pointer 
      typedmetric = MattesMutualInformationImageToImageMetricType::New();
    //typedmetric->SetNumberOfHistogramBins( 40 );
    //typedmetric->SetNumberOfSpatialSamples( 10000 );
    metric = typedmetric;
    symmmetric->SetAsymMetric( metric  );
    isMaximize = false;
    }
  else if(!strcmp(metric_name.c_str(),"normmi"))
    {
    typename NormalizedMutualInformationHistogramImageToImageMetricType::Pointer 
      typedmetric = NormalizedMutualInformationHistogramImageToImageMetricType::New();
    metric = typedmetric;
    symmmetric->SetAsymMetric( metric  );
    isMaximize = true;
    }
  else if(!strcmp(metric_name.c_str(),"leastsq"))
    {
    typename MeanSquaresImageToImageMetricType::Pointer 
      typedmetric = MeanSquaresImageToImageMetricType::New();
    metric = typedmetric;
    symmmetric->SetAsymMetric( metric  );
    isMaximize = false;
    }
  else if(!strcmp(metric_name.c_str(),"normcorr"))
    {
    typename NormalizedCorrelationImageToImageMetricType::Pointer 
      typedmetric = NormalizedCorrelationImageToImageMetricType::New();
    metric = typedmetric;
    symmmetric->SetAsymMetric( metric  );
    isMaximize = false;
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
  
  std::string fn( argv[6] );
  bool isinit = (fn.compare(fn.size()-4,4,"none") != 0);
  std::string inittypestr( argv[7] );
  InitType inittype;
  if (!strcmp(inittypestr.c_str(),"unchanged"))
    inittype = UNCHANGED;
  else if (!strcmp(inittypestr.c_str(),"fullrand"))
    inittype = FULLRAND;
  else if (!strcmp(inittypestr.c_str(),"rotrand"))
    inittype = ROTRAND;
  else
    {
    cerr << "Unknown init type " << inittypestr << endl;
    return EXIT_FAILURE;
    };
    


  std::string transform_name( argv[8] );

  int nParameters;
  itk::Matrix<double, VDim+1, VDim+1> amat;

  typename RegistrationType::ParametersType initialParameters;
  typename OptimizerType::ScalesType paramScales;

  if(!strcmp(transform_name.c_str(),"rig"))
    {
    transform = VersorRigid3DTransformType::New();
    transform->SetIdentity();
    transform->DebugOn();
    if (isinit)
      {
      ReadTransformFile<MOTBType, VDim>( transform, argv[6], inittype, true);
      //cout << "Transform file " << fn << " read" << endl;
      //std::cout << "Just after reading params: " << transform->GetParameters() << std::endl;
      }
    else
      cout << "Initial Transform set to identity" << endl;
    symmmetric->SetTransform(     transform     );
    nParameters = transform->GetNumberOfParameters();
    initialParameters.SetSize( nParameters );
    paramScales.SetSize( nParameters );
    //transform->DebugOn();
    initialParameters = transform->GetParameters();
    itk::Matrix<double, VDim, VDim> tmat = transform->GetMatrix();
    itk::Vector<double, VDim> toff = transform->GetOffset();
    if (transform->GetDebug())
      {
      cout << "calling transform print after setting all"<< endl;
      transform->Print( std::cout );
      }
    Flip_LPS_to_RAS( amat, tmat, toff);
    for (size_t i=0; i < VDim; i++)
      paramScales[i] = 1.0;
    // 1 degree of rotation is set equivalent to 1 voxel translation
    for (size_t i=0; i < VDim ; i++)
      //paramScales[VDim + i] = fixedspacing[i] * (180/3.14159);
      //paramScales[VDim + i] = (100*fixedspacing[i]) ;
      paramScales[VDim + i] = 100*minspacing;

    }
  else if(!strcmp(transform_name.c_str(),"aff"))
    {
    transform = AffineTransformType::New();
    //transform->DebugOn();
    transform->SetIdentity();
    if (isinit)
      {
      ReadTransformFile<MOTBType, VDim>( transform, argv[6], inittype, false);
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
    if (transform->GetDebug())
      {
      cout << "calling transform print after setting all"<< endl;
      transform->Print( std::cout );
      }
    Flip_LPS_to_RAS( amat, tmat, toff);
    paramScales.fill(1.0);
    for (size_t i=0; i < VDim; i++)
      paramScales[i*VDim + i] = 1.0;
    // 1 degree of rotation is set equivalent to 1 voxel translation 
    // Make sure because here the parameters are not quarternion
    for (size_t i=0; i < VDim ; i++)
      //paramScales[VDim*VDim + i] = fixedspacing[i] * (180/3.14159);
      //paramScales[VDim*VDim + i] = (100*fixedspacing[i]);
      paramScales[VDim*VDim + i] = 100*minspacing;

    }
  else 
    {
    cerr << "Unknown transform type " << transform_name << endl;
    return EXIT_FAILURE;
    }; 


  // Create the Command observer and register it with the optimizer.
  //
  typename CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
 
  optimizer->DebugOn();
  typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  typename InterpolatorType::Pointer   fixedinterpolator  = InterpolatorType::New();

  symmmetric->SetInterpolator(  interpolator  );
  symmmetric->SetFixedInterpolator(  fixedinterpolator  );
  symmmetric->SetTransformParameters( initialParameters );
  symmmetric->Initialize();

  optimizer->SetCostFunction( symmmetric );
  
  double steplength = atof(argv[14] );
  OptimizerType::StepsType nsteps( symmmetric->GetNumberOfParameters() );
  int tsteps[3];
  for ( size_t i = 0; i < symmmetric->GetNumberOfParameters(); i++ )
    nsteps[i] = atoi(argv[12]);
  // Last 3 parameters are translation, setting steps for these
  for ( size_t i = symmmetric->GetNumberOfParameters() - 3; i < symmmetric->GetNumberOfParameters(); i++ )
    tsteps[i-3] = atoi(argv[13]);
  for (size_t i=symmmetric->GetNumberOfParameters()-3; i<symmmetric->GetNumberOfParameters(); i++)
    {
    nsteps[i] = tsteps[i - symmmetric->GetNumberOfParameters() + 3];
    }

  optimizer->SetInitialPosition( initialParameters );
  optimizer->SetScales( paramScales );
  
  transform->SetParameters( initialParameters );
  itk::Matrix<double, VDim, VDim> bmat = transform->GetMatrix();
  itk::Vector<double, VDim> boff = transform->GetOffset();
  if (transform->GetDebug())
    {
    cout << "calling transform print after setting optimizer initial condition"<< endl;
    transform->Print( std::cout );
    }
  Flip_LPS_to_RAS( amat, bmat, boff);
  // Output some information about initial transform and parameters
  std::cout << "Initial Transform: " << std::endl;
  PrintMatrix( amat.GetVnlMatrix(), "%12.5f ", "   ");
  std::cout << "N = " << nParameters << std::endl << "params: " << optimizer->GetInitialPosition() << std::endl;
  std::cout << "fixed params: " << transform->GetFixedParameters() << std::endl;
  std::cout << "param scales are: " << optimizer->GetScales() << std::endl;

  optimizer->SetStepLength( steplength );
  optimizer->SetNumberOfSteps( nsteps );

  optimizer->Print(std::cerr);
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

  double bestValue;
  typename OptimizerType::ParametersType finalParameters;
  if ( isMaximize )
    {
    bestValue = optimizer->GetMaximumMetricValue();
    finalParameters = optimizer->GetMaximumMetricValuePosition();
    }
  else
    {
    bestValue = optimizer->GetMinimumMetricValue();
    finalParameters = optimizer->GetMinimumMetricValuePosition();
    }
  
  unsigned int numberOfIterations = observer->m_CurrentIteration;
  
  

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
  resample->DebugOn();
  resample->SetInput(  movingImageReader->GetOutput() );

  itk::Matrix<double, VDim+1, VDim+1> finalmat;
  itk::Matrix<double, VDim, VDim> tmat;
  itk::Vector<double, VDim> toff;
  
  if(!strcmp(transform_name.c_str(),"rig"))
    {
    finaltransform = VersorRigid3DTransformType::New();
    finaltransform->SetParameters( finalParameters );
    finaltransform->SetFixedParameters( transform->GetFixedParameters() );
    tmat = finaltransform->GetMatrix();
    toff = finaltransform->GetOffset();
    Flip_LPS_to_RAS( finalmat, tmat, toff);
    resample->SetTransform( finaltransform );
  
    }
  else
    {
    finaltransform = AffineTransformType::New();
    finaltransform->SetParameters( finalParameters );
    finaltransform->SetFixedParameters( transform->GetFixedParameters() );
    tmat = finaltransform->GetMatrix();
    toff = finaltransform->GetOffset();
    Flip_LPS_to_RAS( finalmat, tmat, toff);

    resample->SetTransform( finaltransform );
    
    }

  std::cout << "Final tmat: " << tmat << std::endl;
  std::cout << "Final toff: " << toff << std::endl;
  std::cout << "Final Paramaters from final transform: " << finaltransform->GetParameters() <<  std::endl;
  std::cout << "Final Transform:" << std::endl;
  PrintMatrix( finalmat.GetVnlMatrix(), "%12.5f ", "   ");
  std::string outname( argv[5] );
  std::string outmatname = outname + ".mat";
  std::string outimname = outname + ".nii.gz";
  std::string outhwname = outname + "hwdef.nii.gz";

  ras_write( finalmat.GetVnlMatrix(), outmatname.c_str() );



  typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();


  resample->SetDefaultPixelValue( 0.0 );
  resample->SetInterpolator( interpolator );
  resample->UseReferenceImageOn();
  resample->SetReferenceImage( fixedImageReader->GetOutput() );

  std::cout << "Resamplefilter transform parameters: " << resample->GetTransform()->GetParameters() << endl;
  std::cout << "Resamplefilter transform fixed parameters: " << resample->GetTransform()->GetFixedParameters() << endl;

  resample->Update();
  // Write output image
  typedef itk::ImageFileWriter< FixedImageType > WriterType;
  typename WriterType::Pointer      writer =  WriterType::New();

  writer->SetFileName( outimname.c_str() );
  writer->SetInput( resample->GetOutput() );
  writer->Update();

  // Debugging why metric is poor outside, although optimization result looks good
  
  typename MOTBType::Pointer idtran;
  idtran = VersorRigid3DTransformType::New();
  idtran->SetIdentity();

  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType> MetricType;
  typename MetricType::Pointer finalmetric;

  if(!strcmp(metric_name.c_str(),"mutualinfo"))
    {
    typename MattesMutualInformationImageToImageMetricType::Pointer
      typedmetric = MattesMutualInformationImageToImageMetricType::New();
    finalmetric = typedmetric;
    }
  else if(!strcmp(metric_name.c_str(),"normmi"))
    {
    typename NormalizedMutualInformationHistogramImageToImageMetricType::Pointer
      typedmetric = NormalizedMutualInformationHistogramImageToImageMetricType::New();
    finalmetric = typedmetric;
    }
  else if(!strcmp(metric_name.c_str(),"leastsq"))
    {
    typename MeanSquaresImageToImageMetricType::Pointer
      typedmetric = MeanSquaresImageToImageMetricType::New();
    finalmetric = typedmetric;
    }
  else if(!strcmp(metric_name.c_str(),"normcorr"))
    {
    typename NormalizedCorrelationImageToImageMetricType::Pointer
      typedmetric = NormalizedCorrelationImageToImageMetricType::New();
    finalmetric = typedmetric;
    }
  else
    {
    cerr << "Unknown metric type " << metric_name << endl;
    return EXIT_FAILURE;
    };


  finalmetric->DebugOn();
  finalmetric->SetFixedImage( fixedImageReader->GetOutput() );
  finalmetric->SetMovingImage( resample->GetOutput() );
  finalmetric->SetFixedImageRegion(
       fixedImageReader->GetOutput()->GetBufferedRegion() );
  finalmetric->SetInterpolator( interpolator );
  finalmetric->SetTransform( idtran );
  finalmetric->Initialize();
  std::cout << "metric transform parameters with identity: " << finalmetric->GetTransform()->GetParameters() << endl;
  std::cout << "metric transform fixed parameters with identity: " << finalmetric->GetTransform()->GetFixedParameters() << endl;
  cout << "Calculating final metric outside of optimization with already resliced image: " << finalmetric->GetValue( idtran->GetParameters() ) <<endl;
  finalmetric->SetMovingImage( movingImageReader->GetOutput() );
  finalmetric->SetTransform( finaltransform );
  finalmetric->Initialize();
  std::cout << "finalmetric transform parameters: " << finalmetric->GetTransform()->GetParameters() << endl;
  std::cout << "finalmetric transform fixed parameters : " << finalmetric->GetTransform()->GetFixedParameters() << endl;
  cout << "Calculating final metric outside of optimization with the transform applied by metric: " << finalmetric->GetValue( finaltransform->GetParameters() ) <<endl;

  // Debug -- write out the fixed and moving images used in the metric -- can we then find out why the mtric is different here than outside ?
/*
    typename WriterType::Pointer      xwriter =  WriterType::New();
    xwriter->SetInput( fixedImageReader->GetOutput() );
    xwriter->SetFileName( "ibn_work_L/evolreg/fixedinmetric.nii.gz" );
    xwriter->Update();
    typename WriterType::Pointer      ywriter =  WriterType::New();
    ywriter->SetInput(  movingImageReader->GetOutput() );
    ywriter->SetFileName( "ibn_work_L/evolreg/movinginmetric.nii.gz" );
    ywriter->Update();
  */

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
    std::cerr << " fixedImageMask  movingImageMask ";
    std::cerr << "outputNaming ";
    std::cerr << "InitialTransformFile(none for no initial) ";
    std::cerr << "InitializationType(unchanged, fullrand or rotrand) ";
    std::cerr << "RegType(rig or aff) ";
    std::cerr << "Metric ";
    std::cerr << "Symmetry ";
    std::cerr << "Master/Slave ";
    std::cerr << "NoOfAngleSteps NoOfTranslationSteps StepLength"<< std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "fixedImageFile   : " << argv[1] << std::endl;
  std::cout << "movingImageFile  : " << argv[2] << std::endl;
  std::cout << "fixedImageMask   : " << argv[3] << std::endl;
  std::cout << "movingImageMask  : " << argv[4] << std::endl;
  std::cout << "outputNaming     : " << argv[5] << std::endl;
  std::cout << "Initial Transform: " << argv[6] << std::endl;
  std::cout << "Init Type        : " << argv[7] << std::endl;
  std::cout << "Registration Type: " << argv[8] << std::endl;
  std::cout << "Metric Type      : " << argv[9] << std::endl;
  std::cout << "Symmetry Type    : " << argv[10] << std::endl;
  std::cout << "Master/Slave     : " << argv[11] << std::endl;
  std::cout << "#anglesteps      : " << argv[12] << std::endl;
  std::cout << "#translatesteps  : " << argv[13] << std::endl;
  std::cout << "Step Length      : " << argv[14] << std::endl;

    return BruteForceRegistration<double,3>(argc, argv);
  
}
