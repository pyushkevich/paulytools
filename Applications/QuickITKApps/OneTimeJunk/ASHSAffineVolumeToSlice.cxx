#include "itkImage.h"
#include "itkImageFileReader.h"
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFactory.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkAffineTransform.h>
#include <iostream>
#include <cstdio>
#include <vnl/vnl_inverse.h>

using namespace std;

int usage()
{
  cout <<
    "ashs_affine_vol2slice: convert 3D affine transform to 2D\n"
    "usage: \n"
    "  ashs_affine_vol2slice atlas_tse.nii target_tse.nii atlas_seg.nii targ2atlas.txt range outname\n"
    "parameters: \n"
    "  atlas_tse.nii       Atlas TSE image\n"
    "  target_tse.nii      TSE image to segment\n"
    "  atlas_seg.nii       Binary segmentation of atlas\n"
    "  targ2atlas.txt      ITK transform (atlas as fixed, target as moving) \n"
    "  range               How many slices to try on both sides\n"
    "  outname             Naming prefix\n";
  return -1;
}

int main(int argc, char *argv[])
{
  typedef itk::Image<float,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::MatrixOffsetTransformBase<double, 3, 3> MOTBType;
  typedef itk::AffineTransform<double, 3> AffTran;

  if(argc < 7)
    return usage();

  const char *fnAtlas      = argv[1];
  const char *fnTarget     = argv[2];
  const char *fnAtlasSeg   = argv[3];
  const char *fnAffine     = argv[4];
  int range = atoi(argv[5]);
  const char *fnOutput     = argv[6];

  // Load the atlas image
  ReaderType::Pointer readAtlas = ReaderType::New();
  readAtlas->SetFileName(fnAtlas);
  readAtlas->Update();
  ImageType::Pointer ia = readAtlas->GetOutput();

  // Load the target image
  ReaderType::Pointer readTarget = ReaderType::New();
  readTarget->SetFileName(fnTarget);
  readTarget->Update();
  ImageType::Pointer it = readTarget->GetOutput();

  // Load the atlas image
  ReaderType::Pointer readAtlasSeg = ReaderType::New();
  readAtlasSeg->SetFileName(fnAtlasSeg);
  readAtlasSeg->Update();
  ImageType::Pointer iaseg = readAtlasSeg->GetOutput();

  // Read the transform
  itk::TransformFactory<MOTBType>::RegisterTransform();
  itk::TransformFactory<AffTran>::RegisterTransform();
  itk::TransformFileReader::Pointer readTransform = itk::TransformFileReader::New();
  readTransform->SetFileName(fnAffine);
  readTransform->Update();

  itk::TransformBase *base = readTransform->GetTransformList()->front();
  MOTBType *motb = dynamic_cast<MOTBType *>(base);

  // Get VOX/VOX transform from reference to atlas, i.e.
  // Q: K_atlas = Q * K_target
  typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
  MatrixType dca, dct, vxa, vxt, tna, tnt, Ma, Mt, A;
  dca.set_identity(); vxa.set_identity(); tna.set_identity();
  dct.set_identity(); vxt.set_identity(); tnt.set_identity();
  A.set_identity();
  for (int i=0; i < 3; i++)
    {
    for (int j=0; j < 3; j++)
      {
      dca[i][j] = ia->GetDirection()[i][j];
      dct[i][j] = it->GetDirection()[i][j];
      A[i][j] = motb->GetMatrix()[i][j];
      }
    vxa[i][i] = ia->GetSpacing()[i];
    vxt[i][i] = it->GetSpacing()[i];
    tna[i][3] = ia->GetOrigin()[i];
    tnt[i][3] = it->GetOrigin()[i];
    dca[3][i] = 0; vxa[3][i] = 0; tna[3][i] = 0;
    dct[3][i] = 0; vxt[3][i] = 0; tnt[3][i] = 0;
    A[i][3] = motb->GetOffset()[i];
    }

  // These are the voxel to LPS mappings
  Ma = tna * dca * vxa;
  Mt = tnt * dct * vxt;
  
  // The Q matrix
  MatrixType Q = vnl_inverse(Ma) * vnl_inverse(A) * Mt;
  MatrixType Qi = vnl_inverse(Q);
  cout << "Q = " << Q << endl;

  // For each slice in the atlas, find matching slice in the target
  int nsa = ia->GetBufferedRegion().GetSize()[2];
  int nst = it->GetBufferedRegion().GetSize()[2];

  // Accumulator
  std::vector<double> tslice(nsa, 0.0), nvox(nsa, 0.0);
  typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
  ItType ps(iaseg, iaseg->GetBufferedRegion());
  for(; !ps.IsAtEnd(); ++ps)
    {
    int islice = ps.GetIndex()[2];
    if(ps.Get() != 0.0)
      {
      vnl_vector_fixed<double, 4> vx, vxt;
      vx[0] = ps.GetIndex()[0]; vx[1] = ps.GetIndex()[1]; vx[2] = ps.GetIndex()[2];
      vx[3] = 1.0;
      vxt = Qi * vx;

      tslice[islice] += vxt[2];
      nvox[islice]++;
      }
    }

  // Now, associate a range of atlas slices with each target slice
  vector<int> asfirst(nst, -1);
  vector<int> aslast(nst, -1);

  for (int i = 0; i < nsa; i++)
    {
    int tscenter = (int)(0.5 + tslice[i] / nvox[i]);
    for (int t = std::max(0, tscenter-range); 
      t <= std::min(nst-1, tscenter+range); t++)
      {
      if (asfirst[t] < 0 || asfirst[t] > i)
        asfirst[t] = i;
      if (aslast[t] < 0 || aslast[t] < i)
        aslast[t] = i;
      }
    }

  // Now, for each target slice, compute 2D affine transform from each
  // atlas to the target
  for (int t = 0; t < nst; t++)
    {
    for (int a = asfirst[t]; a <= aslast[t]; a++)
      {
      if(a < 0) continue;

      // Get origins of slices
      vnl_vector_fixed<double, 4> vsa, vst, soa, sot;
      vsa[0] = 0.0; vsa[1] = 0.0; vsa[2] = a; vsa[3] = 1.0;
      vst[0] = 0.0; vst[1] = 0.0; vst[2] = t; vst[3] = 1.0;
      soa = Ma * vsa; sot = Mt * vst;

      typedef vnl_matrix_fixed<double, 3, 3> MatrixType2D;
      MatrixType2D Ga, Gt, R, Z;
      Ga.set_identity(); Gt.set_identity(); R.set_identity();
      for(size_t q = 0; q < 2; q++)
        {
        Ga[q][q] = vxa[q][q]; Ga[q][2] = soa[q];
        Gt[q][q] = vxt[q][q]; Gt[q][2] = sot[q];
        R[q][0] = Q[q][0]; R[q][1] = Q[q][1]; R[q][2] = Q[q][2] * t + Q[q][3];
        }

      // This is the matrix between physical spaces of the two slices
      cout << "R(" << t << "," << a << ") = " << R << endl;
      Z = Ga * R * vnl_inverse(Gt);

      // Create an ITK affine transform
      typedef itk::MatrixOffsetTransformBase<double, 2, 2> MOTBType2D;
      MOTBType2D::Pointer atran = MOTBType2D::New();

      // Populate its matrix
      MOTBType2D::MatrixType amat = atran->GetMatrix();
      MOTBType2D::OffsetType aoff = atran->GetOffset();
      for(size_t r = 0; r < 2; r++)
        {
        for(size_t c = 0; c < 2; c++)
          {
          amat(r,c) = Z(r,c);
          }
        aoff[r] = Z(r,2);
        }

      atran->SetMatrix(amat);
      atran->SetOffset(aoff);

      // Generate filename
      char buffer[1024];
      sprintf(buffer,"%s_ts%02d_as%02d_itktran.txt", fnOutput, t, a);

      // Write the transform
      itk::TransformFileWriter::Pointer wrt = itk::TransformFileWriter::New();
      wrt->SetInput(atran);
      wrt->SetFileName(buffer);
      wrt->Update();
      }
    }
}
