#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <vector>
#include <iostream>

using namespace std;

// Dimension of the slices
static const unsigned int slice_dim = 2;

// Mapper functions for thresholding body, head, tail
template<class TInput, class TOutput> class PickFunctor 
{
public:
  TInput x;
  TOutput operator()(const TInput &in)
    { if(in == x) return 1; else return 0; }
  bool operator == (const PickFunctor<TInput, TOutput> &dummy)
    { return dummy.x == x; }
  bool operator != (const PickFunctor<TInput, TOutput> &dummy)
    { return dummy.x != x; }
};

int main(int argc, char *argv[])
{
  typedef itk::Image<short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  if(argc < 4)
    {
    cout << "subfield_leveler reference_image source_image output_image" << endl;
    return -1;
    }

  // Read reference segmentation
  ReaderType::Pointer fReader = ReaderType::New();
  fReader->SetFileName(argv[1]);
  fReader->Update();
  ImageType::Pointer ref = fReader->GetOutput();

  // Read target segmentation
  ReaderType::Pointer fReader2 = ReaderType::New();
  fReader2->SetFileName(argv[2]);
  fReader2->Update();
  ImageType::Pointer src = fReader2->GetOutput();

  // For each z-slice, find the number of head, tail and body-specific voxels
  size_t n_slices = ref->GetBufferedRegion().GetSize()[slice_dim];
  vector<long> nHead(n_slices), nTail(n_slices), nBody(n_slices);

  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  for(IteratorType it(ref, ref->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    ImageType::IndexType idx = it.GetIndex();
    short val = it.Value();
    if(val == 5)
      nHead[idx[slice_dim]]++;
    else if(val == 6)
      nTail[idx[slice_dim]]++;
    else if(val == 1 || val == 2 || val == 3 || val == 4)
      nBody[idx[slice_dim]]++;
    }

  // Report statistics
  for(size_t i = 0; i < n_slices; i++)
    {
    cout << "Slice " << i << " belongs to ";
    if(nTail[i] > nHead[i] + nBody[i])
      cout << "TAIL";
    if(nHead[i] > nBody[i] + nTail[i])
      cout << "HEAD";
    if(nBody[i] > nHead[i] + nTail[i])
      cout << "BODY";
    cout << endl;
    }

  // Compute the distance transforms from the source segmentation
  typedef itk::Image<float, 3> FloatImageType;
  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, FloatImageType> DistMapType;
  typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, PickFunctor<short, short> > PickerType;

  // Set up distance maps to each label
  DistMapType::Pointer fDist[11];
  for(size_t i = 0; i < 11; i++)
    {
    PickFunctor<short,short> pf; pf.x = i;
    PickerType::Pointer fPicker = PickerType::New();
    fPicker->SetFunctor(pf);
    fPicker->SetInput(src);
    fDist[i] = DistMapType::New();
    fDist[i]->SetUseImageSpacing(true);
    fDist[i]->SetSquaredDistance(true);
    fDist[i]->SetInput(fPicker->GetOutput());
    fDist[i]->Update();
    }

  // Iterate through the source image
  src->DisconnectPipeline();
  for(IteratorType it(src, src->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    ImageType::IndexType idx = it.GetIndex();
    size_t i = idx[slice_dim];
    short v = it.Value();
    size_t cand_tail[] = {6, 8, 9, 10};
    size_t cand_head[] = {5, 8, 9, 10};
    size_t cand_body[] = {1, 2, 3, 4, 8, 9, 10};
    size_t *mycand = NULL, n_mycand = 0;

    if(nTail[i] > nHead[i] + nBody[i])
      {
      if(v==1 || v==2 || v==3 || v==4 || v==5 || v==7)
        { mycand = cand_tail; n_mycand = 4; cout << "Voxel " << idx << " label " << v << " is in the tail ";}
      }
    if(nHead[i] > nBody[i] + nTail[i])
      {
      if(v==1 || v==2 || v==3 || v==4 || v==6 || v==7)
        { mycand = cand_head; n_mycand = 4; cout << "Voxel " << idx << " label " << v << " is in the head ";}
      }
    if(nBody[i] > nHead[i] + nTail[i])
      {
      if(v==5 || v==6 || v==7)
        { mycand = cand_body; n_mycand = 7; cout << "Voxel " << idx << " label " << v << " is in the body ";}
      }
    
    if(n_mycand)
      {
      double best_dist = 1e100;
      size_t closest = 0;
      for(size_t j = 0; j < n_mycand; j++)
        {
        size_t c = mycand[j];
        double d = fDist[c]->GetOutput()->GetPixel(idx);
        if(d < best_dist)
          { best_dist = d; closest = c; }
        }
      it.Set(closest);
      cout << ", became " << closest << endl;
      }
    }

  // Save the result
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fWriter = WriterType::New();
  fWriter->SetInput(src);
  fWriter->SetFileName(argv[3]);
  fWriter->Update();
}
