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

  if(argc < 3)
    {
    cout << "subfield_mask_slices reference_image output_image" << endl;
    return -1;
    }

  // Read reference segmentation
  ReaderType::Pointer fReader = ReaderType::New();
  fReader->SetFileName(argv[1]);
  fReader->Update();
  ImageType::Pointer ref = fReader->GetOutput();

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
  cout << "Slice Partition: ";
  size_t iLastHead = 0;
  for(size_t i = 0; i < n_slices; i++)
    {
    if(nTail[i] > nHead[i] + nBody[i])
      cout << "T";
    if(nHead[i] > nBody[i] + nTail[i])
      { cout << "H"; iLastHead = i; }
    if(nBody[i] > nHead[i] + nTail[i])
      cout << "B";
    else
      cout << ".";
    }
  cout << endl;

  // Flag slices as ERC or not
  vector<bool> iERC(n_slices,false);
  for(size_t i = iLastHead-1; i < iLastHead+2; i++)
    iERC[i] = true;

  // Report statistics
  cout << "                 ";
  for(size_t i = 0; i < n_slices; i++)
    {
    if(iERC[i])
      cout << "E";
    else
      cout << ".";
    }
  cout << endl;

  ref->DisconnectPipeline();
  for(IteratorType it(ref, ref->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    ImageType::IndexType idx = it.GetIndex();
    size_t i = idx[slice_dim];
    short label = 0;
    if(nTail[i] > nHead[i] + nBody[i])
      label = 3;
    else if(nHead[i] > nBody[i] + nTail[i])
      label = 1;
    else if(nBody[i] > nHead[i] + nTail[i])
      label = 2;
    
    if(iERC[i])
      label += 4;

    it.Set(label);
    }

  // Save the result
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fWriter = WriterType::New();
  fWriter->SetInput(ref);
  fWriter->SetFileName(argv[2]);
  fWriter->Update();
}
