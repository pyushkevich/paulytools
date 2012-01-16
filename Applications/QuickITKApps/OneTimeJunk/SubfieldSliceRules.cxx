#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <vector>
#include <set>
#include <map>
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

int usage()
{
  printf(
    "subfield_slice_rules: generate exclusion priors from segmentation and rules\n"
    "usage:\n"
    "  subfield_slice_rules input_seg.nii rules.txt output_pattern.nii\n"
    "parameters:\n"
    "  input_seg.nii        Input multi-label segmentation\n"
    "  rules.txt            Text file containing rules\n"
    "  output_pattern.nii   Output pattern (printf format, e.g. out_%%02d.nii.gz\n"
    "rules file format:\n"
    "  Each line in the rules file specifies exclusion rules in the form \n"
    "    X : Y\n"
    "  meaning 'if a slice contains any of the labels in the set X, it may \n"
    "  not contain any of the labels in the set Y' \n"
    "  For example: \n"
    "    5: 1 2 3 4\n"
    "    1 2 3 4: 5\n");
  return -1;
}

int main(int argc, char *argv[])
{
  typedef itk::Image<short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  if(argc < 4) return usage();

  // Read reference segmentation
  ReaderType::Pointer fReader = ReaderType::New();
  fReader->SetFileName(argv[1]);
  fReader->Update();
  ImageType::Pointer ref = fReader->GetOutput();

  // Read the rules files
  FILE *f = fopen(argv[2],"rt");
  if(!f)
    {
    cerr << "Can't open rules file" << argv[2] << endl;
    return -1;
    }

  typedef std::set<int> Exclusion;
  typedef std::map<int, Exclusion> Rules;
  Rules rules;
  Exclusion allRHS;

  char buffer[1024];
  while(fgets(buffer, 1024, f))
    {
    cout << "Read rule '" << buffer << endl;
    char *pColon = strchr(buffer,':');
    if(!pColon)
      {
      cerr << "Invalid rule!" << endl;
      return -1;
      }

    // Each entry on the left is treated as a separate rule
    *pColon = 0;

    // Read the right hand side of the rule
    Exclusion ruleRHS;
    char *pTok = strtok(pColon+1," \t\n,");
    while(pTok)
      {
      int rhs = atoi(pTok);
      if(rhs == 0)
        {
        cerr << "Bad rule RHS: " << pTok << endl;
        return -1;
        }
      ruleRHS.insert(rhs);
      allRHS.insert(rhs);
      pTok = strtok(NULL, " \t\n,");
      }
    
    pTok = strtok(buffer, " \t\n,");
    while(pTok)
      {
      int lhs = atoi(pTok);
      if(lhs == 0)
        {
        cerr << "Bad rule LHS: " << pTok << endl;
        return -1;
        }

      for(Exclusion::iterator it = ruleRHS.begin(); it!=ruleRHS.end(); ++it)
        {
        rules[lhs].insert(*it);
        }

      pTok = strtok(NULL, " \t\n,");
      }
    }
  fclose(f);

  // Create output images for all the RHS labels
  typedef std::map<int, ImageType::Pointer> PriorMap;
  PriorMap pm;

  for(Exclusion::iterator it = allRHS.begin(); it!=allRHS.end(); ++it)
    {
    pm[*it] = ImageType::New();
    pm[*it]->CopyInformation(ref);
    pm[*it]->SetRegions(ref->GetLargestPossibleRegion());
    pm[*it]->Allocate();
    pm[*it]->FillBuffer(0);
    }

  int slice_dim = 2;
  std::vector<Exclusion> sliceIdx(ref->GetBufferedRegion().GetSize()[slice_dim]);

  // For each z-slice, determine what unique voxels it has
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  for(IteratorType it(ref, ref->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    sliceIdx[it.GetIndex()[slice_dim]].insert(it.Get());
    }

  // For each slice, determine which labels are excluded
  std::vector<Exclusion> sliceExclude(sliceIdx.size());
  for(int k = 0; k < (int) sliceIdx.size(); k++)
    {
    for(Exclusion::iterator it = sliceIdx[k].begin(); it != sliceIdx[k].end(); ++it)
      {
      int lhs = *it;
      for(Exclusion::iterator qt = rules[lhs].begin(); qt!=rules[lhs].end(); ++qt)
        {
        sliceExclude[k].insert(*qt);
        }
      }
    }

  // Now paint all the slices
  for(IteratorType it(ref, ref->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    int k = it.GetIndex()[slice_dim];
    for(Exclusion::iterator qt = sliceExclude[k].begin(); qt != sliceExclude[k].end(); ++qt)
      {
      pm[*qt]->SetPixel(it.GetIndex(), 1);
      }
    }

  // Save the result
  for(Exclusion::iterator qt = allRHS.begin(); qt!=allRHS.end(); ++qt)
    {
    char fname[1024];
    int rhs = *qt;
    sprintf(fname, argv[3], rhs);

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer fWriter = WriterType::New();
    fWriter->SetInput(pm[rhs]);
    fWriter->SetFileName(fname);
    fWriter->Update();
    }
}
