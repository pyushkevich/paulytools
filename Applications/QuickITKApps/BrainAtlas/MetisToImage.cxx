/**
 * Read metis output, create an image
 */
#include "ReadWriteImage.h"
#include <iostream>

using namespace std;


int main(int argc, char *argv[])
{
  if(argc != 5) 
    {
    cerr << "usage: imginput.hdr metis2image metisout.txt metistemp.txt imageoutput.hdr" << endl;
    return -1;
    }

  // Read the input image image
  typedef itk::Image<short,3> ImageType;
  ImageType::Pointer img = ImageType::New();
  ReadImage(img,argv[1]);

  // Read the metis output files
  FILE *f1 = fopen(argv[2],"rt");
  FILE *f2 = fopen(argv[3],"rt");

  // Create output image
  ImageType::Pointer imgOut = ImageType::New();
  imgOut->SetRegions(img->GetBufferedRegion());
  imgOut->Allocate();
  imgOut->FillBuffer(0);

  while(!feof(f2))
    {
    unsigned int label, x,y,z;
    fscanf(f1,"%d\n",&label);
    fscanf(f2,"%d %d %d\n",&x,&y,&z);
    // cout << "label : " << label << " index : " << index << endl;

    itk::Index<3> idx = {{x,y,z}};
    imgOut->SetPixel(idx,(1+label));
    }
  
  fclose(f1);fclose(f2);

  // Write the image
  WriteImage(imgOut,argv[4]);

  return 0;
}
