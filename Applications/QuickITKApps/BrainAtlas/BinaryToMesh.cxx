/**
 * Program: Binary2Mesh
 * Author: Paul Yushkevich
 * Purpose: Convert a binary image to a VTK mesh
 */
#include "BinaryImageToMeshFilter.h"
#include "vtkPolyDataWriter.h"
#include "ReadWriteImage.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

#include "CommandLineArgumentParser.h"

using namespace std;

int usage()
{
  cout << "binary2mesh: extracts iso-surface in binary image" << endl;
  cout << "usage: " << endl;
  cout << "   binary2mesh [options] input.img output.vtk" << endl;
  cout << "options: " << endl;
  cout << "   -r X.XX     Anti-aliasing parameter (default 0.024, lower->better" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Command line options
  double aaParm = 0.024;

  // Check parameters
  if(argc < 2) return usage();
  
  // Read in the command line parameters
  for(unsigned int i=1;i<argc-2;i++)
    {
    if(0 == strcmp(argv[i],"-r")) 
      {
      aaParm = atof(argv[++i]);
      cout << "using mean RMS error of " << aaParm << endl;
      }
    }
  
  // Read the input image
  typedef itk::Image<char, 3> ImageType;
  ImageType::Pointer imgInput;
  ReadImage(imgInput, argv[argc-2]);

  cout << "Image dimensions: " << imgInput->GetBufferedRegion() << endl;

  // Create the mesh
  typedef BinaryImageToMeshFilter<ImageType> FilterType;
  FilterType::Pointer fltMesh = FilterType::New();
  fltMesh->SetInput(imgInput);
  fltMesh->SetAntiAliasMaxRMSError(aaParm);
  fltMesh->Update();

  // Save the mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(fltMesh->GetMesh());
  writer->SetFileName(argv[argc-1]);
  writer->SetFileTypeToBinary();
  writer->Update();
}
