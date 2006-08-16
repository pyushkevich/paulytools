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
  cout << "   -a X.XX     Anti-aliasing parameter (default 0.024, lower->better" << endl;
  cout << "   -s X.XX     Scale the image by factor X.XX before mesh extraction" << endl;
  cout << "   -ea FILE    Export the anti-aliased image as an image file FILE" << endl;
  cout << "   -es FILE    Export the scaled image as an image file FILE" << endl;
  cout << "   -d X.XX     Decimate the mesh by a factor of X.XX (between 0 and 1)" << endl;
  cout << "   -sm NNN      Smooth the mesh for NNN iterations" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Command line options
  double aaParm = 0.024;
  double xScale = 1.0;
  double xDecimate = 0.0;
  int nSmoothIterations = 0;
  const char *fnOutAA = NULL, *fnOutResample = NULL;

  // Check parameters
  if(argc < 2) return usage();
  
  // Read in the command line parameters
  for(unsigned int i=1;i<argc-2;i++)
    {
    if(0 == strcmp(argv[i],"-a")) 
      {
      aaParm = atof(argv[++i]);
      cout << "using mean RMS error of " << aaParm << endl;
      }
    else if(0 == strcmp(argv[i],"-d")) 
      {
      xDecimate = atof(argv[++i]);
      cout << "will decimate besh by factor of " << xDecimate << endl;
      }
    else if(0 == strcmp(argv[i],"-ea"))
      {
      fnOutAA = argv[++i];
      cout << "will export anti-alias image to " << fnOutAA << endl;
      }
    else if(0 == strcmp(argv[i], "-s"))
      {
      xScale = atof(argv[++i]);
      cout << "scaling the image by " << xScale << ", gonna be slow!" << endl;
      }
    else if(0 == strcmp(argv[i],"-es"))
      {
      fnOutResample = argv[++i];
      cout << "will export resampled image to " << fnOutResample << endl;
      }
    else if(0 == strcmp(argv[i],"-sm"))
      {
      nSmoothIterations = atoi(argv[++i]);
      cout << "will smooth the mesh using " << nSmoothIterations << "iterations"  << endl;
      }
    else
      {
      cerr << "Bad option " << argv[i] << endl;
      return usage();
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
  fltMesh->SetResampleScaleFactor(xScale);
  fltMesh->SetAntiAliasMaxRMSError(aaParm);
  fltMesh->SetDecimateFactor(xDecimate);
  fltMesh->SetSmoothingIterations(nSmoothIterations);

  cout << "==========================================================" << endl;
  fltMesh->Update();
  cout << "==========================================================" << endl;
  
  // If needed save the resampled image
  if(fnOutResample)
    {
    FilterType::FloatImageType::Pointer imgResample = fltMesh->GetResampledImage();
    WriteImage(imgResample, fnOutResample);
    }

  // If needed, save the anti-alias image
  if(fnOutAA)
    {
    FilterType::FloatImageType::Pointer imgAlias = fltMesh->GetAntiAliasImage();
    WriteImage(imgAlias, fnOutAA);
    }

  // Save the mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(fltMesh->GetMesh());
  writer->SetFileName(argv[argc-1]);
  writer->SetFileTypeToBinary();
  writer->Update();
}
