#include "ReadWriteVTK.h"
#include <vtkBYUReader.h>
#include <vtkBYUWriter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

using namespace std;

vtkPolyData *ReadVTKData(string fn)
{
  vtkPolyData *p1 = NULL;

  // Choose the reader based on extension

  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".stl") == fn.length() - 4)
    {
    vtkSTLReader *reader = vtkSTLReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else
    {
    cout << "Could not find a reader for " << fn << endl;
    return NULL;
    }
}

void WriteVTKData(vtkPolyData *data, string fn)
{
  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUWriter *writer = vtkBYUWriter::New();
    writer->SetFileName(fn.c_str());
    writer->SetInput(data);
    writer->Update();
    }
  else if(fn.rfind(".stl") == fn.length() - 4)
    {
    vtkSTLWriter *writer = vtkSTLWriter::New();
    writer->SetFileName(fn.c_str());
    writer->SetInput(data);
    writer->Update();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(fn.c_str());
    writer->SetInput(data);
    writer->Update();
    }
  else
    {
    cout << "Could not find a writer for " << fn << endl;
    cout << fn.rfind(".byu") << endl;
    cout << fn.rfind(".vtk") << endl;
    return;
    }
}
