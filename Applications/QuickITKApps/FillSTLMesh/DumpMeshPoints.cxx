#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkCommand.h>
#include <itkImage.h>

#include "DrawTriangles.h"
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "dumpmeshpoints - Dumps points in a VTK mesh, q-hull compatible" << endl;
  cout << "usage: " << endl;
  cout << "   dumpmeshpoints [options] <input.byu|input.vtk> " << endl;
  cout << "options: " << endl;
  cout << "   -a NAME   Dump the contents of array name as a column" << endl;
  cout << "   -e        Dump the edges in the mesh" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 2) return usage();

  // Read the appropriate mesh
  vtkPolyData *p1 = ReadVTKData(argv[argc-1]);

  // Get the list of arrays that need to be dumped along with points
  std::vector<vtkDataArray *> vdat;
  for(size_t i = 1; i < argc - 1; ++i)
    {
    if(0 == strcmp(argv[i], "-a"))
      {
      vtkDataArray *arr = p1->GetPointData()->GetArray(argv[++i]);
      if(arr)
        vdat.push_back(arr);
      else
        cerr << "No array by name " << argv[i] << " found in the mesh" << endl;
      }
    else if(0 == strcmp(argv[i], "-e"))
      {
      // Get the edge filename
      string fnEdge = argv[++i];

      // Get all the edges in the mesh
      typedef std::pair<int, int> Edge;
      std::set<Edge> eset;

      // Go through the cells and count all edges
      p1->BuildCells();
      p1->BuildLinks();
      for(vtkIdType id = 0; id < p1->GetNumberOfCells(); id++)
        {
        vtkIdList *nbr = vtkIdList::New();
        p1->GetCellPoints(id, nbr);
        for(size_t j = 0; j < nbr->GetNumberOfIds(); j++)
          {
          vtkIdType q1 = nbr->GetId(j);
          vtkIdType q2 = nbr->GetId((j + 1) % nbr->GetNumberOfIds());
          eset.insert(make_pair(std::min(q1, q2), std::max(q1,q2)));
          }
        nbr->Delete();
        }

      // Save all the edges
      cout << eset.size() << endl;
      for(std::set<Edge>::const_iterator it = eset.begin(); it!=eset.end(); ++it)
        {
        cout << it->first << " " << it->second << endl;
        }
      }

    }

  // Get the points from the mesh
  vtkPoints *pts = p1->GetPoints();

  // Scale all the points in the mesh
  cout << pts->GetNumberOfPoints() << endl;
  cout << 3 + vdat.size() << endl;
  for(size_t iPoint = 0; iPoint < pts->GetNumberOfPoints(); iPoint++)
    {
    cout << pts->GetPoint(iPoint)[0] << " " 
      << pts->GetPoint(iPoint)[1] << " " 
      << pts->GetPoint(iPoint)[2];
    for(size_t k = 0; k < vdat.size(); k++)
      cout << " " << vdat[k]->GetTuple1(iPoint);
    cout << endl;
    }
}
