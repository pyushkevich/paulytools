#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGrid.h"
#include "MedialAtom.h"
#include "ITKImageWrapper.h"
#include "itkImage.h"

void ExportMedialMeshToVTK(
  MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, 
  ITKImageWrapper<float> *xImage,
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xGrid->GetNumberOfAtoms());
  
  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xGrid->GetNumberOfAtoms());

  // Allocate the metric tensor array
  vtkFloatArray *lMetric = vtkFloatArray::New();
  lMetric->SetName("Covariant Tensor Determinant");
  lMetric->SetNumberOfComponents(1);
  lMetric->SetNumberOfTuples(xGrid->GetNumberOfAtoms());
  
  // Allocate the metric tensor array
  vtkFloatArray *lRho = vtkFloatArray::New();
  lRho->SetName("Rho Function");
  lRho->SetNumberOfComponents(1);
  lRho->SetNumberOfTuples(xGrid->GetNumberOfAtoms());
  
  // Allocate the metric tensor array
  vtkFloatArray *lRadius = vtkFloatArray::New();
  lRadius->SetName("Radius Function");
  lRadius->SetNumberOfComponents(1);
  lRadius->SetNumberOfTuples(xGrid->GetNumberOfAtoms());

  // Allocate another dummy array
  vtkFloatArray *lDummy1 = vtkFloatArray::New();
  lDummy1->SetName("Dummy1");
  lDummy1->SetNumberOfComponents(1);
  lDummy1->SetNumberOfTuples(xGrid->GetNumberOfAtoms());

  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xGrid->GetNumberOfQuads());
  pMedial->SetPoints(lPoints);
  pMedial->GetPointData()->SetNormals(lNormals);
  pMedial->GetPointData()->AddArray(lMetric);
  pMedial->GetPointData()->AddArray(lRho);
  pMedial->GetPointData()->AddArray(lRadius);
  pMedial->GetPointData()->AddArray(lDummy1);

  // Allocate and add the image intensity array
  bool flagImage = xImage && xImage->IsImageLoaded();
  vtkFloatArray *lImage = vtkFloatArray::New();
  if(flagImage)
    {
    lImage->SetName("Image");
    lImage->SetNumberOfComponents(1);
    lImage->SetNumberOfTuples(xGrid->GetNumberOfAtoms());
    pMedial->GetPointData()->AddArray(lImage);
    }
  
  // Add all the points
  MedialAtomIterator *it = xGrid->NewAtomIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    SMLVec3d X = xAtoms[i].X; SMLVec3d N = xAtoms[i].N;
    lPoints->InsertNextPoint(X[0], X[1], X[2]);
    lNormals->SetTuple3(i, N[0], N[1], N[2]);
    lMetric->SetTuple1(i, xAtoms[i].G.g);
    lRadius->SetTuple1(i, xAtoms[i].R);
    lRho->SetTuple1(i, xAtoms[i].xLapR);

    // Sample the image along the middle
    if(flagImage)
      {
      itk::Index<3> idx;
      idx[0] = xAtoms[i].uIndex;
      idx[1] = xAtoms[i].vIndex;
      idx[2] = xImage->GetImageSize(2) >> 1;
      lImage->SetTuple1(i, xImage->GetInternalImage()->GetPixel(idx));
      }

    int q = (int)(xAtoms[i].u * 64);
    lDummy1->SetTuple1(i, q % 8 == 0 ? 1 : 0);
    }
  delete it;

  // Add all the quads
  MedialQuadIterator *itq = xGrid->NewQuadIterator();
  for(; !itq->IsAtEnd(); ++(*itq))
    {
    vtkIdType xQuad[4];
    xQuad[0] = itq->GetAtomIndex(0, 0);
    xQuad[1] = itq->GetAtomIndex(1, 0);
    xQuad[2] = itq->GetAtomIndex(1, 1);
    xQuad[3] = itq->GetAtomIndex(0, 1);
    pMedial->InsertNextCell(VTK_QUAD, 4, xQuad);
    }
  delete itq;

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lRho->Delete();
  lRadius->Delete();
  lMetric->Delete();
  lNormals->Delete();
  lPoints->Delete();
  lDummy1->Delete();
  lImage->Delete();
  pMedial->Delete();
}

vtkUnstructuredGrid *
ExportVolumeMeshToVTK(MedialAtomGrid *xGrid, MedialAtom *xAtoms, size_t nSamples)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xGrid->GetNumberOfInternalPoints(nSamples));
  
  // Allocate the polydata
  vtkUnstructuredGrid *pMedial = vtkUnstructuredGrid::New();
  pMedial->Allocate(xGrid->GetNumberOfCells(nSamples));
  pMedial->SetPoints(lPoints);

  // Add all the points
  MedialInternalPointIterator *it = xGrid->NewInternalPointIterator(nSamples);
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    SMLVec3d X = GetInternalPoint(it, xAtoms);
    lPoints->InsertPoint(i, X[0], X[1], X[2]);
    }
  delete it;

  // Add all the quads
  MedialInternalCellIterator *itq = xGrid->NewInternalCellIterator(nSamples);
  for(; !itq->IsAtEnd(); ++(*itq))
    {
    vtkIdType xQuad[8];
    xQuad[0] = itq->GetInternalPointIndex(0, 0, 0);
    xQuad[3] = itq->GetInternalPointIndex(1, 0, 0);
    xQuad[1] = itq->GetInternalPointIndex(0, 1, 0);
    xQuad[2] = itq->GetInternalPointIndex(1, 1, 0);
    xQuad[0] = itq->GetInternalPointIndex(0, 0, 1);
    xQuad[3] = itq->GetInternalPointIndex(1, 0, 1);
    xQuad[1] = itq->GetInternalPointIndex(0, 1, 1);
    xQuad[2] = itq->GetInternalPointIndex(1, 1, 1);
    pMedial->InsertNextCell(VTK_HEXAHEDRON, 8, xQuad);
    }
  delete itq;

  /*
  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();
  */

  return pMedial;

}


void ExportBoundaryMeshToVTK(
  MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, 
  ITKImageWrapper<float> *xImage,
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xGrid->GetNumberOfBoundaryPoints());
  
  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xGrid->GetNumberOfBoundaryPoints());
  
  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xGrid->GetNumberOfBoundaryQuads());
  pMedial->SetPoints(lPoints);
  
  // Add the arrays to the poly data
  pMedial->GetPointData()->SetNormals(lNormals);

  // Allocate the image intensity array
  bool flagImage = xImage && xImage->IsImageLoaded();
  vtkFloatArray *lImage = vtkFloatArray::New();
  if(flagImage)
    {
    lImage->SetName("Image");
    lImage->SetNumberOfComponents(1);
    lImage->SetNumberOfTuples(xGrid->GetNumberOfBoundaryPoints());
    pMedial->GetPointData()->AddArray(lImage);
    }

  // Add all the points
  MedialBoundaryPointIterator *it = xGrid->NewBoundaryPointIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    BoundaryAtom &B = GetBoundaryPoint(it, xAtoms);
    lPoints->InsertNextPoint(B.X[0], B.X[1], B.X[2]);
    lNormals->SetTuple3(i, B.N[0], B.N[1], B.N[2]);

    // Sample the image along the middle
    if(flagImage)
      {
      itk::Index<3> idx;
      idx[0] = xAtoms[it->GetAtomIndex()].uIndex;
      idx[1] = xAtoms[it->GetAtomIndex()].vIndex;
      idx[2] = it->GetBoundarySide() ? 0 : xImage->GetImageSize(2) - 1;
      lImage->SetTuple1(i, xImage->GetInternalImage()->GetPixel(idx));
      }

    }
  delete it;

  // Add all the quads
  MedialBoundaryQuadIterator *itq = xGrid->NewBoundaryQuadIterator();
  for(; !itq->IsAtEnd(); ++(*itq))
    {
    vtkIdType xQuad[4];
    xQuad[0] = itq->GetBoundaryIndex(0, 0);
    xQuad[1] = itq->GetBoundaryIndex(0, 1);
    xQuad[2] = itq->GetBoundaryIndex(1, 1);
    xQuad[3] = itq->GetBoundaryIndex(1, 0);
    pMedial->InsertNextCell(VTK_QUAD, 4, xQuad);
    }
  delete itq;

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lNormals->Delete();
  lPoints->Delete();
  pMedial->Delete();
  lImage->Delete();
}

/*
void ExportIntensityFieldToVTK(
  MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, 
  ITKImageWrapper<float> *imgField,
  
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xGrid->GetNumberOfBoundaryPoints());
  
  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xGrid->GetNumberOfBoundaryPoints());

  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xGrid->GetNumberOfBoundaryQuads());
  pMedial->SetPoints(lPoints);
  pMedial->GetPointData()->SetNormals(lNormals);

  // Add all the points
  MedialBoundaryPointIterator *it = xGrid->NewBoundaryPointIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    BoundaryAtom &B = GetBoundaryPoint(it, xAtoms);
    lPoints->InsertNextPoint(B.X[0], B.X[1], B.X[2]);
    lNormals->SetTuple3(i, B.N[0], B.N[1], B.N[2]);
    }
  delete it;

  // Add all the quads
  MedialBoundaryQuadIterator *itq = xGrid->NewBoundaryQuadIterator();
  for(; !itq->IsAtEnd(); ++(*itq))
    {
    vtkIdType xQuad[4];
    xQuad[0] = itq->GetBoundaryIndex(0, 0);
    xQuad[1] = itq->GetBoundaryIndex(0, 1);
    xQuad[2] = itq->GetBoundaryIndex(1, 1);
    xQuad[3] = itq->GetBoundaryIndex(1, 0);
    pMedial->InsertNextCell(VTK_QUAD, 4, xQuad);
    }
  delete itq;

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lNormals->Delete();
  lPoints->Delete();
  pMedial->Delete();
}
*/
