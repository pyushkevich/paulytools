#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGrid.h"
#include "MedialAtom.h"
#include "MedialAtomGrid.h"
#include "GenericMedialModel.h"
#include "ITKImageWrapper.h"
#include "itkImage.h"

vtkFloatArray *AddMedialScalarField(vtkPolyData *target, GenericMedialModel *model, char *name)
{
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetName(name);
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(model->GetNumberOfAtoms());
  target->GetPointData()->AddArray(array);
  return array;
}

vtkFloatArray *AddMedialVectorField(vtkPolyData *target, GenericMedialModel *model, char *name)
{
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetName(name);
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(model->GetNumberOfAtoms());
  target->GetPointData()->AddArray(array);
  return array;
}

void ExportMedialMeshToVTK(
  GenericMedialModel *xModel,
  ITKImageWrapper<float> *xImage,
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xModel->GetNumberOfAtoms());
  
  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xModel->GetNumberOfTriangles());
  pMedial->SetPoints(lPoints);

  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xModel->GetNumberOfAtoms());
  pMedial->GetPointData()->SetNormals(lNormals);

  // Allocate the scalar arrays
  vtkFloatArray *lMetric = AddMedialScalarField(pMedial, xModel, "Covariant Tensor Determinant");
  vtkFloatArray *lRho = AddMedialScalarField(pMedial, xModel, "Rho Function");
  vtkFloatArray *lRadius = AddMedialScalarField(pMedial, xModel, "Radius Function");
  vtkFloatArray *lDummy1 = AddMedialScalarField(pMedial, xModel, "Dummy1");
  vtkFloatArray *lBending = AddMedialScalarField(pMedial, xModel, "Bending Energy");
  vtkFloatArray *lRegularity = AddMedialScalarField(pMedial, xModel, "Regularity Penalty");
  vtkFloatArray *lAngle = AddMedialScalarField(pMedial, xModel, "Metric Angle");
  vtkFloatArray *lCoordU = AddMedialScalarField(pMedial, xModel, "U Coordinate");
  vtkFloatArray *lCoordV = AddMedialScalarField(pMedial, xModel, "V Coordinate");
  vtkFloatArray *lMeanCurv = AddMedialScalarField(pMedial, xModel, "Mean Curvature");
  vtkFloatArray *lGaussCurv = AddMedialScalarField(pMedial, xModel, "Gauss Curvature");
  vtkFloatArray *lKappa1 = AddMedialScalarField(pMedial, xModel, "Kappa1");
  vtkFloatArray *lKappa2 = AddMedialScalarField(pMedial, xModel, "Kappa2");
  vtkFloatArray *lNormal = AddMedialVectorField(pMedial, xModel, "Atom Normal");
  vtkFloatArray *lStretch = AddMedialScalarField(pMedial, xModel, "Stretch Norm");
  vtkFloatArray *lCurvPen = AddMedialScalarField(pMedial, xModel, "Curvature Penalty Feature");

  vtkFloatArray *lContraOffDiag = 
    AddMedialScalarField(pMedial, xModel, "Off Diagonal Term of Contravariant MT");

  // Allocate and add the image intensity array
  bool flagImage = xImage && xImage->IsImageLoaded();
  vtkFloatArray *lImage = vtkFloatArray::New();
  if(flagImage)
    {
    lImage->SetName("Image");
    lImage->SetNumberOfComponents(1);
    lImage->SetNumberOfTuples(xModel->GetNumberOfAtoms());
    pMedial->GetPointData()->AddArray(lImage);
    }
  
  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();

  // Add all the points
  for(MedialAtomIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    MedialAtom &a = xAtoms[i];

    SMLVec3d X = a.X; SMLVec3d N = a.N;
    lPoints->InsertNextPoint(X[0], X[1], X[2]);
    lNormals->SetTuple3(i, N[0], N[1], N[2]);
    lMetric->SetTuple1(i, a.G.g);
    lRadius->SetTuple1(i, a.R);
    lRho->SetTuple1(i, a.xLapR);

    lCoordU->SetTuple1(i, a.u);
    lCoordV->SetTuple1(i, a.v);

    lMeanCurv->SetTuple1(i, a.xMeanCurv);
    lGaussCurv->SetTuple1(i, a.xGaussCurv);
    lKappa1->SetTuple1(i, a.xMeanCurv - sqrt( a.xMeanCurv *  a.xMeanCurv -  a.xGaussCurv));
    lKappa2->SetTuple1(i, a.xMeanCurv + sqrt( a.xMeanCurv *  a.xMeanCurv -  a.xGaussCurv));

    lNormal->SetTuple3(i, a.N(0), a.N(1), a.N(2));

    // Compute the stretch ???
    vnl_matrix<double> J(3,2);
    J.set_column(0, a.Xu);
    J.set_column(1, a.Xv);
    vnl_svd<double> svd(J);
    double v1 = svd.W()(0,0);
    double v2 = svd.W()(1,1);
    lStretch->SetTuple1(i, sqrt(v1*v1 + v2*v2));

    // Compute sum of the squares of principal curvatures
    double k2 = 4 * a.xMeanCurv * a.xMeanCurv - 2 * a.xGaussCurv;
    lCurvPen->SetTuple1(i, k2);

    // Set the bending energy
    lBending->SetTuple1(i,
      dot_product(a.Xuu, a.Xuu) + dot_product(a.Xvv,a.Xvv) + 2.0 * dot_product(a.Xuv, a.Xuv));

    // Set the regularity energy
    double reg1 = a.G.xChristoffelSecond[0][0][0] + a.G.xChristoffelSecond[1][0][1];
    double reg2 = a.G.xChristoffelSecond[0][1][0] + a.G.xChristoffelSecond[1][1][1];
    double reg = reg1 * reg1 + reg2 * reg2;
    lRegularity->SetTuple1(i, reg);

    // Set the angle between Xu and Xv
    double dp = dot_product(a.Xu, a.Xv);

    lAngle->SetTuple1(i,
      (dp * dp) / (a.Xu.squared_magnitude() * a.Xv.squared_magnitude()));
    lContraOffDiag->SetTuple1(i, 
      a.G.xCovariantTensor[0][1] * a.G.xCovariantTensor[0][1] / a.G.g);

    // Sample the image along the middle
    if(flagImage)
      {
      itk::Index<3> idx;
      idx[0] = xAtoms[i].uIndex;
      idx[1] = xAtoms[i].vIndex;
      idx[2] = xImage->GetImageSize(2) >> 1;
      lImage->SetTuple1(i, xImage->GetInternalImage()->GetPixel(idx));
      }

    double du = xAtoms[i].u - 0.25 * floor(4 * xAtoms[i].u + 0.5);
    double dv = xAtoms[i].v - 0.25 * floor(4 * xAtoms[i].v + 0.5);
    double del = std::min(du * du, dv * dv);
    double q = exp(-0.01 * del * del);
    lDummy1->SetTuple1(i, q);
    }

  // Add all the quads
  for(MedialTriangleIterator itt(xGrid); !itt.IsAtEnd(); ++itt)
    {
    vtkIdType xTri[4];
    xTri[0] = itt.GetAtomIndex(0);
    xTri[1] = itt.GetAtomIndex(1);
    xTri[2] = itt.GetAtomIndex(2);
    pMedial->InsertNextCell(VTK_TRIANGLE, 3, xTri);
    }

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lCoordU->Delete();
  lCoordV->Delete();
  lRho->Delete();
  lRadius->Delete();
  lMetric->Delete();
  lNormals->Delete();
  lPoints->Delete();
  lDummy1->Delete();
  lImage->Delete();
  lBending->Delete();
  lRegularity->Delete();
  lContraOffDiag->Delete();
  lAngle->Delete();
  lMeanCurv->Delete();
  lGaussCurv->Delete();
  pMedial->Delete();
}

vtkUnstructuredGrid *
ExportVolumeMeshToVTK(GenericMedialModel *xModel, size_t nSamples)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xModel->GetNumberOfInternalPoints(nSamples));
  
  // Allocate the polydata
  vtkUnstructuredGrid *pMedial = vtkUnstructuredGrid::New();
  pMedial->Allocate(xModel->GetNumberOfCells(nSamples));
  pMedial->SetPoints(lPoints);

  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();

  // Add all the points
  for(MedialInternalPointIterator it(xGrid,nSamples); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    SMLVec3d X = GetInternalPoint(it, xAtoms);
    lPoints->InsertPoint(i, X[0], X[1], X[2]);
    }

  // Add all the quads
  for(MedialInternalCellIterator itc(xGrid,nSamples); !itc.IsAtEnd(); ++itc)
    {
    vtkIdType xWedge[6];
    xWedge[0] = itc.GetInternalPointIndex(0, 0);
    xWedge[1] = itc.GetInternalPointIndex(1, 0);
    xWedge[2] = itc.GetInternalPointIndex(2, 0);
    xWedge[3] = itc.GetInternalPointIndex(0, 1);
    xWedge[4] = itc.GetInternalPointIndex(1, 1);
    xWedge[5] = itc.GetInternalPointIndex(2, 1);
    pMedial->InsertNextCell(VTK_WEDGE, 6, xWedge);
    }

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
  GenericMedialModel *xModel,
  ITKImageWrapper<float> *xImage,
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xModel->GetNumberOfBoundaryPoints());
  
  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xModel->GetNumberOfBoundaryPoints());
  
  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xModel->GetNumberOfBoundaryTriangles());
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
    lImage->SetNumberOfTuples(xModel->GetNumberOfBoundaryPoints());
    pMedial->GetPointData()->AddArray(lImage);
    }

  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();

  // Add all the points
  for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    BoundaryAtom &B = GetBoundaryPoint(it, xAtoms);
    lPoints->InsertNextPoint(B.X[0], B.X[1], B.X[2]);
    lNormals->SetTuple3(i, B.N[0], B.N[1], B.N[2]);

    // Sample the image along the middle
    if(flagImage)
      {
      itk::Index<3> idx;
      idx[0] = xAtoms[it.GetAtomIndex()].uIndex;
      idx[1] = xAtoms[it.GetAtomIndex()].vIndex;
      idx[2] = it.GetBoundarySide() ? 0 : xImage->GetImageSize(2) - 1;
      lImage->SetTuple1(i, xImage->GetInternalImage()->GetPixel(idx));
      }

    }

  // Add all the quads
  for(MedialBoundaryTriangleIterator itt(xGrid); !itt.IsAtEnd(); ++itt)
    {
    vtkIdType xTri[3];
    xTri[0] = itt.GetBoundaryIndex(0);
    xTri[1] = itt.GetBoundaryIndex(1);
    xTri[2] = itt.GetBoundaryIndex(2);
    pMedial->InsertNextCell(VTK_TRIANGLE, 3, xTri);
    }

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
