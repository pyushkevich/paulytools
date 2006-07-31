#ifndef __MedialModelIO_h_
#define __MedialModelIO_h_

#include "Registry.h"
#include "CartesianMedialModel.h"
#include "SubdivisionMedialModel.h"

class vtkPolyData;

class CartesianMedialModelIO {
public:
  static CartesianMedialModel *ReadModel(Registry &R);
  static void WriteModel(CartesianMedialModel *model, const char *file);
};

class SubdivisionMedialModelIO {
public:
  static SubdivisionMedialModel *ReadModel(Registry &R);
  static void WriteModel(SubdivisionMedialModel *model, const char *file);

private:
  static vtkPolyData *ReadMesh(const std::string &name, const std::string &type);
};

class MedialModelIO {
public:
  static GenericMedialModel *ReadModel(const char *file);
  static void WriteModel(GenericMedialModel *model, const char *file);
};

#endif // __MedialModelIO_h_
