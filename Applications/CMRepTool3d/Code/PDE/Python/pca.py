from medialpde import *;
from common import *;

# Create a medial model
m=MedialPDE(6,10,60,120)

# Load the mean
m.LoadFromParameterFile(dirWork + "pca/mean/meanleft.mpde");
m.SaveVTKMesh(dirWork+"pca/mean/meanleft.med.vtk",dirWork+"pca/mean/meanleft.bnd.vtk");

# Load the first PCA mode
for i in range(1,81):
  id=str(1000+i)[-3:];
  print "Processing: ", id
  m.LoadFromParameterFile(dirWork + "pca/left_mode1/pca" + id + ".mpde");
  m.SaveVTKMesh(
    dirWork + "pca/left_mode1/vtk/pca" + id + ".med.vtk", 
    dirWork + "pca/left_mode1/vtk/pca" + id + ".bnd.vtk");
