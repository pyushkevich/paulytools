from medialpde import *;
from common import *;

# Create a medial model
m=MedialPDE(8, 12, 32, 80)

# Load each sample
samples=(
  "1001L", "1002L", "1003L", "1004L", "1005L", "1006L", "1007L", "1008L",
  "1010L", "1011L", "1012L", "1013L", "1014L", "1015L", "1016L", "1017L",
  "1018L", "1019L", "1020L", "1021L", "1022L", "1023L", "1024L", "1025L",
  "1026L", "1027L", "1028L", "1029L", "1030L", "1031L", "1032L", "1033L",
  "1034L", "1035L", "1036L", "1037L", "1038L", "1039L", "1040L", "1041L",
  "1042L", "1043L", "1044L", "1045L", "1046L", "1047L", "1049L", "1050L",
  "1051L", "1052L", "1053L", "1054L", "1055L", "1056L", "1057L", "1058L",
  "1059L", "1060L", "1061L", "1062L", "1063L", "1064L", "1065L", "1066L",
  "1067L", "1068L", "1070L", "1071L", "1072L", "1073L", "1074L", "1075L",
  "1076L", "1077L", "1078L", "1079L", "1080L", "1081L", "1082L", "1083L",
  "1084L", "1085L", "1086L", "1087L", "1088L", "1089L");

# Create the PCA object
pca = MedialPCA()
for file in samples:
  print "Reading " + dirWork + "forpca/st" + file + ".ctf06.mpde"
  m.LoadFromParameterFile(dirWork + "forpca/st" + file + ".ctf06.mpde");
  
  img = FloatImage()
  img.LoadFromFile(dirWork + "t1sample/st" + file + ".hdr");
  m.SetIntensityImage(img)
  
  pca.AddSample(m);

# Compute the PCA proper
pca.ComputePCA();

# Export the PCA matrix in matlab format
pca.ExportShapeMatrix(dirWork + "pca/matrix/shapemat.mat");

# Sample along the first four modes
for i in range(4):
  for j in range(-30,30,1):
    print "Mode %(i)d, Frame %(j)d" % {"i":i,"j":j}      
    pca.SetFSLocationToMean()
    pca.SetFSLocation(i, 0.1 * j)
    pca.GetShapeAtFSLocation(m)

    sub = {'path':dirWork, 'i':i, 'id':j+30}
    fn1 = "%(path)s/pca/pcamed_m%(i)d_%(id)02d.vtk" % sub
    fn2 = "%(path)s/pca/pcabnd_m%(i)d_%(id)02d.vtk" % sub
    m.SaveVTKMesh(fn1, fn2)

    img = FloatImage()
    m.GetIntensityImage(img)
    fn3 = "%(path)s/pca/pcaimg_m%(i)d_%(id)02d.img.gz" % sub
    img.SaveToFile(fn3)
    
    


