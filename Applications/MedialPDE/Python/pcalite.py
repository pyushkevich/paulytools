from medialpde import *;
from common import *;
import sys
import os

# Get the experiment id
expid = sys.argv[1]

# Create a medial model
m=MedialPDE(8, 10, 32, 80)

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

# Add all the samples to the PCA
for file in samples:
  # The filename to read the PDE from
  id = "st" + file;
  fnMPDE = dirWork + "cmrep/" + expid + "/" + id + "/" + id + ".ctf80.mpde";

  # Check if the file exists
  if(os.access(fnMPDE, os.R_OK)):
    # Load the m-rep
    m.LoadFromParameterFile(fnMPDE)
    
    # Add the sample to the PCA  
    print "Added " + id + " to the PCA"
    pca.AddSample(m);

# Compute the PCA proper
pca.ComputePCA();

# Export the PCA matrix in matlab format
fnMatFile = CheckDir(dirWork + "pca/" + expid + "/matrix") + "/shapemat.mat"
pca.ExportShapeMatrix(fnMatFile);

# Export the mean model
pca.SetFSLocationToMean()
pca.GetShapeAtFSLocation(m)
dirVTK = CheckDir(dirWork + "/pca/" + expid +"/vtk")
dirMean = CheckDir(dirWork + "/pca/" + expid +"/mean")
m.SaveVTKMesh(dirVTK + "/pcamed_mean.vtk", dirVTK + "/pcabnd_mean.vtk")
m.SaveToParameterFile(dirMean + "/mean.mpde");
