from medialpde import *;
from common import *;

# Create a medial model
m=MedialPDE(6,10,60,120)

# Load each sample
samples=(
  "1000L", "1001L", "1002L", "1003L", "1004L", "1006L", "1007L", 
  "1008L", "1010L", "1011L", "1012L", "1013L", "1015L", "1016L", 
  "1017L", "1018L", "1019L", "1021L", "1022L", "1023L", "1024L", 
  "1025L", "1026L", "1027L", "1028L", "1029L", "1030L", "1031L", 
  "1032L", "1033L", "1034L", "1035L", "1036L", "1039L", "1040L", 
  "1041L", "1042L", "1043L", "1044L", "1045L", "1046L", "1049L", 
  "1050L", "1051L", "1052L", "1053L", "1054L", "1055L", "1056L", 
  "1057L", "1058L");

# Create the PCA object
pca = MedialPCA()
for file in samples:
  print "Reading " + dirWork + "forpca/st" + file + ".ctf06.mpde"
  m.LoadFromParameterFile(dirWork + "forpca/st" + file + ".ctf06.mpde");
  pca.AddSample(m);

# Compute the PCA proper
pca.ComputePCA();
