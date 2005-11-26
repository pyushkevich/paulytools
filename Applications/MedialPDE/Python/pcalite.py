#! /usr/bin/env python
#
from medialpde import *;
from common import *;
import sys
import os

# Usage command
def usage():
  print """
  Usage: pcalite expid
    expid              The experiment ID code (assumes standard filesystem)
  """
  os._exit(-1)

# Check the number of parameters
if(len(sys.argv) <> 2):
  usage()

# Get the experiment id
expid = sys.argv[1]

# Create a medial model
m=MedialPDE(8, 10, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);

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

# A method to save PCA to a set of files
def SaveMedialPCA(pca, mpde, expid, nModes, nSamples):
  """Save the medial PCA to a matrix, mean and mode files"""

  # Create a substitution rule
  sub = {"path":dirWork, "expid":expid}

  # Check all directories
  CheckDir("%(path)s/pca/%(expid)s/matrix" % sub)
  CheckDir("%(path)s/pca/%(expid)s/mean" % sub)
  CheckDir("%(path)s/pca/%(expid)s/vtk" % sub)
  CheckDir("%(path)s/pca/%(expid)s/vtk/mean" % sub)
  CheckDir("%(path)s/pca/%(expid)s/vtk/modes" % sub)

  # First, save the shape matrix
  pca.ExportShapeMatrix("%(path)s/pca/%(expid)s/matrix/shapemat.mat" % sub)

  # Next, save the mean shape
  pca.SetFSLocationToMean();
  pca.GetShapeAtFSLocation(mpde);
  mpde.SaveVTKMesh(
      "%(path)s/pca/%(expid)s/vtk/mean/pcamed.mean.vtk" % sub,
      "%(path)s/pca/%(expid)s/vtk/mean/pcabnd.mean.vtk" % sub);
  mpde.SaveToParameterFile("%(path)s/pca/%(expid)s/mean/mean.mpde" % sub);

  # Now save the models along the PCA modes
  for i in range(nModes):
    sub["i"] = i;
    CheckDir("%(path)s/pca/%(expid)s/vtk/modes/mode%(i)d" % sub)

    # Iterate over the samples in this mode
    for j in range(-nSamples/2, nSamples/2 + 1, 1):
      sub["j"] = j; sub["id"] = j + nSamples/2; sub["z"] = j * 6.0 / nSamples;
      print ("Mode %(i)d, Frame %(j)d, z = %(z)g" % sub);

      # Go to the PCA location
      pca.SetFSLocationToMean()
      pca.SetFSLocation(i, sub["z"]);
      pca.GetShapeAtFSLocation(mpde)

      # Save the VTK mesh
      mpde.SaveVTKMesh(
          "%(path)s/pca/%(expid)s/vtk/modes/mode%(i)d/pcamed_m%(i)d_%(id)02d.vtk" % sub,
          "%(path)s/pca/%(expid)s/vtk/modes/mode%(i)d/pcabnd_m%(i)d_%(id)02d.vtk" % sub);

      
# Create the PCA object
pca = MedialPCA()

# Add all the samples to the PCA
for file in samples:
  # The filename to read the PDE from
  id = "st" + file;
  fnMPDE = dirWork + "cmrep/" + expid + "/" + id + "/" + id + ".ctf80.mpde";

  # Check if the file exists
  if(os.access(fnMPDE, os.R_OK)):
    
    # Load the m-rep if we can
    m1=MedialPDE(8, 10, 
        sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);

    if(m1.LoadFromParameterFile(fnMPDE)):

      # Check the Jacobian penalty
      m1.ComputeBoundaryJacobianPenalty(True)
      
      # Add the sample to the PCA  
      print "Added " + id + " to the PCA"
      pca.AddSample(m1);

# Compute the PCA proper
print "Computing PCA!"
pca.ComputePCA();

# Export the PCA matrix in matlab format
SaveMedialPCA(pca, m, expid, 4, 24)
