#! /usr/bin/env python
#
from medialpde import *;
from common import *;
import sys
import os

# Usage command
def usage():
  print """
  Usage: pcalite expid side
    expid              The experiment ID code (assumes standard filesystem)
    side               The side on which to compute pca (L or R)
  """
  os._exit(-1)

# Check the number of parameters
if(len(sys.argv) <> 3):
  usage()

# Get the experiment id
expid = sys.argv[1]

# Get the side (L or R)
side = sys.argv[2]
if(side != 'L' and side != 'R'):
  usage()

sidelong = "none"
if(side == 'L'):
  sidelong = "left"
else:
  sidelong = "right"

# Create a medial model
m=MedialPDE(8, 10, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);

# Set the sample names
samples=(
  "1001", "1002", "1003", "1004", "1005", "1006", "1007", "1008",
  "1010", "1011", "1012", "1013", "1014", "1015", "1016", "1017",
  "1018", "1019", "1020", "1021", "1022", "1023", "1024", "1025",
  "1026", "1027", "1028", "1029", "1030", "1031", "1032", "1033",
  "1034", "1035", "1036", "1037", "1038", "1039", "1040", "1041",
  "1042", "1043", "1044", "1045", "1046", "1047", "1049", "1050",
  "1051", "1052", "1053", "1054", "1055", "1056", "1057", "1058",
  "1059", "1060", "1061", "1062", "1063", "1064", "1065", "1066",
  "1067", "1068", "1070", "1071", "1072", "1073", "1074", "1075",
  "1076", "1077", "1078", "1079", "1080", "1081", "1082", "1083",
  "1084", "1085", "1086", "1087", "1088", "1089");

# A method to save PCA to a set of files
def SaveMedialPCA(pca, mpde, expid, sidelong, nModes, nSamples):
  """Save the medial PCA to a matrix, mean and mode files"""

  # Create a substitution rule
  sub = {"path":dirWork, "expid":expid, "sidelong":sidelong}

  # Check all directories
  CheckDir("%(path)s/pca/%(expid)s/%(sidelong)s/matrix" % sub)
  CheckDir("%(path)s/pca/%(expid)s/%(sidelong)s/mean" % sub)
  CheckDir("%(path)s/pca/%(expid)s/%(sidelong)s/vtk" % sub)
  CheckDir("%(path)s/pca/%(expid)s/%(sidelong)s/vtk/mean" % sub)
  CheckDir("%(path)s/pca/%(expid)s/%(sidelong)s/vtk/modes" % sub)

  # First, save the shape matrix
  pca.ExportShapeMatrix("%(path)s/pca/%(expid)s/%(sidelong)s/matrix/shapemat.mat" % sub)

  # Next, save the mean shape
  pca.SetFSLocationToMean();
  pca.GetShapeAtFSLocation(mpde);
  mpde.SaveVTKMesh(
      "%(path)s/pca/%(expid)s/%(sidelong)s/vtk/mean/pcamed.mean.vtk" % sub,
      "%(path)s/pca/%(expid)s/%(sidelong)s/vtk/mean/pcabnd.mean.vtk" % sub);
  mpde.SaveToParameterFile("%(path)s/pca/%(expid)s/%(sidelong)s/mean/mean.mpde" % sub);

  # Now save the models along the PCA modes
  for i in range(nModes):
    sub["i"] = i;
    CheckDir("%(path)s/pca/%(expid)s/%(sidelong)s/vtk/modes/mode%(i)d" % sub)

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
          "%(path)s/pca/%(expid)s/%(sidelong)s/vtk/modes/mode%(i)d/pcamed_m%(i)d_%(id)02d.vtk" % sub,
          "%(path)s/pca/%(expid)s/%(sidelong)s/vtk/modes/mode%(i)d/pcabnd_m%(i)d_%(id)02d.vtk" % sub);

      
# Create the PCA object
pca = MedialPCA()

# Add all the samples to the PCA
for file in samples:
  # The filename to read the PDE from
  id = "st" + file + side;
  fnMPDE = dirWork + "cmrep/" + expid + "/" + id + "/" + id + ".bdm92.mpde";

  # Check if the file exists
  if(os.access(fnMPDE, os.R_OK)):
    
    # Load the m-rep if we can
    m1=MedialPDE(9, 12, 
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
SaveMedialPCA(pca, m, expid, sidelong, 4, 24)
