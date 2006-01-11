#! /usr/bin/env python
#
from pipeline2 import *;

# Usage command
def usage():
  print """
  Usage: pca_segment input_image pca_path output_dir
    pca_path                The location of the shape statistics directory
    input_image             The location of the binary image volume to fit
    output_dir              The directory where the output will go to 
  """
  os._exit(-1)

# Check the number of parameters
if(len(sys.argv) <> 4):
  usage()

# Define the program parameters
pcaData = {
  "mean":   os.path.join(sys.argv[2], "mean/mean.mpde"),
  "matrix": os.path.join(sys.argv[2], "matrix/shapemat.mat"),
  "ncu":8, 
  "ncv":10
}

# Check that the input image exists
if(not os.path.isfile(sys.argv[1])):
  usage()

# Check that the PCA data is correctly specified
if(not os.path.isfile(pcaData["mean"]) or not os.path.isfile(pcaData["matrix"])):
  print "Missing PCA data (" + pcaData["mean"] + " or " + pcaData["matrix"] + ")"
  usage()

# Define the output paths
dirOutput = sys.argv[3]
expdata = {
  "id"      : "pcafit",
  "root"    : dirOutput,
  "blur"    : os.path.join(dirOutput,"imgblur"),
  "cmrep"   : os.path.join(dirOutput,"cmrep"),
  "vtk"     : os.path.join(dirOutput,"vtk"),
  "dump"    : os.path.join(dirOutput,"meshdump")
}

# Create the output directory and subdirectories
for i in "root", "blur", "cmrep", "vtk", "dump":
  CheckDir(expdata[i])

# Crate blur image if necessary
MakeImages(expdata, sys.argv[1], False)

# Align template by moments of inertia
Stage_MOInertia(expdata, pcaData["mean"], "align","med", 8, 10)

# Affine stage
Stage_AFF_CG_VO(expdata, "align",  "affine", "med", 8, 10, 600);

# PCA Stage
Stage_PCA_CG_VO(expdata, "affine", "pca04",  "med", 8, 10, 600, pcaData, 4);

# Second PCA Stage
Stage_PCA_CG_VO(expdata, "pca04",  "pca12",  "med", 8, 10, 600, pcaData, 12);

# Deformation stages
Stage_CTF_CG_VO(expdata, "pca12",  "ctf24",  "med", 8, 10, 2, 4, 1600);
Stage_CTF_CG_VO(expdata, "ctf24",  "ctf35",  "med", 8, 10, 3, 5, 600);
Stage_CTF_CG_VO(expdata, "ctf35",  "ctf46",  "med", 8, 10, 4, 6, 600);
Stage_CTF_CG_VO(expdata, "ctf46",  "ctf57",  "med", 8, 10, 5, 7, 600);
Stage_CTF_CG_VO(expdata, "ctf57",  "ctf68",  "med", 8, 10, 6, 8, 600);
Stage_CTF_CG_VO(expdata, "ctf68",  "ctf79",  "med", 8, 10, 7, 9, 600);
Stage_CTF_CG_VO(expdata, "ctf79",  "ctf80",  "med", 8, 10, 8, 10, 600);

# Boundary image match stage
Stage_CTF_CG_BM(expdata, "ctf80",  "bdm80",  "low", 8, 10, 8, 10, 600);
Stage_CTF_CG_BM(expdata, "bdm80",  "bdm92",  "low", 9, 12, 9, 12, 600);
