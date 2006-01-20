#! /usr/bin/env python
#
from pipeline2 import *;

# Usage command
def usage():
  print """
  Usage: segment_new input_image input_template output_dir
    input_image             The location of the binary image volume to fit
    input_template          The location of the template MPDE
    output_dir              The directory where the output will go to 
  """
  os._exit(-1)

# Check the number of parameters
if(len(sys.argv) <> 4):
  usage()

# Check that the input image exists
for i in range(1,3):
  if(not os.path.isfile(sys.argv[i])):
    print "Missing file " , sys.argv[i]
    usage()

# Define the output paths
dirOutput = sys.argv[3]
expdata = {
  "id"      : "segment",
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
Stage_MOInertia(expdata, sys.argv[2], "align", "med", 9, 12)

# Affine stage
Stage_AFF_CG_VO(expdata, "align",  "affine", "med", 9, 12, 600);

# Deformation stages
Stage_CTF_CG_VO(expdata, "affine",  "ctf24",  "med", 9, 12, 2, 4, 1600);
Stage_CTF_CG_VO(expdata, "ctf24",  "ctf35",  "med", 9, 12, 3, 5, 600);
Stage_CTF_CG_VO(expdata, "ctf35",  "ctf46",  "med", 9, 12, 4, 6, 600);
Stage_CTF_CG_VO(expdata, "ctf46",  "ctf57",  "med", 9, 12, 5, 7, 600);
Stage_CTF_CG_VO(expdata, "ctf57",  "ctf68",  "med", 9, 12, 6, 8, 600);
Stage_CTF_CG_VO(expdata, "ctf68",  "ctf79",  "med", 9, 12, 7, 9, 600);
Stage_CTF_CG_VO(expdata, "ctf79",  "ctf80",  "med", 9, 12, 8, 10, 600);

# Boundary image match stage
Stage_CTF_CG_BM(expdata, "ctf80",  "bdm80",  "low", 9, 12, 8, 10, 600);
Stage_CTF_CG_BM(expdata, "bdm80",  "bdm92",  "low", 9, 12, 9, 12, 600);
