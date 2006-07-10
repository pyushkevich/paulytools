#! /usr/bin/env python
#
from common2 import *;
import sys
import os


# Usage comand
def usage():
  print """
  Usage: save_mesh.py model_mpde out_base_name
    model_mpde        File with the cm-rep that will be saved as a VTK graphic
    out_base_name     The destination name, to which _bnd.vtk and _med.vtk will
                      be appended.
  """
  sys.exit(-1)


# Check the number of parameters
if len(sys.argv) <> 3:
  usage()

# Define the filenames
[fnModel, fnOutput] = sys.argv[1:3]

# Try accessing all the files
for fn in fnModel:
  if not os.path.isfile(fn):
    print "Can't access file " + fn
    sys.exit(-1)

# Construct and load the medial PDE
mp = MedialPDE(8, 10, sampling["nu"], sampling["nv"], 
    sampling["cut"], sampling["ncu"], sampling["ncv"])
mp.LoadFromParameterFile(fnModel)

# Save the VTK mesh
mp.SaveVTKMesh(fnOutput + "_med.vtk", fnOutput + "_bnd.vtk")
