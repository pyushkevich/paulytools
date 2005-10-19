#! /usr/bin/env python
#
from common2 import *;
import sys
import os


# Usage comand
def usage():
  print """
  Usage: vis_sample model_mpde ref_field out_base_name
    model_mpde        File with the cm-rep that will be saved as a VTK graphic
    ref_field         The field in the medial coordinate frame
    out_base_name     The destination name, to which _bnd.vtk and _med.vtk will
                      be appended.
  """
  sys.exit(-1)


# Check the number of parameters
if len(sys.argv) <> 4:
  usage()

# Define the filenames
[fnModel, fnImage, fnOutput] = sys.argv[1:4]

# Try accessing all the files
for fn in fnModel, fnImage:
  if not os.path.isfile(fn):
    print "Can't access file " + fn
    sys.exit(-1)

# Construct and load the medial PDE
mp = MedialPDE(8, 10, sampling["nu"], sampling["nv"], 
    sampling["cut"], sampling["ncu"], sampling["ncv"])
mp.LoadFromParameterFile(fnModel)

# Load the field image 
i1 = FloatImage()
i1.LoadFromFile(fnImage)
i1.SetOutsideValue(0.0)

# Save the VTK mesh
mp.SetIntensityImage(i1)
mp.SaveVTKMesh(fnOutput + "_med.vtk", fnOutput + "_bnd.vtk")
