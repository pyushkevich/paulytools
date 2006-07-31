#!/usr/bin/python

# Import the symbols from the MPDE library
from medialpde import *
import sys
import os

# Usage stuff
def usage():
  print """
  Script: fit cmrep template to a 3D binary image
  Usage:  cmrep_binary template.cmrep image.img result.cmrep parameters.reg nSteps
  """
  os._exit(-1)

# Check that the parameters are accurate
if(len(sys.argv) <> 4):
  usage()

# Get the filenames
fnTemplate, fnImage, fnOutput, fnParams = sys.argv[1:5]
nSteps = int(sys[5])

# Load the template from file
cmrep = MedialPDE()
cmrep.LoadFromFile(fnTemplate)

# Load the image and its gradients from file
img = FloatImage()
img.LoadFromPath(fnImage)
img.SetOutsideValue(-1)

# Run the optimization
cmrep.RunOptimization(img, fnParams, nSteps);

# Save the result
cmrep.SaveToFile(fnOutput)
