#!@PYTHON_BINARY@

# Import the symbols from the MPDE library
from cmrep_common import *

# Usage stuff
def usage():
  print """
  Script: fit cmrep template to a 3D binary image
  Usage:  cmrep_binary template.cmrep image.img result.cmrep parameters.reg nSteps
  """
  os._exit(-1)

# Check that the parameters are accurate
if(len(sys.argv) <> 6):
  usage()

# Get the filenames
fnTemplate, fnImage, fnOutput, fnParams = sys.argv[1:5]
nSteps = int(sys.argv[5])

# Load the template from file
cmrep = MedialPDE(fnTemplate)

# Load the image and its gradients from file
img = LoadImageJet(fnImage)

# Run the optimization
cmrep.RunOptimization(img, nSteps, fnParams);

# Save the result
cmrep.SaveToParameterFile(fnOutput)
