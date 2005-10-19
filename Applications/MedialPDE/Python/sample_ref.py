#! /usr/bin/env python
#
from common2 import *;
import sys
import os


# Usage comand
def usage():
  print """
  Usage: sample_ref model_mpde ref_field mask_img output_img
    model_mpde        File with the cm-rep along which to sample
    ref_field         The field in the medial coordinate frame
    mask_img          A mask image to which the output will be mapped
    output_img        The destination image
    nsamples          Number of samples in the Z direction
  """
  sys.exit(-1)


# Check the number of parameters
if len(sys.argv) <> 6:
  usage()

# Define the filenames
[fnModel, fnImage, fnMask, fnOutput] = sys.argv[1:5]

# Get the number of samples
nSamples = int( sys.argv[5] );

# Try accessing all the files
for fn in fnModel, fnImage, fnMask:
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

# Load the mask image 
i2 = FloatImage()
i2.LoadFromFile(fnMask)
i2.SetOutsideValue(0.0)

# Perform the sampling
mp.SampleReferenceFrameImage(i1, i2, nSamples)

# Save the result of the sampling
i2.SaveToFile(fnOutput)


