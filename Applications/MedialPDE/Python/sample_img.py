#! /usr/bin/env python
#
from common2 import *;
import sys
import os

# Usage comand
def usage():
  print """
  Usage: sample_img model_mpde tfield_img output_img ox oy oz nsamples 
    model_mpde        File with the cm-rep along which to sample
    tfield_img        The t-field for that individual
    output_img        Place to store the result
    ox oy oz          The origin of the t-field relative to the origin of the 
                      image from which model_mpde was extracted
    nsamples          Number of samples in the Z direction
  """
  sys.exit(-1)

# Check the number of parameters
if len(sys.argv) <> 8:
  usage()

# Define the filenames
[fnModel, fnImage, fnOutput] = sys.argv[1:4]

# Get the origin
xOrigin = map(float, sys.argv[4:7])

# Get the number of samples
nSamples = int( sys.argv[7] );

# Try loading all the files
for fn in fnModel, fnImage:
  if not os.path.isfile(fn):
    print "Can't access file " + fn
    sys.exit(-1)

# Construct and load the medial PDE
mp = MedialPDE(8, 10, sampling["nu"], sampling["nv"], 
    sampling["cut"], sampling["ncu"], sampling["ncv"])
mp.LoadFromParameterFile(fnModel)

# Load the GLM image
i1 = FloatImage()
i1.LoadFromFile(fnImage)
i1.SetOutsideValue(0.0)
i1.SetOrigin(xOrigin[0], xOrigin[1], xOrigin[2])

# Sample the image
i2 = FloatImage()
mp.SampleImage(i1, i2, nSamples)

# Save the output
i2.SaveToFile(fnOutput)
