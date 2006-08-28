# This is a common include for the cm-rep library

# Import all the symbols from the C++ MPDE library
from medialpde import *

# Import other libraries
import sys
import os
from optparse import OptionParser

# Helper function for loading image gradients
def LoadImageJet(image):
  # Split image filename into name and extension
  (imgstem, imgext) = os.path.splitext(image)

  # Load the image
  jet = FloatImage()
  jet.LoadFromPath(imgstem, imgext[1:])
  jet.SetOutsideValue(-1.0) 
  
  # Return the image
  return jet

