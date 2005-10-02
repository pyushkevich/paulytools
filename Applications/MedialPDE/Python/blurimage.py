# Load all the symbols in the library
from medialpde import *
from common import *
import sys
import os

# Choose one of the hippocampi
id = sys.argv[1]

# Bail out if the data does not exist
if(not os.access(dirWork + "hippo/imgiso/" + id + ".mha", os.R_OK) ):
  print "Can not find appropriate image file!"
  sys.exit(1)

# Create the necessary images
print "Creating blurred images for " + dirWork + "hippo/imgiso/" + id + ".mha"
MakeImages(id)

