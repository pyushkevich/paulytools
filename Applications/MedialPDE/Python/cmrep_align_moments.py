#!/usr/bin/python
from cmrep_common import *

# Get the command line arguments
(template, image, output) = sys.argv[1:4]

# Load the image from file
img = LoadImageJet(image)

# Load the template
mrep = MedialPDE(template)

# Peform the moments alignment
mrep.MatchImageByMoments(img, 5)

# Save the model
mrep.SaveToParameterFile(output)
