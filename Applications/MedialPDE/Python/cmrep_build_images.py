#!@PYTHON_BINARY@
from cmrep_common import *

# Get the command line arguments
(inimg, outimg, blur) = sys.argv[1:4]

# Split the output into extension and stem
(outstem, outext) = os.path.splitext(outimg)

# Create the image
bin = BinaryImage()
bin.LoadFromFile(inimg)

# Convert image to floating point and blur
img = FloatImage()
img.SetToBlurredBinary(bin, float(blur))
img.SetOutsideValue(-1.0)
img.SaveToPath(outstem, outext[1:])

