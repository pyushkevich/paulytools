# Load all the symbols in the library
from common import *
from medialpde import *

# Load a binary image
bin = BinaryImage()
bin.LoadFromFile(dirWork + "avg/average_hippoL.mha")

# Convert into the floating point image with gradient
img = FloatImage()
img.SetToBlurredBinary(bin, 1.0)
img.SetOutsideValue(-1.0)

# Save the image and the gradient
img.SaveToPath(dirWork + "avg/average_hippo_blurred","mha")
