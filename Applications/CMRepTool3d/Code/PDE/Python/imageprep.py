# Load all the symbols in the library
from common import *
from medialpde import *

# Load a binary image
bin = BinaryImage()
bin.LoadFromFile(dirWork + "avg/average_hippoL.mha")

# Convert into the floating point image with gradient
img = FloatImage()

# Apply different levels of blurring to the image
img.SetToBlurredBinary(bin, 2.0)
img.SetOutsideValue(-1.0)
img.SaveToPath(dirWork + "avg/average_hippo_blurred_hi","mha")

# Apply different levels of blurring to the image
img.SetToBlurredBinary(bin, 1.0)
img.SetOutsideValue(-1.0)
img.SaveToPath(dirWork + "avg/average_hippo_blurred_med","mha")

# Apply different levels of blurring to the image
img.SetToBlurredBinary(bin, 0.6)
img.SetOutsideValue(-1.0)
img.SaveToPath(dirWork + "avg/average_hippo_blurred_low","mha")
