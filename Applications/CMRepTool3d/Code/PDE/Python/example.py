# Load all the symbols in the library
from medialpde import *

# Create a medial PDE object
mp = MedialPDE(5, 5, 10)
mp.GenerateSampleModel()
mp.Solve()

# Load an image from file
img = Image3f()
img.LoadFromFile("hippo_speed.hdr")

# Compute the match between the image and the model
xMatch = mp.ComputeImageMatch(img)
print "Image Match = ", xMatch
