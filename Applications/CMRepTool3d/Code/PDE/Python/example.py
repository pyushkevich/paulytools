# Load all the symbols in the library
from common import *
from medialpde import *

# Create a medial PDE object
mp = MedialPDE(2, 4, 20, 40)
mp.LoadFromDiscreteMRep("/tmp/surf01.txt", -0.5)
mp.Solve()

# Load the image and gradients
img = FloatImage()
img.LoadFromFile(dirWork + "img/st1006L_med.mha")

# Match the volume to the image
mp.MatchImageByMoments(img, 5)

# Save the m-rep as a mesh and an m-rep file
SaveMRep(mp, "st1006L", "align")

# Compute the boundary image match
mp.SetMatchToVolumeOverlap()
print "Image Match = ", mp.ComputeImageMatch(img)
