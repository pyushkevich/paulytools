# Load all the symbols in the library
from common import *
from medialpde import *

# Create a medial PDE object
mp = MedialPDE(3, 5, 10)
mp.LoadFromDiscreteMRep("/tmp/surf01.txt", -0.25, 10)
mp.Solve()

# Load the image and gradients
img = FloatImage()
img.LoadFromFile(dirWork + "avg/average_hippo_blurred.mha")

# Match the volume to the image
mp.MatchImageByMoments(img, 10)

# Save the m-rep as a mesh and an m-rep file
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_01.byu");
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_01.mpde");

# Compute the boundary image match
print "Image Match = ", mp.ComputeImageMatch(img)

# Draw the m-rep
RenderMedialPDE(mp)
