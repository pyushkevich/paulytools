# Load all the symbols in the library
from common import *
from medialpde import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 20, 40)
mp.LoadFromDiscreteMRep("/tmp/surf01.txt", -0.5)
mp.Solve()

# Load the image and gradients
img = FloatImage()
img.LoadFromFile(dirWork + "avg/average_hippo_blurred_hi.mha")

# Match the volume to the image
mp.MatchImageByMoments(img, 5)

# Save the m-rep as a mesh and an m-rep file
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_align.byu");
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_align.mpde");

# Compute the boundary image match
print "Image Match = ", mp.ComputeImageMatch(img)

# Draw the m-rep
RenderMedialPDE(mp)
