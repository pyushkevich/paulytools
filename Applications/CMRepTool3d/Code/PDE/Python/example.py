# Load all the symbols in the library
from medialpde import *

# Set the working directory
dirWork = "/home/pauly/data2005/Stanley/data/"
 
# Create a medial PDE object
mp = MedialPDE(3, 5, 10)
mp.LoadFromDiscreteMRep("/tmp/surf01.txt", 10)
mp.Solve()

# Compute a moment match
bin = BinaryImage()
bin.LoadFromFile(dirWork + "avg/average_hippoL.mha");
mp.MatchImageByMoments(bin)

# Save the m-rep
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_01.byu");
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_01.mpde");

# Convert the binary image into a floating point and compute match
img = FloatImage();
img.SetToBlurredBinary(bin, 1.0);
img.SetOutsideValue(-1.0);

# Save the image and the gradient
img.SaveToFile(dirWork + "avg/average_hippo_blurred.mha");
img.SaveGradientToFile(0, dirWork + "avg/average_hippo_blurred.dx.mha");
img.SaveGradientToFile(1, dirWork + "avg/average_hippo_blurred.dy.mha");
img.SaveGradientToFile(2, dirWork + "avg/average_hippo_blurred.dz.mha");

print "Image Match = ", mp.ComputeImageMatch(img)

RenderMedialPDE(mp)
