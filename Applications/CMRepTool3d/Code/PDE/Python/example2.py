# Load all the symbols in the library
from medialpde import *

# Set the working directory
dirWork = "/home/pauly/data2005/Stanley/data/"
 
# Create a medial PDE object
mp = MedialPDE(3, 5, 16)
mp.LoadFromParameterFile(dirWork + "avg/average_mrepL_01.mpde");
mp.Solve()

# Convert the binary image into a floating point and compute match
img = FloatImage();
img.LoadFromFile(dirWork + "avg/average_hippo_blurred.mha");
img.LoadGradientFromFile(0, dirWork + "avg/average_hippo_blurred.dx.mha");
img.LoadGradientFromFile(1, dirWork + "avg/average_hippo_blurred.dy.mha");
img.LoadGradientFromFile(2, dirWork + "avg/average_hippo_blurred.dz.mha");
img.SetOutsideValue(-1.0);

# Compute the gradient of the match
print "Image Match = ", mp.ComputeImageMatch(img)
mp.ConjugateGradientOptimization(img, 5000);
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_02.mpde");
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_02.byu");

# Render the result
# RenderMedialPDE(mp)


