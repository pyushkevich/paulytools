# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(3, 5, 16)
mp.LoadFromParameterFile(dirWork + "avg/average_mrepL_01.mpde");
mp.Solve()

# Convert the binary image into a floating point and compute match
img = FloatImage();
img.LoadFromPath(dirWork + "avg/average_hippo_blurred","mha");
img.SetOutsideValue(-1.0);

# Run the optimization loop
print "Image Match = ", mp.ComputeImageMatch(img)
#mp.GradientDescentOptimization(img, 100, 0.05);
mp.ConjugateGradientOptimization(img, 500, 4.0);

mp.SaveToParameterFile(dirWork + "avg/average_mrepL_02.mpde");
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_02.byu");

# Render the result
# RenderMedialPDE(mp)


