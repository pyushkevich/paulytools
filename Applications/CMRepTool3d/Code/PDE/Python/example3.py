# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 24, 72)
mp.LoadFromParameterFile(dirWork + "avg/average_mrepL_affine.mpde");
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromPath(dirWork + "avg/average_hippo_blurred_hi","mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToDeformable(0.5, 0.0);
mp.SetOptimizerToConjugateGradientDescent(0.1);
mp.RunOptimization(img, 400);

# Save the results
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_ctf01.mpde");
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_ctf01.byu");
