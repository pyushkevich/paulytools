# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 32, 80)
mp.LoadFromParameterFile(dirWork + "avg/average_mrepL_align.mpde");
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromPath(dirWork + "avg/average_hippo_blurred_low","mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToAffine();
mp.SetOptimizerToConjugateGradientDescent(1.0);
mp.RunOptimization(img, 400);

# Save the results
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_affine.mpde");
mp.SaveBYUMesh(dirWork + "avg/average_mrepL_affine.byu");

RenderMedialPDE(mp);
