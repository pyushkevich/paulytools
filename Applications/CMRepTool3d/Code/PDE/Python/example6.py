# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(6, 8, 40, 80)
mp.LoadFromParameterFile(dirWork + "avg/average_mrepL_ctf03.mpde");
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromPath(dirWork + "avg/average_hippo_blurred_low","mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToDeformable(1.0, 1.0);
mp.SetOptimizerToConjugateGradientDescent(4.0);
mp.SetMatchToBoundaryGradient();
mp.EnableMeshDump("/tmp/meshdump/step5",0.001);
mp.RunOptimization(img, 800);

# Save the results
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_ctf04.mpde");
mp.SaveVTKMesh(
    dirWork + "avg/average_mrepL_ctf04.med.vtk",
    dirWork + "avg/average_mrepL_ctf04.bnd.vtk");
