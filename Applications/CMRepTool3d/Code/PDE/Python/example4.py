# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 40, 80)
mp.LoadFromParameterFile(dirWork + "avg/average_mrepL_ctf01.mpde");
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromPath(dirWork + "avg/average_hippo_blurred_low","mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToDeformable(0.0, 1.0);
mp.SetOptimizerToConjugateGradientDescent(1.0);
mp.SetMatchToVolumeOverlap();
mp.EnableMeshDump("/tmp/meshdump/step3",0.001);
mp.RunOptimization(img, 400);

# Save the results
mp.SaveToParameterFile(dirWork + "avg/average_mrepL_ctf02.mpde");
mp.SaveVTKMesh(
    dirWork + "avg/average_mrepL_ctf01.med.vtk",
    dirWork + "avg/average_mrepL_ctf01.bnd.vtk");
