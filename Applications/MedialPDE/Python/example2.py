# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 32, 80)
mp.LoadFromParameterFile(dirWork + "mreps/st1006L.align.mpde");
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromFile(dirWork + "img/st1006L_med.mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToAffine();
mp.SetOptimizerToConjugateGradientDescent(0.1);
mp.SetMatchToVolumeOverlap();
mp.EnableMeshDump("/tmp/meshdump/step1",0.001);
mp.RunOptimization(img, 400);

# Save the results
mp.SaveToParameterFile(dirWork + "mreps/st1006L.affine.mpde");
mp.SaveVTKMesh(
    dirWork + "mreps/st1006L.affine.med.vtk",
    dirWork + "mreps/st1006L.affine.bnd.vtk");
