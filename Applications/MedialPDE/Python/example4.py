# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 40, 80)
mp.LoadFromParameterFile(dirWork + "mreps/st1006L.ctf01.mpde")
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromFile(dirWork + "img/st1006L_med.mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToDeformable(0.5, 0.5);
mp.SetOptimizerToConjugateGradientDescent(1.0);
mp.SetMatchToVolumeOverlap();
mp.EnableMeshDump("/tmp/meshdump/step3",0.001);
mp.RunOptimization(img, 400);

# Save the results
mp.SaveToParameterFile(dirWork + "mreps/st1006L.ctf02.mpde");
mp.SaveVTKMesh(
    dirWork + "mreps/st1006L.ctf02.med.vtk",
    dirWork + "mreps/st1006L.ctf02.bnd.vtk");
