# Load all the symbols in the library
from medialpde import *
from common import *

# Create a medial PDE object
mp = MedialPDE(4, 8, 40, 80)
mp.LoadFromParameterFile(dirWork + "mreps/st1006L.affine.mpde");
mp.Solve()

# Load the highly blurred image
img = FloatImage();
img.LoadFromFile(dirWork + "img/st1006L_med.mha");
img.SetOutsideValue(-1.0);

# Compute the affine transform
mp.SetOptimizationToDeformable(0.5, 0.0);
mp.SetOptimizerToConjugateGradientDescent(0.1);
mp.SetMatchToVolumeOverlap();
mp.EnableMeshDump("/tmp/meshdump/step2",0.001);
mp.RunOptimization(img, 800);

# Save the results
mp.SaveToParameterFile(dirWork + "mreps/st1006L.ctf01.mpde");
mp.SaveVTKMesh(
    dirWork + "mreps/st1006L.ctf01.med.vtk",
    dirWork + "mreps/st1006L.ctf01.bnd.vtk");
