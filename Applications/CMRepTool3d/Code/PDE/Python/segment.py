# Load all the symbols in the library
from medialpde import *
from common import *

# Choose one of the hippocampi
id = "st1006L"

# Create the necessary images
MakeImages(dirWork, id)

# Create a medial PDE object
mp = MedialPDE(2, 4, 20, 40)
mp.LoadFromDiscreteMRep(dirWork + "init/surf01.txt", -0.5)
mp.Solve()

# Load the coarsest image
img = FloatImage()
img.LoadFromFile(dirWork + "img/st1006L_med.mha")
img.SetOutsideValue(-1.0);

# Match the m-rep by moments
mp.MatchImageByMoments(img, 5)
SaveMRep(mp, id, "align")

# Compute affine transform
mp.SetOptimizationToAffine()
mp.SetOptimizerToConjugateGradientDescent(0.1)
mp.SetMatchToVolumeOverlap()
mp.EnableMeshDump("/tmp/meshdump/step1",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "affine")

# Compute low frequency match
mp.SetOptimizationToDeformable(1.0, 0.0);
mp.EnableMeshDump("/tmp/meshdump/step2",0.01)
mp.RunOptimization(img, 800)
SaveMRep(mp, id, "ctf01")

# Compute add the rho component to the match
mp.SetOptimizationToDeformable(1.0, 1.0);
mp.EnableMeshDump("/tmp/meshdump/step3",0.01)
mp.RunOptimization(img, 800)
SaveMRep(mp, id, "ctf02")

# Scale up to 4 by 8 optimization
mp = MedialPDE(4, 6, 20, 40)
mp.LoadFromParameterFile(dirWork + "mreps/" + id + ".ctf02.mpde")
mp.SetOptimizerToConjugateGradientDescent(0.1)
mp.SetMatchToVolumeOverlap()

# Compute match with rho
mp.SetOptimizationToDeformable(1.0, 1.0);
mp.EnableMeshDump("/tmp/meshdump/step4",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "ctf03")

# Scale up to 6 by 10 optimization
mp = MedialPDE(6, 8, 20, 40)
mp.LoadFromParameterFile(dirWork + "mreps/" + id + ".ctf03.mpde")
mp.SetOptimizerToConjugateGradientDescent(0.1)
mp.SetMatchToVolumeOverlap()

# Compute match with rho
mp.SetOptimizationToDeformable(1.0, 1.0);
mp.EnableMeshDump("/tmp/meshdump/step5",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "ctf04")

# Scale up to 6 by 10 optimization
mp = MedialPDE(6, 10, 20, 40)
mp.LoadFromParameterFile(dirWork + "mreps/" + id + ".ctf04.mpde")
mp.SetOptimizerToConjugateGradientDescent(0.1)
mp.SetMatchToVolumeOverlap()

# Compute match with rho
mp.SetOptimizationToDeformable(1.0, 1.0);
mp.EnableMeshDump("/tmp/meshdump/step6",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "ctf05")











