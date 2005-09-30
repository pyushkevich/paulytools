# Load all the symbols in the library
from medialpde import *
from common import *
import sys
import os

# Choose one of the hippocampi
id = sys.argv[1]

# Bail out if the data does not exist
#if(not os.access(dirWork + "hippo/imgiso/" + id + ".mha", os.R_OK) ):
#  print "Can not find appropriate image file!"
#  sys.exit(1)

# Set the mesh dump directory
dirMesh = dirWork + "tmp/meshdump/" + id;

# Create the necessary images
# MakeImages(id)

# Create a medial PDE object
mp = MedialPDE(2, 4, 32, 80)
print "Loading template from " + dirInput + "init/init.mpde" 
mp.LoadFromParameterFile(dirInput + "init/init.mpde")

# Load the coarsest image
img = LoadBlurImage(id, "med")

# Match the m-rep by moments
mp.MatchImageByMoments(img, 5)
SaveMRep(mp, id, "align")

# Compute affine transform
mp.SetOptimizationToAffine()
mp.SetOptimizerToConjugateGradientDescent(0.1)
mp.SetMatchToVolumeOverlap()
mp.EnableMeshDump(dirMesh + "/step1",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "affine")

# Compute low frequency match
mp.SetOptimizationToDeformable(1.0, 0.0);
mp.EnableMeshDump(dirMesh + "/step2",0.01)
mp.RunOptimization(img, 800)
SaveMRep(mp, id, "ctf01")

# Compute add the rho component to the match
mp.SetOptimizationToDeformable()
mp.EnableMeshDump(dirMesh + "/step3",0.01)
mp.RunOptimization(img, 800)
SaveMRep(mp, id, "ctf02")

# Compute match with rho
mp.SetNumberOfCoefficients(4, 6)
mp.SetMatchToBoundaryGradient()
mp.EnableMeshDump(dirMesh + "/step4",0.01)
mp.RunOptimization(img, 600)
SaveMRep(mp, id, "ctf03")

# Load the low-resolution image
img = LoadBlurImage(id, "low")

# Scale up to 6 by 10 optimization
mp.SetNumberOfCoefficients(6, 8)
mp.EnableMeshDump(dirMesh + "/step5",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "ctf04")

# Scale up to 6 by 10 optimization
mp.SetNumberOfCoefficients(6, 10)
mp.EnableMeshDump(dirMesh + "/step6",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "ctf05")

# Scale up to 6 by 10 optimization
mp.SetNumberOfCoefficients(8, 12)
mp.EnableMeshDump(dirMesh + "/step7",0.01)
mp.RunOptimization(img, 400)
SaveMRep(mp, id, "ctf06")
