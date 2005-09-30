# Load all the symbols in the library
from medialpde import *
from common import *
import sys
import os

# Choose one of the hippocampi
id = sys.argv[1]

# Set the mesh dump directory
dirMesh = dirWork + "tmp/meshdump/" + id;
fnSource = dirWork + "cmrep/" + id + "/" + id + ".ctf02_5.mpde"

# Create a medial PDE object
mp = MedialPDE(2, 4, 32, 80)
print "Loading template from " + fnSource
mp.LoadFromParameterFile(fnSource);

# Load the coarsest image
img = LoadBlurImage(id, "low")

# Run conj grad optimization
mp.SetOptimizerToConjugateGradientDescent(0.1)

# Compute match with rho
mp.SetNumberOfCoefficients(3, 5)
mp.SetMatchToVolumeOverlap()
mp.EnableMeshDump(dirMesh + "/step4",0.01)
mp.RunOptimization(img, 600)
SaveMRep(mp, id, "ctf03")

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
