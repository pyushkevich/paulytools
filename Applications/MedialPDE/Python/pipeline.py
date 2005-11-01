# Load all the symbols in the library
from medialpde import *
from common import *
import sys
import os

# Define the temporary directory
dirMesh = "${MPDE_TEMP_ROOT}/meshdump/"

# Check whether the output of a step already exists
def CheckOutput(id, expid, name):
  target = dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + name + ".mpde";
  if(os.access(target, os.R_OK)):
    print "File " + target + " already exists."
    return True;
  else:
    return False;

# Define a the affine stage
def Stage_MOInertia(id, expid, nmTemplate, nmEnd, imgType, nu, nv):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(id, expid, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(nmTemplate);
  mp.MatchImageByMoments(img, 5);
    
  # Run the optimization
  SaveMRep(mp, id, expid, nmEnd);

# Define a the affine stage
def Stage_AFF_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(id, expid, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToAffine();
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);

# Define a coarse-to-fine stage of the optimization
def Stage_XYZ_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(id, expid, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToCoarseToFine(nu, nv, 0, 0);
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);

# Define a coarse-to-fine stage of the optimization
def Stage_CTF_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(id, expid, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToDeformable();
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);


def Stage_CTF_CG_BM(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(id, expid, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImageWithGradient(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToBoundaryGradient();
  mp.SetOptimizationToDeformable();
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);


# PCA-based optimization
def Stage_PCA_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter, pcaData, nModes):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(id, expid, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);

  # Load the PCA data
  mp.SetPCAMatrix(pcaData["ncu"], pcaData["ncv"], pcaData["matrix"]);

  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToPCA(nModes);
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);

