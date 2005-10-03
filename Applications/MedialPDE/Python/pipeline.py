# Load all the symbols in the library
from medialpde import *
from common import *
import sys
import os

# Define the temporary directory
dirMesh = "${MPDE_TEMP_ROOT}/meshdump/"

# Define a the affine stage
def Stage_MOInertia(id, expid, nmTemplate, nmEnd, imgType, nu, nv):
  """Run an optimization stage"""
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, 32, 80);
  mp.LoadFromParameterFile(nmTemplate);
  mp.MatchImageByMoments(img, 5);
    
  # Run the optimization
  SaveMRep(mp, id, expid, nmEnd);

# Define a the affine stage
def Stage_AFF_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, 32, 80);
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
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, 32, 80);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToDeformable(1.0, 0.0);
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);

# Define a coarse-to-fine stage of the optimization
def Stage_CTF_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, 32, 80);
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


# PCA-based optimization
def Stage_PCA_CG_VO(id, expid, nmStart, nmEnd, imgType, nu, nv, nIter, pcaData, nModes):
  """Run an optimization stage"""
    
  # Load the image from file
  img = LoadBlurImage(id, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, 32, 80);
  mp.LoadFromParameterFile(
      dirWork + "cmrep/" + expid + "/" + id + "/" + id + "." + nmStart + ".mpde");
  mp.SetNumberOfCoefficients(nu, nv);

  # Load the PCA data
  mp.SetPCAMatrix(pcaData["matrix"], pcaData["ncu"], pcaData["ncv"]);

  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToPCA(nModes);
  mp.EnableMeshDump(dirMesh + expid + "/" + id + "/" + nmEnd, 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, id, expid, nmEnd);

