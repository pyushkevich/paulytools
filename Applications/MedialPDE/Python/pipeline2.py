# Load all the symbols in the library
from medialpde import *
from common2 import *
import sys
import os

# Define the temporary directory
dirMesh = "${MPDE_TEMP_ROOT}/meshdump/"

# Check whether the output of a step already exists
def CheckOutput(expdata, name):
  target = os.path.join(expdata["cmrep"],expdata["id"] + "." +  name + ".mpde");
  if(os.path.isfile(target)):
    print "File " + target + " already exists."
    return True;
  else:
    return False;

# Define a the affine stage
def Stage_MOInertia(expdata, nmTemplate, nmEnd, imgType, nu, nv):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(nmTemplate);
  mp.MatchImageByMoments(img, 5);
    
  # Run the optimization
  SaveMRep(mp, expdata, nmEnd);

# Define a the affine stage
def Stage_AFF_CG_VO(expdata, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
      os.path.join(expdata["cmrep"],expdata["id"] + "." + nmStart + ".mpde"));
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToAffine();
  mp.EnableMeshDump(os.path.join(expdata["dump"], expdata["id"] + "_" + nmEnd), 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, expdata, nmEnd);

# Define a coarse-to-fine stage of the optimization
def Stage_XYZ_CG_VO(expdata, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
    os.path.join(expdata["cmrep"],expdata["id"] + "." + nmStart + ".mpde"));
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToCoarseToFine(nu, nv, 0, 0);
  mp.EnableMeshDump(os.path.join(expdata["dump"], expdata["id"] + "_" + nmEnd), 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, expdata, nmEnd);

# Define a coarse-to-fine stage of the optimization
def Stage_DEF_CG_VO(expdata, nmStart, nmEnd, imgType, nu, nv, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
    os.path.join(expdata["cmrep"],expdata["id"] + "." + nmStart + ".mpde"));
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToDeformable();
  mp.EnableMeshDump(os.path.join(expdata["dump"], expdata["id"] + "_" + nmEnd), 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, expdata, nmEnd);

# Define a coarse-to-fine stage of the optimization
def Stage_CTF_CG_VO(expdata, nmStart, nmEnd, imgType, nu, nv, nuctf, nvctf, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
    os.path.join(expdata["cmrep"],expdata["id"] + "." + nmStart + ".mpde"));
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToCoarseToFine(nuctf, nvctf, nuctf, nvctf);
  mp.EnableMeshDump(os.path.join(expdata["dump"], expdata["id"] + "_" + nmEnd), 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, expdata, nmEnd);



def Stage_CTF_CG_BM(expdata, nmStart, nmEnd, imgType, nu, nv, nuctf, nvctf, nIter):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
    os.path.join(expdata["cmrep"],expdata["id"] + "." + nmStart + ".mpde"));
  mp.SetNumberOfCoefficients(nu, nv);
    
  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToBoundaryGradient();
  mp.SetOptimizationToCoarseToFine(nuctf, nvctf, nuctf, nvctf);
  mp.EnableMeshDump(os.path.join(expdata["dump"], expdata["id"] + "_" + nmEnd), 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, expdata, nmEnd);


# PCA-based optimization
def Stage_PCA_CG_VO(expdata, nmStart, nmEnd, imgType, nu, nv, nIter, pcaData, nModes):
  """Run an optimization stage"""

  # Check if the output already exists
  if(CheckOutput(expdata, nmEnd)):
    return;
    
  # Load the image from file
  img = LoadBlurImage(expdata, imgType);
  
  # Load the m-rep
  mp = MedialPDE(nu, nv, sampling["nu"], sampling["nv"], sampling["cut"], sampling["ncu"], sampling["ncv"]);
  mp.LoadFromParameterFile(
    os.path.join(expdata["cmrep"],expdata["id"] + "." + nmStart + ".mpde"));
  mp.SetNumberOfCoefficients(nu, nv);

  # Load the PCA data
  mp.SetPCAMatrix(pcaData["ncu"], pcaData["ncv"], pcaData["matrix"]);

  # Set up the optimizer
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToPCA(nModes);
  mp.EnableMeshDump(os.path.join(expdata["dump"], expdata["id"] + "_" + nmEnd), 0.01);
  
  # Run the optimization
  mp.RunOptimization(img, nIter);
  SaveMRep(mp, expdata, nmEnd);

