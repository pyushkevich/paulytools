# Load all the symbols in the library
from medialpde import *
import sys
import os

# Set the working directory
dirInput = "${MPDE_INPUTDATA_ROOT}/"
dirWork = "${MPDE_DATA_ROOT}/"

# Define the sampling for the cm-rep experiments
sampling = {"nu":32, "nv":80, "cut":0.5, "ncu":2, "ncv":2}

# Function to generate working images
def MakeImages(id):
  """ Blur binary image to floating point """
  # Create a binary image
  bin = BinaryImage()
  bin.LoadFromFile(dirWork + "hippo/imgiso/" + id + ".mha")
  
  # Convert into the floating point image with gradient
  img = FloatImage()

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 1.2)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(dirWork + "hippo/imgblur/" + id + "_hi","mha")

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 0.6)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(dirWork + "hippo/imgblur/" + id + "_med","mha")

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 0.2)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(dirWork + "hippo/imgblur/" + id + "_low","mha")

# Function to load a blurred image, type is hi, med or low
def LoadBlurImage(id, type):
  """ Load blurred image from file """
  img = FloatImage()
  img.LoadFromFile(dirWork + "hippo/imgblur/" + id + "_" + type + ".mha")
  img.SetOutsideValue(-1.0)
  return img

# Function to load a blurred image and gradient, type is hi, med or low
def LoadBlurImageWithGradient(id, type):
  """ Load blurred image from file """
  img = FloatImage()
  img.LoadFromPath(dirWork + "hippo/imgblur/" + id + "_" + type,"mha")
  img.SetOutsideValue(-1.0)
  return img

# Function to save an m-rep and the meshes too
def SaveMRep(mp, id, expid, spec):
  subpath = expid + "/" + id + "/" + id + "." + spec;
  mp.SaveToParameterFile(dirWork + "cmrep/" + subpath + ".mpde")
  mp.SaveVTKMesh(
    dirWork + "vtk/" + subpath + ".med.vtk", 
    dirWork + "vtk/" + subpath + ".bnd.vtk")

# Ensure that a filename can be created
def CheckDir(dir):
  if(not os.access(dir,os.X_OK)):
    os.makedirs(dir)
  return dir
  
