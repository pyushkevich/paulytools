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
def MakeImages(exp, imgInput, force = True):
  """ Blur binary image to floating point """
  # Check if the images exist
  if(not force):
    if( os.path.isfile(os.path.join(exp["blur"],exp["id"]+"_hi.mha")) and
        os.path.isfile(os.path.join(exp["blur"],exp["id"]+"_med.mha")) and
        os.path.isfile(os.path.join(exp["blur"],exp["id"]+"_low.mha"))):
      print "Blurred images already exist"
      return
      
  # Create a binary image
  bin = BinaryImage()
  bin.LoadFromFile(imgInput);
  
  # Convert into the floating point image with gradient
  img = FloatImage()

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 1.2)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(os.path.join(exp["blur"],exp["id"]+"_hi"),"mha");

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 0.6)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(os.path.join(exp["blur"],exp["id"]+"_med"),"mha");

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 0.2)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(os.path.join(exp["blur"],exp["id"]+"_low"),"mha");

# Function to load a blurred image, type is hi, med or low
def LoadBlurImage(exp, type):
  """ Load blurred image from file """
  img = FloatImage()
  img.LoadFromFile( 
      os.path.join(exp["blur"], exp["id"] + "_" + type + ".mha"));
  img.SetOutsideValue(-1.0)
  return img

# Function to load a blurred image and gradient, type is hi, med or low
def LoadBlurImageWithGradient(exp, type):
  """ Load blurred image from file """
  img = FloatImage()
  img.LoadFromPath(
      os.path.join(exp["blur"], exp["id"] + "_" + type), "mha");
  img.SetOutsideValue(-1.0)
  return img

# Function to save an m-rep and the meshes too
def SaveMRep(mp, exp, spec):
  mp.SaveToParameterFile(
      os.path.join(exp["cmrep"], exp["id"] + "." + spec +  ".mpde"))
  mp.SaveVTKMesh(
      os.path.join(exp["vtk"], exp["id"] + "." + spec +  ".med.vtk"),
      os.path.join(exp["vtk"], exp["id"] + "." + spec +  ".bnd.vtk"))

# Ensure that a filename can be created
def CheckDir(dir):
  if(not os.access(dir,os.X_OK)):
    os.makedirs(dir)
  return dir
  
