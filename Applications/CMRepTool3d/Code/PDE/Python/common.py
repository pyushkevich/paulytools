# Load all the symbols in the library
from medialpde import *

# Set the working directory
dirWork = "/home/pauly/data2005/Stanley/data/"

# Function to generate working images
def MakeImages(dirWork, id):
  """ Blur binary image to floating point """
  # Create a binary image
  bin = BinaryImage()
  bin.LoadFromFile(dirWork + "img/" + id + ".mha")
  
  # Convert into the floating point image with gradient
  img = FloatImage()

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 1.2)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(dirWork + "img/" + id + "_hi","mha")

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 0.6)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(dirWork + "img/" + id + "_med","mha")

  # Apply different levels of blurring to the image
  img.SetToBlurredBinary(bin, 0.2)
  img.SetOutsideValue(-1.0)
  img.SaveToPath(dirWork + "img/" + id + "_low","mha")

# Function to save an m-rep and the meshes too
def SaveMRep(mp, id, spec):
  mp.SaveToParameterFile(dirWork + "mreps/" + id + "." + spec + ".mpde")
  mp.SaveVTKMesh(dirWork + "mreps/" + id + "." + spec + ".med.vtk",
      dirWork + "mreps/" + id + "." + spec + ".bnd.vtk")
