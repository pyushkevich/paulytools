#!@PYTHON_BINARY@
from cmrep_common import *;
import sys
import os

# A method to save PCA to a set of files
def SaveMedialPCA(pca, mpde, path, nModes, nSamples):
  """Save the medial PCA to a matrix, mean and mode files"""

  # Create a substitution rule
  sub = {"path":path}

  # Check all directories
  CheckDir("%(path)s/matrix" % sub)
  CheckDir("%(path)s/mean" % sub)
  CheckDir("%(path)s/vtk" % sub)
  CheckDir("%(path)s/vtk/mean" % sub)
  CheckDir("%(path)s/vtk/modes" % sub)

  # First, save the shape matrix
  pca.ExportShapeMatrix("%(path)s/matrix/shapemat.mat" % sub)

  # Next, save the mean shape
  pca.SetFSLocationToMean();
  pca.GetShapeAtFSLocation(mpde);
  mpde.SaveVTKMesh(
      "%(path)s/vtk/mean/pcamed.mean.vtk" % sub,
      "%(path)s/vtk/mean/pcabnd.mean.vtk" % sub);
  mpde.SaveToParameterFile("%(path)s/mean/mean.cmrep" % sub);

  # Now save the models along the PCA modes
  for i in range(nModes):
    sub["i"] = i;
    CheckDir("%(path)s/vtk/modes/mode%(i)d" % sub)

    # Iterate over the samples in this mode
    # for j in range(-nSamples/2, nSamples/2 + 1, 1):
    for j in range(0, nSamples/2 + 1, 1):
      sub["j"] = j; sub["id"] = j + nSamples/2; sub["z"] = j * 6.0 / nSamples;
      print ("Mode %(i)d, Frame %(j)d, z = %(z)g" % sub);

      # Go to the PCA location
      pca.SetFSLocationToMean()
      pca.SetFSLocation(i, sub["z"]);
      pca.GetShapeAtFSLocation(mpde)

      # Save the VTK mesh
      mpde.SaveVTKMesh(
          "%(path)s/vtk/modes/mode%(i)d/pcamed_m%(i)d_%(id)02d.vtk" % sub,
          "%(path)s/vtk/modes/mode%(i)d/pcabnd_m%(i)d_%(id)02d.vtk" % sub);

# ========================= MAIN CODE ==============================      
# Usage command
parser = OptionParser(
    usage="usage %prog [options] sample1.cmrep ... sampleN.cmrep",
    description=
    '''Compute PCA using a set of cmreps. The cm-reps must have the
    same number of parameters''')

parser.add_option("-o", "--output",
    nargs=1, type="string", metavar="DIR",
    help="Directory where to store PCA data.");

parser.add_option("-r", "--reference",
    nargs=1, type="string", metavar="FILE",
    help="Specify a cm-rep model to use as a reference");

# Check parameters
(opts, args) = parser.parse_args()
if len(args) < 2:
  parser.error("too few samples for PCA")

if not opts.output:
  parser.error("output directory not specified!")

# Create a medial PCA object
pca = MedialPCA();

# Add each of the samples to the PCA object
for i in range(0, len(args)):
  if os.access(args[i], os.R_OK):
    sample = MedialPDE(args[i])
    pca.AddSample(sample);

# One sample is needed as a reference
if opts.reference:
  mpde = MedialPDE(opts.reference)
else:
  mpde = MedialPDE(args[0]);

# Next compute the PCA
print "About to compute PCA!"
pca.ComputePCA(mpde);

# Create the PCA object
SaveMedialPCA(pca, mpde, opts.output, 0, 0)
