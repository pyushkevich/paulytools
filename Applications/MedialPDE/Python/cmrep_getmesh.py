#!@PYTHON_BINARY@
from cmrep_common import *

# Define the options to this problem
parser = OptionParser(
    usage="usage: %prog [options] input.cmrep medial.vtk boundary.vtk",
    description=
    '''This program creates a VTK mesh corresponding to a medial model.
    It writes a VTK mesh for the medial surface and another for the boundary
    surface.''')

# Process the command line options
(options,args) = parser.parse_args()

# Check that the parameters are accurate
if(len(args) == 2) :
  fnMed = args[1] + ".med.vtk";
  fnBnd = args[1] + ".bnd.vtk";
elif len(args) == 3 :
  fnMed = args[1];
  fnBnd = args[2];
else :
  parser.error("incorrect number of positional arguments")

# Load the cm-rep and save the mesh
print "Saving mesh to " + fnMed + " and " + fnBnd;
cmrep = MedialPDE(args[0])
cmrep.SaveVTKMesh(fnMed, fnBnd)
