#!@PYTHON_BINARY@
# =========================================================
# PROGRAM: cmrep_brute2pde
# PURPOSE: Convert a brute-force medial mesh to a PDE-based
#          medial mesh. 
# DESCRIPTION:
from cmrep_common import *;

# Define the command line arguments
parser = OptionParser(
    usage = "usage: %prog [options] input.cmrep output.cmrep",
    description=
    '''This program converts a BruteForce cm-rep to a PDE-based one''');

# Parse the command line
opts,args = parser.parse_args()

# Make sure that there are two meshes passed in
if(len(args) <> 2):
  parser.error("incorrect number of positional arguments")

# Load a cmrep from file
cmrep = SubdivisionMPDE(args[0])

# Subdivide the cm-rep as requested
cmrep.BruteForceToPDE();

# Save the cm-rep to a new file
cmrep.SaveToParameterFile(args[1])
