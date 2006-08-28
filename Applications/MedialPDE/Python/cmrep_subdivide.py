#!/usr/bin/python
# =========================================================
# PROGRAM: cmrep_subdivide
# PURPOSE: Subdivide subdivision surface cm-rep mesh into 
#          a finer level mesh.
# DESCRIPTION:
#   A subdivision surface cm-rep is defined by a coefficient-level
# mesh which is then subdivided n levels to generate an atom-level
# mesh. The purpose of this code is to change the subdividion level
# either on the coefficient level mesh or in the atom level mesh
from cmrep_common import *;

# Define the command line arguments
parser = OptionParser(
    usage = "usage: %prog [options] input.cmrep output.cmrep",
    description=
    '''A subdivision surface cm-rep is defined by a coefficient-level
    mesh which is then subdivided n levels to generate an atom-level
    mesh. The purpose of this code is to change the subdividion level
    either on the coefficient level mesh or in the atom level mesh''')

parser.add_option('-c', nargs=1, type="int", metavar="NLEVELS", default="0",
    help="Number of levels by which to subdivide the coefficient-level mesh")

parser.add_option('-a', nargs=1, type="int", metavar="NLEVELS", default="0",
    help="Number of levels by which to subdivide the atom-level mesh");

# Parse the command line
opts,args = parser.parse_args()

# Make sure that there are two meshes passed in
if(len(args) <> 2):
  parser.error("incorrect number of positional arguments")

# Load a cmrep from file
cmrep = SubdivisionMPDE(args[0])

# Subdivide the cm-rep as requested
print "Subdividing coeffs by", opts.c , "levels; atoms by", opts.a, "levels"
cmrep.SubdivideMeshes(opts.c, opts.a)

# Save the cm-rep to a new file
cmrep.SaveToParameterFile(args[1])
