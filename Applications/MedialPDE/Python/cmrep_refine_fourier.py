#!/usr/bin/python
from cmrep_common import *

# Define the options in this program
parser = OptionParser(
    usage="usage: %prog [options] input.cmrep output.cmrep",
    description=
    '''This program will refine a cmrep defined using Fourier surfaces. It can
    change the number of coefficients defined in a cm-rep and it can also
    change the size of the sampling grid at which image match functions are
    integrated over the boundary/interior of the cm-rep.''')

parser.add_option("-c", "--coefficients", 
    nargs=2, type="int", metavar="NU NV",
    help="Change the number of coefficients in U and V")

parser.add_option("-a", "--atoms",
    nargs=2, type="int", metavar="NU NV",
    help="Change the number of atoms in U and V")

parser.add_option("-e", "--edges", metavar="N",
    type=int,
    help='''Number of additional atoms inserted into the grid
         near the edges of the cm-rep''')

# Process the command line options
(options,args) = parser.parse_args()

# Check that the parameters are accurate
if(len(args) <> 2):
  parser.error("incorrect number of positional arguments")

# Load the input cm-rep
cmrep = CartesianMPDE(args[0])

# Change the number of Fourier coefficients if requested
if(parser.has_option("-c")):
  (nu, nv) = parser.get_option("-c")
  cmrep.SetNumberOfCoefficients(nu, nv)

# Change the grid size if requested
if(parser.has_option("-a")):
  (nu, nv) = parser.get_option("-a")
  cmrep.SetGridSize(nu, nv)

# Save the cmrep
cmrep.SaveToFile(args[1])
