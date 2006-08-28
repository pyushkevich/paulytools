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

parser.add_option("-e", "--edges", 
    type=int, nargs=2, metavar="NU NV", default=(0, 0),
    help='''Number of additional atoms inserted into the grid
         near the edges of the medial sheet (2 is a good amount) ''')

parser.add_option("-f", "--edgefactor",
    type=float, metavar="X", default=0.5,
    help='''Factor by which to subdivide the atom grid near edges''')

# Process the command line options
(options,args) = parser.parse_args()

# Check that the parameters are accurate
if(len(args) <> 2):
  parser.error("incorrect number of positional arguments")

# Load the input cm-rep
cmrep = CartesianMPDE(args[0])

# Change the number of Fourier coefficients if requested
if options.coefficients:
  (nu, nv) = options.coefficients
  print "Setting number of coefficients to", nu, ",", nv
  cmrep.SetNumberOfCoefficients(nu, nv)

# Change the grid size if requested
if options.atoms:
  (nu, nv) = options.atoms
  (eu, ev) = options.edges
  efac = options.edgefactor
  print "Setting number of atoms to", nu + eu * 2, ",", nv + ev * 2
  cmrep.SetGridSize(nu, nv, eu, ev, efac)

# Save the cmrep
cmrep.SaveToParameterFile(args[1])
