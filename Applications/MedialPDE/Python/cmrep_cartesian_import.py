#!@PYTHON_BINARY@
from cmrep_common import *

# Define the options in this program
parser = OptionParser(
    usage="usage: %prog [options] input.txt output.cmrep",
    description=
    '''This program creates a Fourier-based Cartesian MPDE model fitted to 
    a set of points from a file. The points have [u v xi x y z rho r] value
    tuples. Optionally, optimization will be used to try to ensure that the
    output model is valid.''')

parser.add_option("-c", "--coefficients", 
    nargs=2, type="int", metavar="NU NV", default=(6,6),
    help="Set the number of coefficients in U and V")

parser.add_option("-a", "--atoms",
    nargs=2, type="int", metavar="NU NV", default=(32,32),
    help="Set the number of atoms in U and V")

parser.add_option("-e", "--edges", 
    type=int, nargs=2, metavar="NU NV", default=(0, 0),
    help='''Number of additional atoms inserted into the grid
         near the edges of the medial sheet (2 is a good amount) ''')

parser.add_option("-f", "--edgefactor",
    type=float, metavar="X", default=0.5,
    help="Factor by which to subdivide the atom grid near edges")

parser.add_option("-r", "--rho",
    type=float, metavar="X", default=0.0,
    help="Set the value of rho to a constant, omit to compute optimal rho");

parser.add_option("-k", "--clean", default=False, action="store_true",
    help="Clean the cm-rep to remove irregularities after fitting")

parser.add_option("-i", "--inflate", default=False, action="store_true",
    help="Inflate the model to optimally fit the radius field");

# Process the command line options
(options,args) = parser.parse_args()

# Check that the parameters are accurate
if(len(args) <> 2):
  parser.error("incorrect number of positional arguments")

# Get the coefficient values
(ncu, ncv) = options.coefficients;
(nau, nav) = options.atoms;
(eu, ev) = options.edges
efac = options.edgefactor;

# Create the initial cm-rep
cmrep = CartesianMPDE(ncu, ncv, nau, nav, efac, eu, ev);
cmrep.ImportFromPointFile(args[0], options.rho, options.clean, options.inflate);
cmrep.Solve();
cmrep.SaveToParameterFile(args[1]);

