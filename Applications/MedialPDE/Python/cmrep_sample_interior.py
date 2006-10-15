#!@PYTHON_BINARY@
from cmrep_common import *

# Define the options in this program
parser = OptionParser(
    usage="usage: %prog [options] input.cmrep output_points.txt",
    description=
    '''This program extracts points on the interior of a cm-rep and saves
    them to a text file. The points are taken in the u, v, xi coordinate
    system. The u and v samples are the same as the grid in the cm-rep and the
    xi coordinate is sampled according to the -x option. By default, the
    sampling is from -1 to 1 with 0.1 step size. The points are saved to a
    file. Each line in this file has six entries separated by spaces: u v xi x y z''')

parser.add_option("-x", "--xi", 
    nargs=3, type="float", metavar="X0 X1 DX", default=(-1.0, 1.0, 0.1),
    help="Set the sampling in xi (from X0 to X1, step size DX).")

parser.add_option("-i", "--image",
    nargs=1, type="string", metavar="IMAGE",
    help="Image from which to sample intensity values");

# Process the command line options
(options,args) = parser.parse_args()

# Check that the parameters are accurate
if(len(args) <> 2):
  parser.error("incorrect number of positional arguments")

# Get the xi sampling
(xiStart, xiEnd, xiStep) = options.xi

# Load the cm-rep and perform the sampling
cmrep = MedialPDE(args[0])

# See if the image has been supplied
if options.image:
  image = FloatImage(options.image)
  cmrep.SampleInterior(args[1], xiStep, xiStart, xiEnd, image);
else
  cmrep.SampleInterior(args[1], xiStep, xiStart, xiEnd);


