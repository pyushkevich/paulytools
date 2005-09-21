from medialpde import *;
from common import *;
import sys
import os

# Choose one of the hippocampi
id = sys.argv[1]

# Set up the filenames
fnMREP = dirWork + "/cmrep/" + id + "/" + id + ".ctf02.mpde"

# Check that the filenames exist
if not os.access(fnMREP, os.R_OK):
  print "Can't access file " + fnMREP
  sys.exit(1)

# The work
m=MedialPDE(8,12,32,80)
m.LoadFromParameterFile(fnMREP);
for i in range(16):
  print "Testing for component", i
  m.TestDerivativeComputation( i );


