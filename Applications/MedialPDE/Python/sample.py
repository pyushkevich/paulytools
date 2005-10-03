from medialpde import *;
from common import *;
import sys
import os

# Choose one of the hippocampi
id = sys.argv[1]
expid = sys.argv[2]

# Set up the filenames
fnT1 = dirWork + "/t1/proc/" + id[:-1] + ".hdr"
fnMREP = dirWork + "/cmreps/" + expid + "/" + id + ".ctf80.mpde"
fnOut = dirWork + "/t1/sample/" + expid + "/" + id;

# Check that the filenames exist
for fn in fnT1, fnMREP:
  if not os.access(fn, os.R_OK):
    print "Can't access file " + fn
    sys.exit(1)

# The work
m=MedialPDE(8, 12, 64, 160)
m.LoadFromParameterFile(fnMREP);

i1=FloatImage()
i1.LoadFromFile(fnT1)
i1.SetOutsideValue(0.0)

i2=FloatImage()
m.SampleImage(i1, i2, 16)

i2.SaveToPath(fnOut,"img.gz")



