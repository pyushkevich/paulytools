# Load from common file
from common import *
import sys
import os

# Visualize the MPDE in the work path
fnMesh = dirWork + sys.argv[1];
print "Loading MPDE " + fnMesh

# Load the template
mp = MedialPDE(2, 4, 32, 80, 0.5, 2, 2);
mp.LoadFromParameterFile(fnMesh);

# Render
RenderMedialPDE(mp);


