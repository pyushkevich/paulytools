#!@PYTHON_BINARY@
from cmrep_common import *

# Load the cm-rep from file
mp = MedialPDE(sys.argv[1])

# Visualize
RenderMedialPDE(mp)
