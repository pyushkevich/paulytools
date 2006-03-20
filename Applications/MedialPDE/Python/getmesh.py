#! /usr/bin/python
from common import *
import sys
import os

# Visualize the MPDE in the work path
fnMesh = sys.argv[1];
print "Loading MPDE " + fnMesh

# Get the number of samples
su = int(sys.argv[2]);
sv = int(sys.argv[3]);

# Load the template
mp = MedialPDE(2, 4, su, sv, 0.5, 2, 2);
mp.LoadFromParameterFile(fnMesh);

# Render
mp.SaveVTKMesh(sys.argv[4] + ".med.vtk", sys.argv[4] + ".bnd.vtk");


