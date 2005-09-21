from medialpde import *;
from common import *;

# Choose a sample
sample="1042L"

# Create a medial model and load it from a file
m=MedialPDE(8, 12, 32, 80)
m.LoadFromParameterFile(dirWork + "forpca/st" + file + ".ctf06.mpde");

# Load the input image from file


