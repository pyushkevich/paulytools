# Load all the symbols in the library
from pipeline import *

# Choose one of the hippocampi
id = sys.argv[1]

# Set the mesh dump directory
dirMesh = dirWork + "tmp/meshdump/" + id;

# Align template by moments of inertia
Stage_MOInertia(id, "align",            "med", 2, 4)

# Affine stage
Stage_AFF_CG_VO(id, "align",  "affine", "med", 2, 4, 400);

# Subsequent stages: 2x4 grid
Stage_XYZ_CG_VO(id, "affine", "xyz24",  "med", 2, 4, 600);
Stage_CTF_CG_VO(id, "xyz24",  "ctf24",  "med", 2, 4, 600);
Stage_CTF_CG_VO(id, "ctf24",  "low24",  "low", 2, 4, 300);

# Subsequent stages: finer grids
Stage_CTF_CG_VO(id, "low24",  "ctf35",  "low", 3, 5, 600);
Stage_CTF_CG_VO(id, "ctf35",  "ctf46",  "low", 4, 6, 600);
Stage_CTF_CG_VO(id, "ctf46",  "ctf57",  "low", 5, 7, 600);
Stage_CTF_CG_VO(id, "ctf57",  "ctf68",  "low", 6, 8, 600);
