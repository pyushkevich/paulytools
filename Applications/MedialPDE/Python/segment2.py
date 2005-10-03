# Load all the symbols in the library
from pipeline import *

# Choose one of the hippocampi
id = sys.argv[1]
expid = sys.argv[2]

# Set the mesh dump directory
dirMesh = dirWork + "tmp/meshdump/" + expid + "/" + id;

# Set the template to be the initial shape
fnTemplate = dirInput + "init/init.mpde"

# Align template by moments of inertia
Stage_MOInertia(id, expid, fnTemplate, "align","med", 2, 4)

# Align template by moments of inertia
Stage_MOInertia(id, expid, "align",            "med", 2, 4)

# Affine stage
Stage_AFF_CG_VO(id, expid, "align",  "affine", "med", 2, 4, 400);

# Subsequent stages: 2x4 grid
Stage_XYZ_CG_VO(id, expid, "affine", "xyz24",  "med", 2, 4, 800);
Stage_CTF_CG_VO(id, expid, "xyz24",  "med24",  "med", 2, 4, 800);
Stage_CTF_CG_VO(id, expid, "med24",  "ctf24",  "low", 2, 4, 600);

# Subsequent stages: finer grids
Stage_CTF_CG_VO(id, expid, "ctf24",  "ctf35",  "low", 3, 5, 600);
Stage_CTF_CG_VO(id, expid, "ctf35",  "ctf46",  "low", 4, 6, 600);
Stage_CTF_CG_VO(id, expid, "ctf46",  "ctf57",  "low", 5, 7, 600);
Stage_CTF_CG_VO(id, expid, "ctf57",  "ctf68",  "low", 6, 8, 600);
Stage_CTF_CG_VO(id, expid, "ctf68",  "ctf79",  "low", 7, 9, 600);
Stage_CTF_CG_VO(id, expid, "ctf79",  "ctf80",  "low", 8, 10, 600);
