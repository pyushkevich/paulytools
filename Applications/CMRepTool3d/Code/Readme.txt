##############################################
# CM-Reps Command Reference (F2 for console) #
##############################################


help                            This command

load bin FILE                   Load binary image
load [gray|grey] FILE           Load grey image
load [distance|dt] FILE         Load distance transform
load [spline|sp] FILE           Load m-rep spline
load edge FILE                  Load boundary mesh

chdir DIR                       Change current directory

save [image|im] FILE            Save current image
save [spline|sp] FILE           Save currnent m-rep spline
save [pov|povray] FILE          Save spline as povray scene
save [obj|wave] FILE            Save spline as OBJ file

[distanceMap|distmap]           Compute distance map on binary image

convert image [binary|bin]      Threshold image

mask [inside|in] ADD SCALE      Linear transform on DT inside
mask [outside|out] ADD SCALE    Linear transform on DT outside 

blur SCALE                      Blur with a gaussian

list [images|i]                 List images in curr. directory
list [s|spl|splines|spline]     List splines in curr. directory

show slice off                  Hide image slicing
show slice pause                Pause moving slice
show slice last                 Show last slice
show slice next                 Show next slice
show slice axis                 Choose next slice axis
show slice [xy|xz|yz]           Choose particular axis

show cloud                      View DT level set

show medial seethru             Transparent medial surface view
show medial [solid|flat]        Opaque medial surface view

show boundary top               Show top implied boundary toggle
show boundary [bot|bottom]      Show bottom implied boundary toggle

show control                    Show control polygon toggle

show map gauss                  Show gaussian curvature overlat
show map off                    Hide overlays
show map next                   Show next available overlay

show edge                       Show loaded edge mesh

new spline                      Create a new spline

scale spline X Y Z              Affine scale spline
scale image X Y Z               Affine scale image voxel

rotate spline X Y Z Angle       Rotate about the origin 
shift spline X Y Z              Translate spline

crop image PAD RESAMPLE         Crop image around non-zero, resample
threshold image LEVEL           Threshold image at intensity level

[match|m]                       Compute spline/image match
rigid                           Run rigid segmentation
[optimize|opt]                  Run non-rigid segmentation
sleep MS                        Wait some time

set PROP VALUE                  Set property
set set                         List properties

test flow
test fastflow
test flowmove
test gl
test trilerp
test crash
test paste

quit                            Exit program
