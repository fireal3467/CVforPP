#!/bin/bash                                                                  
# 
# Create all output results 
#

# Useful shell settings:

# abort the script if a command fails
set -e

# abort the script if an unitialized shell variable is used
set -u

# make sure the code is up to date

pushd src
make
popd

# generate the result pictures

#src/imgpro input/testpattern.jpg output/testpattern_harris.jpg \
#-log

#src/imgpro input/test_C_bridge01.jpg output/test_C_bridge_blend.jpg \
#-matchTranslation input/test_C_bridge02.jpg

#src/imgpro input/test_D_face01.jpg output/test_D_face_match_detection_dlt.jpg \
#-matchTranslation input/test_D_face02.jpg
#
#src/imgpro input/test_E_sitting01.jpg output/test_E_sitting_match_detection_dlt.jpg \
#-matchTranslation input/test_E_sitting02.jpg
#
src/imgpro input/colorA_lowres.jpg output/color_lowres_blend.jpg \
-matchTranslation input/colorB_lowres.jpg

#src/imgpro input/hw6photo01.jpg output/hw6photo_match_detection.jpg \
#-matchTranslation input/hw6photo02.jpg









