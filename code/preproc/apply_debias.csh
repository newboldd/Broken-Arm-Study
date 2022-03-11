#!/bin/csh
#Calculate bias field and apply correction

set mri_vol = $1

#Estimate bias field
fast -b -B -t 1 -v --nopve ${mri_vol}

fslmaths ${mri_vol} -div ${mri_vol}_bias ${mri_vol}_debias

