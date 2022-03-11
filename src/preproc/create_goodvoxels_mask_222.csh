#!/bin/csh
#set basedir = /data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC03/
set subject = $1
set fsdir = $3 #$basedir/7112b_fs_LR/Ribbon/
set neighsmooth = 5
set factor = 0.5
set session = $2
set subdir = $4
set funcdir = $5
set sesdir = $6


pushd $funcdir
set preproc_runfunc = ${session}_b1_faln_dbnd_xr3d_uwrp_atl_222
set outputdir = $funcdir/goodvoxels_indiv
echo ${outputdir}
mkdir ${outputdir}
set format = `cat $sesdir/${session}_b1_faln_dbnd_xr3d_norm.format`
actmapf_4dfp "${format}" -amean ${preproc_runfunc}
var_4dfp -f"${format}" -s ${preproc_runfunc} 
mv ${preproc_runfunc}_mean.4dfp.* ${outputdir}
mv ${preproc_runfunc}_sd1.4dfp.* ${outputdir}

pushd ${outputdir}
niftigz_4dfp -n ${preproc_runfunc}_mean ${session}_mean
niftigz_4dfp -n ${preproc_runfunc}_sd1 ${session}_sd1
rm ${preproc_runfunc}_mean.4dfp.*
rm ${preproc_runfunc}_sd1.4dfp.*


fslmaths ${outputdir}/${session}_sd1 -div ${outputdir}/${session}_mean ${outputdir}/${session}_cov

fslmaths ${outputdir}/${session}_cov -mas ${fsdir}/ribbon_222.nii.gz ${outputdir}/${session}_cov_ribbon

fslmaths ${outputdir}/${session}_cov_ribbon -div `fslstats ${outputdir}/${session}_cov_ribbon -M` ${outputdir}/${session}_cov_ribbon_norm
fslmaths ${outputdir}/${session}_cov_ribbon_norm -bin -s $neighsmooth ${outputdir}/${session}_SmoothNorm
fslmaths ${outputdir}/${session}_cov_ribbon_norm -s $neighsmooth -div ${outputdir}/${session}_SmoothNorm -dilD ${outputdir}/${session}_cov_ribbon_norm_s${neighsmooth}
fslmaths ${outputdir}/${session}_cov -div `fslstats ${outputdir}/${session}_cov_ribbon -M` -div ${outputdir}/${session}_cov_ribbon_norm_s${neighsmooth} -uthr 1000 ${outputdir}/${session}_cov_norm_modulate
fslmaths ${outputdir}/${session}_cov_norm_modulate -mas ${fsdir}/ribbon_222.nii.gz ${outputdir}/${session}_cov_norm_modulate_ribbon

set STD = `fslstats ${outputdir}/${session}_cov_norm_modulate_ribbon -S`
set MEAN = `fslstats ${outputdir}/${session}_cov_norm_modulate_ribbon -M`

set Lower = `echo "${MEAN} - (${STD} * ${factor})" | bc -l`

set Upper = `echo "${MEAN} + (${STD} * ${factor})" | bc -l`

fslmaths ${outputdir}/${session}_mean -bin ${outputdir}/${session}_mask
fslmaths ${outputdir}/${session}_cov_norm_modulate -thr $Upper -bin -sub ${outputdir}/${session}_mask -mul -1 ${outputdir}/${session}_goodvoxels
popd
popd
