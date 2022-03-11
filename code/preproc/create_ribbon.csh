#!/bin/csh
echo "\n START: CreateRibbon"

set CARET7DIR = /data/heisenberg/data1/mario/wb_dir/workbench/bin_rh_linux64/
set basedir = `pwd`
set freedir = $basedir/fs5.3_native_default/
set pipedir_config = /data/heisenberg/data1/mario/scripts/Pipelines-master/global/config/
set freelabels = ${pipedir_config}/FreeSurferAllLut.txt
set subject = $1

set LeftGreyRibbonValue="3"
set LeftWhiteMaskValue="2"
set RightGreyRibbonValue="42"
set RightWhiteMaskValue="41"

foreach hem ( L R )
  if  ( $hem == "L" ) then
    set GreyRibbonValue="$LeftGreyRibbonValue"
    set WhiteMaskValue="$LeftWhiteMaskValue"
  endif

  if  ( $hem == "R" ) then
    set GreyRibbonValue="$RightGreyRibbonValue"
    set WhiteMaskValue="$RightWhiteMaskValue"
  endif
   
  set T1image = ${basedir}/T1/${subject}_mpr_debias_avgT_333_t88.nii.gz
  niftigz_4dfp -n $T1image:r:r $T1image:r:r
  set nativedir = $basedir/7112b_fs_LR/Native/
  set outputdir = $basedir/7112b_fs_LR/Ribbon/
  echo "${outputdir}"
  mkdir "${outputdir}"
  echo ${CARET7DIR}/wb_command -create-signed-distance-volume ${nativedir}/${subject}.${hem}.white.native.surf.gii ${T1image} ${outputdir}/${subject}.${hem}.white.native.nii.gz 
  ${CARET7DIR}/wb_command -create-signed-distance-volume ${nativedir}/${subject}.${hem}.white.native.surf.gii ${T1image} ${outputdir}/${subject}.${hem}.white.native.nii.gz 
  ${CARET7DIR}/wb_command -create-signed-distance-volume ${nativedir}/${subject}.${hem}.pial.native.surf.gii ${T1image} ${outputdir}/${subject}.${hem}.pial.native.nii.gz 

  fslmaths ${outputdir}/${subject}.${hem}.white.native.nii.gz -thr 0 -bin -mul 255 ${outputdir}/${subject}.${hem}.white_thr0.native.nii.gz
  fslmaths ${outputdir}/${subject}.${hem}.white_thr0.native.nii.gz -bin ${outputdir}/${subject}.${hem}.white_thr0.native.nii.gz
  fslmaths ${outputdir}/${subject}.${hem}.pial.native.nii.gz -uthr 0 -abs -bin -mul 255 ${outputdir}/${subject}.${hem}.pial_uthr0.native.nii.gz
  fslmaths ${outputdir}/${subject}.${hem}.pial_uthr0.native.nii.gz -bin ${outputdir}/${subject}.${hem}.pial_uthr0.native.nii.gz
  fslmaths ${outputdir}/${subject}.${hem}.pial_uthr0.native.nii.gz -mas ${outputdir}/${subject}.${hem}.white_thr0.native.nii.gz -mul 255 ${outputdir}/${subject}.${hem}.ribbon_333.nii.gz
  fslmaths ${outputdir}/${subject}.${hem}.ribbon_333.nii.gz -bin -mul $GreyRibbonValue ${outputdir}/${subject}.${hem}.ribbon_333.nii.gz
  #fslmaths ${outputdir}/${subject}.${hem}.white.native.nii.gz -uthr 0 -abs -bin -mul 255 ${outputdir}/${subject}.${hem}.white_uthr0.native.nii.gz
  #fslmaths ${outputdir}/${subject}.${hem}.white_uthr0.native.nii.gz -bin ${outputdir}/${subject}.${hem}.white_uthr0.native.nii.gz
  #fslmaths ${outputdir}/${subject}.${hem}.white_uthr0.native.nii.gz -mul $WhiteMaskValue ${outputdir}/${subject}.${hem}.white_mask.native.nii.gz
  #fslmaths ${outputdir}/${subject}.${hem}.ribbon.nii.gz -add ${outputdir}/${subject}.${hem}.white_mask.native.nii.gz ${outputdir}/${subject}.${hem}.ribbon.nii.gz
  rm ${outputdir}/${subject}.${hem}.white.native.nii.gz ${outputdir}/${subject}.${hem}.white_thr0.native.nii.gz ${outputdir}/${subject}.${hem}.pial.native.nii.gz ${outputdir}/${subject}.${hem}.pial_uthr0.native.nii.gz ${outputdir}/${subject}.${hem}.white_uthr0.native.nii.gz ${outputdir}/${subject}.${hem}.white_mask.native.nii.gz
end

fslmaths ${outputdir}/${subject}.L.ribbon_333.nii.gz -add ${outputdir}/${subject}.R.ribbon_333.nii.gz ${outputdir}/ribbon_333.nii.gz
#rm ${outputdir}/${subject}.L.ribbon.nii.gz ${outputdir}/${subject}.R.ribbon.nii.gz
${CARET7DIR}/wb_command -volume-label-import ${outputdir}/ribbon_333.nii.gz ${freelabels} ${outputdir}/ribbon_333.nii.gz -discard-others -unlabeled-value 0

echo -e "\n END: CreateRibbon"

