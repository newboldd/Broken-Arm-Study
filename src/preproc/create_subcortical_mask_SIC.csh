#!/bin/csh
# create_subcortical_mask_SIC
set subject = $1
set freesurfdir = $2
set t4file = $3
set maskdir = $4

#mkdir subcortical_mask
cp $freesurfdir/$subject/mri/wmparc.mgz $maskdir
cd $maskdir
cp /data/nil-bluearc/GMT/Laumann/MSC/MSC02/subcortical_mask/FreeSurferSubcorticalLabelTableLut* $maskdir
cp /data/nil-bluearc/GMT/Laumann/MSC/MSC02/subcortical_mask/*.atlasroi.32k_fs_LR.shape.gii $maskdir

mri_convert -rl $freesurfdir/$subject/mri/rawavg.mgz -rt nearest wmparc.mgz wmparc.nii
niftigz_4dfp -4 wmparc.nii wmparc -N
t4img_4dfp $t4file wmparc wmparc_333 -O333 -n
niftigz_4dfp -n wmparc_333 wmparc_333

wb_command -volume-label-import wmparc_333.nii.gz FreeSurferSubcorticalLabelTableLut_nobrainstem_LR.txt subcortical_mask_LR_333.nii -discard-others -unlabeled-value 0

wb_command -volume-label-import wmparc_333.nii.gz FreeSurferSubcorticalLabelTableLut_nobrainstem_sub_L_cbll_R.txt subcortical_mask_sub_L_cbll_R.nii -discard-others -unlabeled-value 0

wb_command -volume-label-import wmparc_333.nii.gz FreeSurferSubcorticalLabelTableLut_nobrainstem_sub_R_cbll_L.txt subcortical_mask_sub_R_cbll_L.nii -discard-others -unlabeled-value 0

rm wmparc.mgz
