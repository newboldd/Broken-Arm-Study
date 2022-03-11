#!/bin/csh


set subj = $1
set basedir = $4 #/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/$subj
set workbenchdir = /data/heisenberg/data1/mario/TIM_CN4/workbench/bin_linux64
set maskdir = $basedir/subcortical_mask/
set surfdir = $basedir/7112b_fs_LR/
set volume_mask = $maskdir/subcortical_mask_LR_333.nii
set vol_smooth = 2
set surf_smooth = 2.55
set L_mask = $maskdir/L.atlasroi.32k_fs_LR.shape.gii #L.atlasroi_noproj.func.gii
set R_mask = $maskdir/R.atlasroi.32k_fs_LR.shape.gii #R.atlasroi_noproj.func.gii 


set TR = 1.1
set session = $2
set timename = $3
set L_surface = ${session}_L_dil10_32k_fs_LR_smooth2.55.func.gii
set R_surface = ${session}_R_dil10_32k_fs_LR_smooth2.55.func.gii
set cifti_name = ${timename}_LR_surf_subcort_333_32k_fsLR_surfsmooth${surf_smooth}_volsmooth${vol_smooth}.dtseries.nii

niftigz_4dfp -n $timename $timename
set HEMS = "L R"
foreach hem ($HEMS)
        set midsurf = ${surfdir}/Native/${subj}.${hem}.midthickness.native.surf.gii
        set midsurf_LR32k = ${surfdir}/fsaverage_LR32k/${subj}.${hem}.midthickness.32k_fs_LR.surf.gii
        set whitesurf = ${surfdir}/Native/${subj}.${hem}.white.native.surf.gii
        set pialsurf = ${surfdir}/Native/${subj}.${hem}.pial.native.surf.gii
        set nativedefsphere = ${surfdir}/Native/${subj}.${hem}.sphere.reg.reg_LR.native.surf.gii
        set outsphere = ${surfdir}/fsaverage_LR32k/${subj}.${hem}.sphere.32k_fs_LR.surf.gii
	set submask = goodvoxels_indiv/${session}_mask.nii.gz
        
        set surfname = ${timename}_${hem}
        echo 'Map volume to surface'
        wb_command -volume-to-surface-mapping ${timename}.nii.gz $midsurf ${surfname}.func.gii -ribbon-constrained $whitesurf $pialsurf -volume-roi $submask
        echo 'Dilate surface timecourse'
        wb_command -metric-dilate ${surfname}.func.gii $midsurf 10 ${surfname}_dil10.func.gii
        echo 'Deform timecourse to 32k fs_LR'
        wb_command -metric-resample ${surfname}_dil10.func.gii $nativedefsphere $outsphere ADAP_BARY_AREA ${surfname}_dil10_32k_fs_LR.func.gii -area-surfs $midsurf $midsurf_LR32k
        echo 'Smooth surface timecourse'
        wb_command -metric-smoothing $midsurf_LR32k ${surfname}_dil10_32k_fs_LR.func.gii $surf_smooth ${session}_${hem}_dil10_32k_fs_LR_smooth2.55.func.gii
         
        rm ${surfname}.func.gii
        rm ${surfname}_dil10.func.gii
	rm ${surfname}_dil10_32k_fs_LR.func.gii
echo $hem
end

wb_command -volume-smoothing ${timename}.nii.gz $vol_smooth ${timename}_smooth${vol_smooth}_wROI255.nii.gz -roi $volume_mask

caret_command64 -file-convert -format-convert XML_BASE64 $L_surface
caret_command64 -file-convert -format-convert XML_BASE64 $R_surface

wb_command -cifti-create-dense-timeseries ${cifti_name} -volume ${timename}_smooth${vol_smooth}_wROI255.nii.gz $volume_mask -left-metric $L_surface -roi-left $L_mask -right-metric $R_surface -roi-right $R_mask -timestep $TR -timestart 0

rm $L_surface
rm $R_surface
rm ${timename}_smooth${vol_smooth}_wROI255.nii.gz

exit
