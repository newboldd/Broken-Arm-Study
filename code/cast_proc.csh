#!/bin/csh
###############
# DJN, 07/2018
###############

set subject = BAS001
set CNDA_project_name = NP1141 # e.g., NP1083

# Output Paths
set basedir = /data/perlman/moochie/analysis/Broken-Arm-Study/ # TODO would love to make this path relative
# so script can be run from anywhere
set origdir = $basedir/orig_data/$subject/ # download destination for CNDA files
set FSdir = $basedir/freesurfer7.2/ # FS outputs
set procdir = $basedir/proc/ # 

set subdir = $procdir/$subject/
set surfdir = $subdir/Surfaces/ # freesurfer outputs
set funcdir = $subdir/Functionals/ # output folder for DCM sort, etc.

# basedir:
# ├── data_orig
#    ├──subdir1, 2, 3, ...
# ├── proc
#    ├──subdir ...
#       ├── T1
#       ├── T2
#       ├── Surfaces
#       └── Functionals
# ├── freesurfer7.2
#    ├──subdir ...
# ├── results
#    ├──subdir ...
# ├── src
# └── code

# Param files
set instruction_file = $basedir/code/instructions.txt # Edit manually
set structparams_file = $subdir/${subject}.structparams # Will be made automatically
set origlist = $basedir/code/${subject}_sessions_orig.txt # CNDA session names
set seslist = $basedir/code/${subject}_sessions.txt # new (pretty) session names
set sesnums_orig = `cat $origlist`
set sesnums = `cat $seslist`

# Environment
set procSRC = $basedir/src/processing_scripts/
set AVIDIR = /data/nil-bluearc/raichle/lin64-tools/
set FSswdir = /usr/local/pkg/freesurfer/bin/
set FSLdir = /usr/local/pkg/fsl6.0.3/bin/
set path = ($path $procSRC)
set path = ($path $AVIDIR)
set path = ($path $FSswdir)
set path = ($path $FSswdir)

set REFDIR = /data/petsun43/data1/atlas/
set FREESURFER_HOME = /usr/local/pkg/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.csh

# Additional dependencies
set CaretAtlasFolder = /data/heisenberg/data1/mario/FSAVG2FSLR_SCRIPTS/global/templates/standard_mesh_atlases
set caret_cmd = ${procSRC}/FreeSurfer2CaretConvertAndRegisterNonlinear.sh
set caret5_cmd = /data/cn/data1/linux/bin/caret_command64
set caret7_cmd = /data/heisenberg/data1/mario/wb_dir/workbench/bin_rh_linux64/wb_command
set workbenchdir = /data/heisenberg/data1/mario/wb_dir/workbench/exe_rh_linux64/
set caretdir = /data/cn/data1/linux/bin/
set global_scripts = /data/heisenberg/data1/mario/FSAVG2FSLR_SCRIPTS/global/scripts


#######################################################
# CONTENTS
#######################################################
# Uncomment a goto statement to skip to that step
set keep_going = 0

#########################
# Data download + sorting
	#goto DOWNLOAD
	#goto DCM_SORT

#########################
# Structural processing
	#goto STRUCT_PARAMS
	goto STRUCT_DCM
	#goto T1_PROC
	#goto T2_PROC
	#goto ATLAS_LINKS
	#goto FREESURFER
	#goto POSTFS
	#goto SUBCORT_MASK

#########################
# Functional processing
	#goto FUNC_PARAMS
	#goto GENERIC_PREPROCESS
	#goto RUN_DVAR_4dfp
	#goto FC_PROCESS
	#goto SURF_PROJECTION

#######################################################

DOWNLOAD:
####################
# Get data from CNDA
####################
#Get CNDA username and password
echo -n "Enter CNDA username: "
set cnda_username = $<

stty -echo
echo -n "Enter password: "
set cnda_password = $<
stty echo
echo

if (! -e $origdir) mkdir $origdir
pushd $origdir

set k = 1
while ( $k <= $#sesnums_orig)
	curl -k -u ${cnda_username}:${cnda_password} "https://cnda.wustl.edu/REST/projects/$CNDA_project_name/experiments/${sesnums_orig[$k]}/DIR/SCANS?format=zip&recursive=true" > temp.zip
	unzip temp.zip
	/bin/rm temp.zip
	@ k++
end
popd
if (! $keep_going) exit

# TODO
# call downloadCNDASession2.sh with username from above
# then run python cleanup conversion script

DCM_SORT:
################
# Sort dcm files
################

# create initial folders if necessary
if (! -e $procdir) mkdir $procdir
if (! -e $subdir) mkdir $subdir
if (! -e $surfdir) mkdir $surfdir
if (! -e $funcdir) mkdir $funcdir

set k = 1
while ( $k <= $#sesnums)
	set ses_orig = $sesnums_orig[$k]
	set ses = $sesnums[$k]
	mkdir ${funcdir}/$ses
	pushd ${funcdir}/$ses
	pseudo_dcm_sort.csh ${origdir}/${ses_orig}/scans
	mv scans.studies.txt ${ses}.studies.txt
	popd
	@ k++
end
if (! $keep_going) exit

STRUCT_PARAMS:
################
# create struct params
################
pushd $funcdir

set T1_label =
set T2_label =
foreach k ( $sesnums )
	pushd $k
	set T1dcm = `cat $k.studies.txt | grep 'T1MPRAGECor' | grep ' 208' | awk -F " " '{print $1}'`
	set T1dcm = `echo $T1dcm | tr -d '\n'`
	set T1num = $#T1dcm
	set t = 1
	while ( $t <= $T1num )
		if ( $T1dcm[$t] > 0 ) then
			set T1_label = `echo ${T1_label} $k/study$T1dcm[$t]`
		endif
		@ t++
	end
	set T2dcm = `cat $k.studies.txt | grep 'T2w' | grep ' 208' | awk -F " " '{print $1}'`
	set T2dcm = `echo $T2dcm | tr -d '\n'`
	set T2num = $#T2dcm
	set t = 1
	while ( $t <= $T2num )
		if ( $T2dcm[$t] > 0 ) then
			set T2_label = `echo ${T2_label} $k/study$T2dcm[$t]`
		endif
		@ t++
	end
	popd
end
echo "set T1    = ( ${T1_label} )" > $structparams_file
echo "set t2wdirs    = ( ${T2_label} )" >> $structparams_file
cat $structparams_file
if (! $keep_going) exit

STRUCT_DCM:
#####################
# Convert dcm to 4dfp
#####################
set structtype = ( T1 T2 )

pushd ${subdir}
source $structparams_file
foreach struct ( $structtype )
	if (! -e $struct) mkdir $struct
end

set k = 1
pushd $structtype[1]
while ( $k <= $#T1 )
	set structscan = $T1[$k]
	dcm_to_4dfp -b ${subject}_mpr$k $funcdir/$structscan
	@ k++
end
popd

set k = 1
pushd $structtype[2]
while ( $k <= $#T2 )
	set structscan = $T2[$k]
	dcm_to_4dfp -b ${subject}_t2w$k $funcdir/$structscan
	@ k++
end
popd
popd
if (! $keep_going) exit

T1_PROC:
###############
# Register T1 to atlas, debias, and average
###############
source $structparams_file
pushd ${subdir}/T1/
set T1num = $#T1

# Transform from Sagittal to Transverse
set k = 1
while ( $k <= $T1num )
	S2T_4dfp ${subject}_mpr${k} ${subject}_mpr${k}T
	niftigz_4dfp -n ${subject}_mpr${k}T ${subject}_mpr${k}T
	@ k++
end

# Debias and convert back to 4dfp
set k = 1
while ( $k <= $T1num )
	echo apply_debias.csh ${subject}_mpr${k}T
	apply_debias.csh ${subject}_mpr${k}T
	niftigz_4dfp -4 ${subject}_mpr${k}T_debias ${subject}_mpr${k}T_debias -N
	@ k++
end

# Mask first T1 for registration
echo bet2 ${subject}_mpr1T_debias ${subject}_mpr1T_debias_bet
bet2 ${subject}_mpr1T_debias.nii.gz ${subject}_mpr1T_debias_bet
niftigz_4dfp -4 ${subject}_mpr1T_debias_bet ${subject}_mpr1T_debias_bet -N

# Register first T1 to atlas
set modes	= (0 0 0 0 0)
@ modes[1]	= 1024 + 256 + 3
@ modes[2]	= 1024 + 256 + 3
@ modes[3]	= 3072 + 256 + 7
@ modes[4]	= 2048 + 256 + 7
@ modes[5]	= 2048 + 256 + 7
set t4file = ${subject}_mpr1T_to_TRIO_Y_NDC_t4
set ref = /data/petsun43/data1/atlas/TRIO_Y_NDC
set mpr = ${subject}_mpr1T_debias
set mpr_mask = ${subject}_mpr1T_debias_bet
set log = ${subject}_mpr1T_to_TRIO_Y_NDC.log
@ k = 1
while ( $k <= $#modes )
	imgreg_4dfp $ref none $mpr $mpr_mask $t4file $modes[$k] >> $log
	@ k++
end

# Average T1s
set k = 1
set T1scans = ()
while ( $k <= $T1num )
	set T1scans = ( ${T1scans} ${subject}_mpr${k}T_debias )
	@ k++
end
avgmpr_4dfp ${T1scans} ${subject}_mpr_debias_avgT useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
t4imgs_4dfp -s ${subject}_mpr_debias_avgT.lst ${subject}_mpr_debias_avgT -O${subject}_mpr1T_debias

if (! $keep_going) exit

T2_PROC:
###############
# Register T2 to T1 and average
###############
source $structparams_file
pushd ${subdir}/T2/
set T2num = $#T2

ln -s ../T1/${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4 .
foreach e ( img img.rec ifh hdr )
	ln -s ../T1/${subject}_mpr1T_debias.4dfp.$e .
end
# Transform from sagittal to transverse
@ k = 1
while ( $k <= $T2num )
	S2T_4dfp ${subject}_t2w${k} ${subject}_t2w${k}T
	niftigz_4dfp -n ${subject}_t2w${k}T ${subject}_t2w${k}T
	@ k++
end

# Debias and convert back to 4dfp
set k = 1
while ( $k <= $T2num )
	echo apply_debias.csh ${subject}_t2w${k}T
	apply_debias.csh ${subject}_t2w${k}T
	niftigz_4dfp -4 ${subject}_t2w${k}T_debias ${subject}_t2w${k}T_debias -N
	@ k++
end

# Register T2 to T1
t2w2mpr_4dfp ${subject}_mpr1T_debias ${subject}_t2w1T_debias -T/data/petsun43/data1/atlas/TRIO_Y_NDC
set k = 1
set T2scans = ()
while ( $k <= $T2num )
	set T2scans = ( ${T2scans} ${subject}_t2w${k}T_debias )
	@ k++
end
avgmpr_4dfp ${T2scans} ${subject}_t2w_debias_avgT useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
t4imgs_4dfp -s ${subject}_t2w_debias_avgT.lst ${subject}_t2w_debias_avgT -O${subject}_t2w1T_debias
cp ${subject}_t2w1T_debias_to_TRIO_Y_NDC_t4 ${subject}_t2w_debias_avgT_to_TRIO_Y_NDC_t4
cp ${subject}_t2w1T_debias_to_${subject}_mpr1T_debias_t4 ${subject}_t2w_debias_avgT_to_${subject}_mpr1T_debias_t4
t4img_4dfp ${subject}_t2w_debias_avgT_to_${subject}_mpr1T_debias_t4 ${subject}_t2w_debias_avgT ${subject}_t2w_debias_avgT_on_${subject}_mpr1T_debias -O${subject}_mpr1T_debias.4dfp.ifh

# Create masked T2
niftigz_4dfp -n ${subject}_t2w_debias_avgT ${subject}_t2w_debias_avgT
bet2 ${subject}_t2w_debias_avgT ${subject}_t2w_debias_avgT_bet
niftigz_4dfp -4 ${subject}_t2w_debias_avgT_bet ${subject}_t2w_debias_avgT_bet
if (! $keep_going) exit

ATLAS_LINKS:
####################
# create atlas links
####################
foreach k ($sesnums)
	pushd $funcdir/$k
	mkdir atlas
	pushd atlas
	ln -sf ${subdir}/T1/${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4 ./${k}_mpr1_to_TRIO_Y_NDC_t4
	ln -sf ${subdir}/T2/${subject}_t2w_debias_avgT_to_${subject}_mpr1T_debias_t4 ./${k}_t2w_to_${k}_mpr1_t4
	ln -sf ${subdir}/T2/${subject}_t2w_debias_avgT_to_TRIO_Y_NDC_t4 ./${k}_t2w_to_TRIO_Y_NDC_t4
	foreach e ( img ifh img.rec hdr )
		ln -sf ${subdir}/T1/${subject}_mpr_debias_avgT_111_t88.4dfp.$e ./${k}_mpr_n1_111_t88.4dfp.$e
		ln -sf ${subdir}/T1/${subject}_mpr_debias_avgT_222_t88.4dfp.$e ./${k}_mpr_n1_222_t88.4dfp.$e
		ln -sf ${subdir}/T1/${subject}_mpr_debias_avgT_333_t88.4dfp.$e ./${k}_mpr_n1_333_t88.4dfp.$e
		ln -sf ${subdir}/T2/${subject}_t2w_debias_avgT_111_t88.4dfp.$e ./${k}_t2w_111_t88.4dfp.$e
		ln -sf ${subdir}/T2/${subject}_t2w_debias_avgT_222_t88.4dfp.$e ./${k}_t2w_222_t88.4dfp.$e
		ln -sf ${subdir}/T2/${subject}_t2w_debias_avgT_333_t88.4dfp.$e ./${k}_t2w_333_t88.4dfp.$e
		ln -sf ${subdir}/T1/${subject}_mpr_debias_avgT.4dfp.$e ./${k}_mpr1.4dfp.$e
		ln -sf ${subdir}/T2/${subject}_t2w_debias_avgT.4dfp.$e ./${k}_t2w.4dfp.$e
	end
	popd
	popd
end
if (! $keep_going) exit

FREESURFER:
#################
# Run freesurfer
#################
niftigz_4dfp -n $subdir/T1/${subject}_mpr_debias_avgT $subdir/T1/${subject}_mpr_debias_avgT
mkdir ${FSdir}
recon-all -all -sd ${FSdir} -s ${subject} -i ${basedir}/${subject}/T1/${subject}_mpr_debias_avgT.nii.gz

POSTFS:
############################
# Post-freesurfer processing
############################
mkdir $surfdir
pushd $surfdir

#Create 32k surfaces
set atlas_name = 'native_fs_LR'
set T1dir = $subdir/T1/
set T2dir = $subdir/T2/
set T1name = ${subject}_mpr_debias_avgT_111_t88

if (! -e $atlas_name) mkdir $atlas_name
cp -R $FSdir/$subject $atlas_name
cp $T1dir/$T1name.4dfp* $atlas_name
chmod 644 $atlas_name/$T1name*
niftigz_4dfp -n $atlas_name/$T1name $atlas_name/$T1name

fslmaths $atlas_name/$T1name.nii.gz -sub $atlas_name/$T1name.nii.gz $atlas_name/zero.nii.gz
fslmerge -t $atlas_name/zero_.nii.gz $atlas_name/zero.nii.gz $atlas_name/zero.nii.gz $atlas_name/zero.nii.gz
mv -f $atlas_name/zero_.nii.gz $atlas_name/zero.nii.gz

$caret_cmd $basedir $subject $subdir/$atlas_name $subdir/$atlas_name Native $FSdir/$subject/ $T1name $subdir/$atlas_name/$T1name.nii.gz $T1name $T1name $CaretAtlasFolder 32000 32 $caret5_cmd $caret7_cmd $subdir/$atlas_name/zero $subdir/$atlas_name/zero $T1name $T1name brainmask_fs $caret_cmd:h $global_scripts

rm -r $atlas_name/$subject

# Resample surfaces to atlas space
set nativedir = $surfdir/native_fs_LR/
set resampledir = $surfdir/7112b_fs_LR/
set Nativevol = ${resampledir}/${subject}_mpr_debias_avgT.nii.gz
set Atlasvol = ${T1dir}/${subject}_mpr_debias_avgT_111_t88.nii.gz

#Copy native surfaces over to start
mkdir ${resampledir}
pushd ${nativedir}
cp -r * ${resampledir}
rm ${Nativevol}
popd
pushd $T1dir
niftigz_4dfp -n $Atlasvol:r:r $Atlasvol:r:r
cp ${Atlasvol} $resampledir
popd

set t4file = ${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4

#Convert t4 to world matrix
set source_vol_4dfp = ${T1dir}/${subject}_mpr1T_debias.4dfp.ifh
set target_vol_4dfp = /data/petsun43/data1/atlas/TRIO_Y_NDC.4dfp.ifh
set source_vol_nii = ${T1dir}/${subject}_mpr1T_debias.nii
pushd $T1dir
nifti_4dfp -n $source_vol_nii:r $source_vol_nii:r
popd
set target_vol_nii = /data/petsun43/data1/atlas/TRIO_Y_NDC_111.nii
set world_file = ${subject}_mpr1T_debias_to_TRIO_Y_NDC.world
set mat_file = ${subject}_mpr1T_debias_to_TRIO_Y_NDC.mat

pushd ${resampledir}
cp ${T1dir}/${t4file} .
aff_conv 4w ${source_vol_4dfp} ${target_vol_4dfp} ${t4file} ${source_vol_nii} ${target_vol_nii} ${world_file}
aff_conv 4f ${source_vol_4dfp} ${target_vol_4dfp} ${t4file} ${source_vol_nii} ${target_vol_nii} ${mat_file}

set matrix = `cat ${mat_file}`
set surfaces = ( midthickness white pial inflated very_inflated )

#Fsaverage 164k
foreach surface ( $surfaces )
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.L.${surface}.164k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.L.${surface}.164k_fs_LR.surf.gii
	$caretdir/caret_command64 -surface-apply-transformation-matrix ${subject}.L.${surface}.164k_fs_LR.coord.gii ${subject}.L.164k_fs_LR.topo.gii ${subject}.L.${surface}.164k_fs_LR.coord.gii -matrix $matrix
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.R.${surface}.164k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.R.${surface}.164k_fs_LR.surf.gii
	$caretdir/caret_command64 -surface-apply-transformation-matrix ${subject}.R.${surface}.164k_fs_LR.coord.gii ${subject}.R.164k_fs_LR.topo.gii ${subject}.R.${surface}.164k_fs_LR.coord.gii -matrix $matrix
end

#Native surface
pushd Native
foreach surface ( $surfaces )
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.L.${surface}.native.surf.gii ${resampledir}/${world_file} ${subject}.L.${surface}.native.surf.gii
	$caretdir/caret_command64 -surface-apply-transformation-matrix ${subject}.L.${surface}.native.coord.gii ${subject}.L.native.topo.gii ${subject}.L.${surface}.native.coord.gii -matrix $matrix
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.R.${surface}.native.surf.gii ${resampledir}/${world_file} ${subject}.R.${surface}.native.surf.gii
	$caretdir/caret_command64 -surface-apply-transformation-matrix ${subject}.R.${surface}.native.coord.gii ${subject}.R.native.topo.gii ${subject}.R.${surface}.native.coord.gii -matrix $matrix
end
popd

#Fsaverage 32k
pushd fsaverage_LR32k
foreach surface ( $surfaces )
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.L.${surface}.32k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.L.${surface}.32k_fs_LR.surf.gii
	$caretdir/caret_command64 -surface-apply-transformation-matrix ${subject}.L.${surface}.32k_fs_LR.coord.gii ${subject}.L.32k_fs_LR.topo.gii ${subject}.L.${surface}.32k_fs_LR.coord.gii -matrix $matrix
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.R.${surface}.32k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.R.${surface}.32k_fs_LR.surf.gii
	$caretdir/caret_command64 -surface-apply-transformation-matrix ${subject}.R.${surface}.32k_fs_LR.coord.gii ${subject}.R.32k_fs_LR.topo.gii ${subject}.R.${surface}.32k_fs_LR.coord.gii -matrix $matrix
end
popd

#Remove unchanged coord files
pushd fsaverage
rm *coord*
popd

if (! $keep_going) exit

SUBCORT_MASK:
##################################
# CREATE SUBCORTICAL MASK & RIBBON
##################################

# TODO
# May need to change script here to accomodate moving surface atlas directories into surfdir
# may make more sense for maskdir to be in surfaces (FS postprocessing -> cifti)
set maskdir = $subdir/subcortical_mask
set t4file = $subdir/T1/${subject}_mpr1T_to_TRIO_Y_NDC_t4
create_subcortical_mask_SIC.csh $subject $subdir/fs5.3_native_default/ $t4file $maskdir
create_ribbon.csh $subject
if (! $keep_going) exit

FUNC_PARAMS:
################
# create params
################
set runtypes = ( movie movie Motor )
set runlengths = ( 5115 2832 205 )
set runnames = ( movie movie motor )
set runnum = $#runtypes

foreach k ($sesnums)
	pushd $funcdir/$k

	set r = 0
	set fstd_dcms =
	set irun_label =
	set matsize =

	while ( $r <= $runnum )
		set bolddcm = `cat $k.studies.txt | grep -E "$runtypes[$r]" | grep " $runlengths[$r]" | awk -F " " '{print $1}'`
		set bolddcm = `echo $bolddcm | tr -d '\n'`
		set boldnum = $#bolddcm
		if ( $boldnum > 1 ) then
			set t = 1
			while ( $t <= $boldnum )
				if ( $bolddcm[$t] > 0) then
					set fstd_dcms = `echo $fstd_dcms $bolddcm[$t]`
					set irun_label = `echo $irun_label "$runnames[$r]${t}"`

				endif
				@ t++
			end
		else
			if ( $bolddcm > 0) then
				set fstd_dcms = `echo $fstd_dcms $bolddcm`
				set irun_label = `echo $irun_label "$runnames[$r]"`
			endif
		endif
		@ r++
	end
	echo $fstd_dcms
	foreach fs ( $fstd_dcms )
		set temp = `grep -E "$fs.?   epf" $k.studies.txt | awk '{print $2}' | awk -F "_" '{print $2}'`
		set matsize = `echo $matsize $temp`
	end

	echo "set patid    = $k" > $k.params
	echo "set irun     = ($irun_label)" >> $k.params
	echo "set fstd     = ($fstd_dcms)" >> $k.params
	echo "set matrix   = ($matsize)" >> $k.params
	# cat $k.studies.txt | awk 'BEGIN{n=0;};$3~/Field/{s[n]=$1;n++;}END{printf("set sefm\t\t= (");for(i=0;i<n;i++)printf(" %d",s[i]);printf(")\n");}' >> $k.params
	cat $k.studies.txt | awk 'BEGIN{n=0;};$3~/fmap/{s[n]=$1;n++;}END{printf("set sefm\t\t= (");for(i=0;i<n;i++)printf(" %d",s[i]);printf(")\n");}' >> $k.params

	# Sequence string for slice-time correction
	set boldnum = $fstd_dcms[1]
	set F = `ls study$boldnum/*dcm | head -1`
	set seqstr = `strings $F | gawk -f parse_strings.awk | gawk '{print NR, $1}' | sort -n -k 2,2 | gawk '{printf("%d,", $1);}'`
	echo "set seqstr	= " $seqstr >> $k.params

	# FC params
	echo "set boldruns = (movie)" > $k.fcparams

	popd
end
if (! $keep_going) exit

GENERIC_PREPROCESS:
###############################################
# Generic preprocessing for dcm_to_4dfp etc...
###############################################

# TODO
# see if we can include params from .json (dcm2niix output) 
# to write params into instruction_file to be read by cross_bold
foreach k ( $sesnums )
	pushd $funcdir/$k
	cross_bold_dn_180706.csh ${k}.params $instruction_file
	popd
end
if (! $keep_going) exit

RUN_DVAR_4dfp:
#######################################
# run_dvar_4dfp individually on each run
#######################################

foreach k ( $sesnums )

	set patid = $k
	pushd ${subdir}/${patid}
	source ${patid}.params
	foreach	r ( $irun )
		pushd ./bold$r/
		echo ${patid}_b${r}_faln_dbnd_xr3d_norm > ${patid}_b${r}.lst
		conc_4dfp ${patid}_b${r}_faln_dbnd_xr3d_norm -l${patid}_b${r}.lst
		run_dvar_4dfp ${patid}_b${r}_faln_dbnd_xr3d_norm.conc -m../atlas/${patid}_func_vols_ave -n0 -b10 -x8
		popd
	end
	popd
end
if (! $keep_going) exit

FCPROCESS:
####################
# RSFC Processing
####################
foreach patid ($sesnums)
	pushd $patid
	echo fcMRI_preproc_180730.csh $patid.params $instruction_file
	fcMRI_preproc_180730.csh $patid.params $instruction_file
	if ($status) exit $status
	popd
end
if (! $keep_going) exit

SURF_PROJECTION:
###########################################
# combined good_voxels + surface projection
###########################################
# May need to change script here to accomodate moving surface atlas directories into surfdir
foreach patid ($sesnums)
	pushd $patid
	source ${patid}.params

	# good voxels
	echo $patid > temp_patid.txt
	create_goodvoxels_mask_mod.csh ${subject} `pwd`/temp_patid.txt $subdir/7112b_fs_LR/Ribbon $subdir
	rm temp_patid.txt

	# tasks
	@ k = 2
	while ($k <= $#irun)
		set run = $irun[$k]
		pushd bold$run
		create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b${run}_faln_dbnd_xr3d_uwrp_atl
		popd
		@ k++
	end

	#rest
	pushd bold1
	surface_projection/create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b1_faln_dbnd_xr3d_uwrp_atl_bpss_resid
	popd
	popd
end
exit
