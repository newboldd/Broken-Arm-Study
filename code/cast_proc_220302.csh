#!/bin/csh
###############
# DJN, 07/2018
###############

set basedir = /data/perlman/moochie/analysis/Broken_Arm_Study/
# basedir:
# ├── subject001
# ├── subject002
# ├── ...
# └── code

set subject = BAS001 # subject of interest
set subdir = $basedir/$subject/
# ├── instructions.params
# ├── orig
# ├── preproc
# ├── sessions_orig.txt
# └── sessions.txt

# set laudir = /data/nil-bluearc/GMT/Laumann/MSC/
set REFDIR = `readlink -f ~/refdir`
set AVIDIR = `readlink -f ~/avi_release_dir`
# set FSdir = ${subdir}/Surfaces/fs7.2/

set FSswdir = /usr/local/pkg/freesurfer/bin/
set FREESURFER_HOME = /usr/local/pkg/freesurfer
set FSLdir = /usr/local/pkg/fsl6.0.3/bin/

# dependencies
# add ref dir to path if not alreay on it
set path = ($path $REFDIR)
set path = ($path $AVIDIR)
set path = ($path $FSswdir)
set path = ($path $FSLdir)
source $FREESURFER_HOME/SetUpFreeSurfer.csh

# set scriptdir = /data/nil-bluearc/GMT/Dillan/scripts/

set instruction_file = ${basedir}/code/preproc/instructions.csh
set funcdir = $subdir/Functionals/ # output folder for DCM sort, etc.
set origdir = $subdir/orig/ # contains unzipped CNDA sessions
set seslist = $subdir/sessions.txt # name of new (pretty) session names
set sesnums = `cat $seslist`
set origlist = $subdir/sessions_orig.txt # name of OG (CNDA) session names
set sesnums_orig = `cat $origlist`
set FSdir = ${subdir}/Surfaces/
set struct_proc_script = ${basedir}/code/preproc/Structural_pp_220302.csh

#####
#goto DCM_SORT
goto STRUCT_PARAMS
#goto STRUCT_PROC
#goto FUNC_PARAMS
#goto ATLAS_LINKS

#goto RUN_DVAR_4dfp
#goto SURFACE_PROJECTION_222
#goto FC_GOOD_SURF
#goto FC_GOOD_SURF_TASK
#####

DOWNLOAD:
####################
# Get data from CNDA
####################
set cnda_username = newboldd
set cnda_password =
set project = NP1083

if (! -e $origdir) mkdir $origdir
pushd $origdir

set k = 1
while ( $k <= $#sesnums_orig)
	curl -k -u ${cnda_username}:${cnda_password} "https://cnda.wustl.edu/REST/projects/$project/experiments/${sesnums_orig[$k]}/DIR/SCANS?format=zip&recursive=true" > temp.zip
	unzip temp.zip
	/bin/rm temp.zip
	@ k++
end
popd #out of origdir
exit

DCM_SORT:
################
# Sort dcm files
################
set k = 1

while ( $k <= $#sesnums)

	set ses_orig = $sesnums_orig[$k]
	set ses = $sesnums[$k]
	mkdir ${funcdir}/$ses
	pushd ${funcdir}/$ses
	# pseudo_dcm_sort.csh -i -s ${origdir}/${ses_orig}/scans/
	pseudo_dcm_sort.csh ${origdir}/${ses_orig}/scans
	mv scans.studies.txt ${ses}.studies.txt
	popd
	@ k++
end

exit


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

	# set bolddcm = `cat $k.studies.txt | grep -E "$runtypes[$r]" | grep " $runlengths[$r]" | awk -F " " '{print $1}'`
	# set bolddcm = `echo $bolddcm | tr -d '\n'`
	# if ( $bolddcm > 0) then
	# 	set fstd_dcms = `echo $bolddcm`
	# 	set irun_label = `echo "${runnames[$r]}"`
	# endif
	@ r++

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
	set seqstr = `strings $F | gawk -f $scriptdir/preproc/parse_strings.awk | gawk '{print NR, $1}' | sort -n -k 2,2 | gawk '{printf("%d,", $1);}'`
	echo "set seqstr	= " $seqstr >> $k.params

	# FC params
	echo "set boldruns = (1)" > $k.fcparams

	popd
end
exit

STRUCT_PARAMS:
################
# create struct params
################

set T1_label =
set T2_label =

pushd $funcdir

foreach k ( $sesnums )
	pushd $k
	set T1dcm = `cat $k.studies.txt | grep 'T1MPRAGECor' | grep ' 208' | awk -F " " '{print $1}'`
	set T1dcm = `echo $T1dcm | tr -d '\n'`
	set T1num = $#T1dcm
	set t = 1
	while ( $t <= $T1num )
		if ( $T1dcm[$t] > 0 ) then
			set T1_label = `echo ${T1_label} $funcdir/$k/study$T1dcm[$t]`
		endif
		@ t++
	end

	set T2dcm = `cat $k.studies.txt | grep 'T2w' | grep ' 208' | awk -F " " '{print $1}'`
	set T2dcm = `echo $T2dcm | tr -d '\n'`
	set T2num = $#T2dcm
	set t = 1
	while ( $t <= $T2num )
		if ( $T2dcm[$t] > 0 ) then
			set T2_label = `echo ${T2_label} $funcdir/$k/study$T2dcm[$t]`
		endif
		@ t++
	end
	popd
end

if ( ! -e ${subdir}/Structurals ) mkdir ${subdir}/Structurals

# static struct params
echo "set patid = $subject" > $subdir/Structurals/${subject}.structparams
echo "set structid = ${subject}_struct" >> $subdir/Structurals/${subject}.structparams
echo "set studydir = ${basedir}" >> $subdir/Structurals/${subject}.structparams
echo "set FSdir = ${FSdir}" >> $subdir/Structurals/${subject}.structparams
echo "set PostFSdir = ${FSdir}/FREESURFER_fs_LR" >> $subdir/Structurals/${subject}.structparams

# dynamic struct params
echo "set mprdirs    = ( ${T1_label} )" >> $subdir/Structurals/${subject}.structparams
echo "set t2wdirs    = ( ${T2_label} )" >> $subdir/Structurals/${subject}.structparams

cat $subdir/Structurals/${subject}.structparams

# exit

STRUCT_PROC:
#####################
# Run Laumann proc script
#####################



$struct_proc_script $subdir/Structurals/${subject}.structparams $instruction_file SURFACE_CREATION
# ./Structural_pp_220302.csh /data/perlman/moochie/analysis/Broken_Arm_study/BAS001/Structurals/BAS001.structparams

STRUCT_DCM:
#####################
# Convert dcm to 4dfp
#####################
set structtype = ( T1 T2 )

pushd ${subdir}
source ${subdir}/Structurals/${subject}.structparams
foreach struct ( $structtype )
	mkdir $struct
end

set k = 1
pushd $structtype[1]
while ( $k <= $#T1 )
	set structscan = $T1[$k]
	dcm_to_4dfp -b ${subject}_mpr$k ../$structscan
	@ k++
end
popd
set k = 1
pushd $structtype[2]
while ( $k <= $#T2 )
	set structscan = $T2[$k]
	dcm_to_4dfp -b ${subject}_t2w$k ../$structscan
	@ k++
end
popd
popd

exit

T1_PROC:
###############
# Register T1 to atlas, debias, and average
###############

set structdir = ${subdir}/T1
source ${subdir}/Structurals/${subject}.structparams
pushd ${structdir}
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
	echo ${scriptdir}/apply_debias.csh ${subject}_mpr${k}T
	${scriptdir}/apply_debias.csh ${subject}_mpr${k}T
	niftigz_4dfp -4 ${subject}_mpr${k}T_debias ${subject}_mpr${k}T_debias -N
	@ k++
end

# Mask first T1 for registration
echo bet2 ${subject}_mpr1T_debias ${subject}_mpr1T_debias_bet
bet2 ${subject}_mpr1T_debias.nii.gz ${subject}_mpr1T_debias_bet
niftigz_4dfp -4 ${subject}_mpr1T_debias_bet ${subject}_mpr1T_debias_bet -N

#HERE1:
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

#exit

T2_PROC:
###############
# Register T2 to T1 and average
###############

set structdir = ${basedir}/${subject}/T2
source ${subdir}/Structurals/${subject}.structparams
pushd ${structdir}
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
	echo ${scriptdir}/apply_debias.csh ${subject}_t2w${k}T
	${scriptdir}/apply_debias.csh ${subject}_t2w${k}T
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
exit


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
#exit

GENERIC_PREPROCESS:
###############################################
# Generic preprocessing for dcm_to_4dfp etc...
###############################################
foreach k ( $sesnums )
	pushd $funcdir/$k
	$scriptdir/preproc/cross_bold_dn_180706.csh ${k}.params $instruction_file
	popd
end
exit


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
#exit


APPLY_REGISTER_UWRP_FIRST_SESSION:
###########################################################
# Register func vol ave unwarp to first session epi and resample bold data
###########################################################
set T_epi      = ${subdir}/$sesnums[1]/unwarp/$sesnums[1]_func_vols_ave_uwrp
set T_epi_mask = ${subdir}/$sesnums[1]/unwarp/$sesnums[1]_func_vols_ave_uwrp_mskt
set U          = ${subdir}/$sesnums[1]/unwarp/$sesnums[1]_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
source $instruction_file

# generate mask for first session
pushd ${subdir}/$sesnums[1]/unwarp
echo msktgen_4dfp $sesnums[1]_func_vols_ave_uwrp_mskt -T$REFDIR/TRIO_Y_NDC
msktgen_4dfp $sesnums[1]_func_vols_ave_uwrp -T$REFDIR/TRIO_Y_NDC
popd

# remake single resampled 333 atlas space fMRI volumetric timeseries for first session
set patid = $sesnums[1]
pushd ${subdir}/${patid}
source ${patid}.params
$rsam_cmnd ${patid}.params $instruction_file

set MBstr = _faln_dbnd
set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ k = 1
while ($k <= $#irun)
	echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
	@ k++
end
conc_4dfp ${lst:r}.conc -l$lst
if ($status) exit $status
set format = atlas/${patid}_func_vols.format
if ($status) exit $status
actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
if ($status) exit $status
ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp
var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp
mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp

# Register epi to first session epi, and resample BOLD to atlas
set modes = (0 0 0 0)
@ modes[1] = 2048 + 3 + 256
@ modes[2] = 2048 + 3 + 256 + 4
@ modes[3] = 2048 + 3 + 256 + 4
@ modes[4] = $modes[3]

@ n = $#sesnums
@ i = 2
while ( $i <= $n )
	set patid = $sesnums[$i]
	cd ${subdir}/${patid}
	source ${patid}.params
	pushd unwarp	# into unwarp
	set t4file = ${patid}_func_vols_ave_uwrp_to_$sesnums[1]_func_vols_ave_uwrp_t4
	if ($status) exit $status
	set log =    ${patid}_func_vols_ave_uwrp_to_$sesnums[1]_func_vols_ave_uwrp.log
	date >! $log
	@ k = 1
	while ($k <= $#modes)
	echo	imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		if ($status) exit $status
		@ k++
	end
	t4_mul $t4file $U ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_111 -O111
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_333 -O333
	if ($status) exit $status
	popd		# out of unwarp
	$rsam_cmnd ${patid}.params $instruction_file

	# remake single resampled 333 atlas space fMRI volumetric timeseries
	set MBstr = _faln_dbnd
	set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
	if (-e $lst) /bin/rm $lst
	touch $lst
	@ k = 1
	while ($k <= $#irun)
		echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
		@ k++
	end
	conc_4dfp ${lst:r}.conc -l$lst
	if ($status) exit $status
	set format = atlas/${patid}_func_vols.format
	if ($status) exit $status
	actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
	if ($status) exit $status
	ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
	mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp
	var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
	ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
	mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp
	mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp

	@ i++
end
#goto FC_GOOD_SURF
exit

MEAN_FIELD_MAP_MAKER:
###########################################################
# Make mean field map
###########################################################
${scriptdir}/meanfield_maker_SIC.csh ${subject} ${seslist}
exit

APPLY_REGISTER_UWRPMEAN_FIRST_SESSION:
###########################################################
# Apply mean distortion correction, register func vol ave unwarp to first session and resample BOLD data
###########################################################
source $instruction_file
set FMmean	= $basedir/${subject}/meanfield/phase_rad_unwrap_secs_resolved_on_TRIO_Y_NDC_111

# Apply mean distortion correction to first session, register to t2w, and resample BOLD to atlas
set patid = $sesnums[1]
pushd ${subdir}/${patid}
#if ( -d unwarp_map ) /bin/rm -r unwarp_map
#cp -pr unwarp unwarp_map
source ${patid}.params

set T		= $basedir/${subject}/T2/${subject}_t2w_debias_avgT
set T_mask	= $basedir/${subject}/T2/${subject}_t2w_debias_avgT_bet
set U		= $basedir/${subject}/T2/${subject}_t2w_debias_avgT_to_TRIO_Y_NDC_t4

$uwrp_cmnd -mean atlas/${patid}_func_vols_ave $FMmean atlas/${patid}_func_vols_ave_to_TRIO_Y_NDC_t4 ${dwell} ${ped}
if ($status) exit $status
pushd unwarp	# into unwarp
set t4file = ${patid}_func_vols_ave_uwrpmean_to_${patid}_t2w_t4
cp  ../atlas/${patid}_func_vols_ave_to_${patid}_t2w_t4 $t4file
if ($status) exit $status
set log =    ${patid}_func_vols_ave_uwrpmean_to_t2w_avgT.log
date >! $log
set modes = (0 0 0 0)
@ modes[1] = 3072 + 3 + 256
@ modes[2] = 2048 + 3 + 256
@ modes[3] = 2048 + 3 + 256
@ modes[4] = $modes[3]
@ j = 1
while ($j <= $#modes)
echo	imgreg_4dfp $T $T_mask ${patid}_func_vols_ave_uwrp none $t4file $modes[$j] >> $log
	imgreg_4dfp $T $T_mask ${patid}_func_vols_ave_uwrp none $t4file $modes[$j] >> $log
	if ($status) exit $status
	@ j++
end
t4_mul $t4file $U 	${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
t4img_4dfp 		${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_333 -O333
t4img_4dfp 		${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_111 -O111
msktgen_4dfp ${patid}_func_vols_ave_uwrp -T$REFDIR/TRIO_Y_NDC # Create masked func_vols_ave_unwarp for subsequent registration
if ($status) exit $status
popd		# out of unwarp
$rsam_cmnd ${patid}.params $instruction_file
if ($status) exit $status
if ( -d unwarp_mean ) /bin/rm -r unwarp_mean
mv unwarp unwarp_mean

# remake single resampled 333 atlas space fMRI volumetric timeseries
set MBstr = _faln_dbnd
set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ k = 1
while ($k <= $#irun)
	echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
	@ k++
end
conc_4dfp ${lst:r}.conc -l$lst
if ($status) exit $status
set format = `cat atlas/${patid}_func_vols.format`
if ($status) exit $status
actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
if ($status) exit $status
ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp_mean
var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp_mean
mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp_mean

## Apply mean distortion correction to each session, register epi to first session epi, and resample BOLD to atlas
set T_epi      = ${subdir}/$sesnums[1]/unwarp_mean/$sesnums[1]_func_vols_ave_uwrp
set T_epi_mask = ${subdir}/$sesnums[1]/unwarp_mean/$sesnums[1]_func_vols_ave_uwrp_mskt
set U          = ${subdir}/$sesnums[1]/unwarp_mean/$sesnums[1]_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4

set modes = (0 0 0 0)
@ modes[1] = 2048 + 3 + 256
@ modes[2] = 2048 + 3 + 256 + 4
@ modes[3] = 2048 + 3 + 256 + 4
@ modes[4] = $modes[3]

@ n = $#sesnums
@ i = 2
while ( $i <= $n )
	set patid = $sesnums[$i]
	cd ${subdir}/${patid}
	source ${patid}.params
	$uwrp_cmnd -mean atlas/${patid}_func_vols_ave $FMmean atlas/${patid}_func_vols_ave_to_TRIO_Y_NDC_t4 ${dwell} ${ped}
	if ($status) exit $status
	pushd unwarp	# into unwarp
	set t4file = ${patid}_func_vols_ave_uwrp_to_$sesnums[1]_func_vols_ave_uwrp_t4
	if ($status) exit $status
	set log =    ${patid}_func_vols_ave_uwrpmean_to_$sesnums[1]_func_vols_ave_uwrp.log
	date >! $log
	@ k = 1
	while ($k <= $#modes)
	echo	imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		if ($status) exit $status
		@ k++
	end
	t4_mul $t4file $U ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_111 -O111
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_333 -O333
	if ($status) exit $status
	popd		# out of unwarp
	$rsam_cmnd ${patid}.params $instruction_file
	if ($status) exit $status
	if ( -d unwarp_mean ) /bin/rm -r unwarp_mean
	mv unwarp unwarp_mean

	# remake single resampled 333 atlas space fMRI volumetric timeseries
	set MBstr = _faln_dbnd
	set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
	if (-e $lst) /bin/rm $lst
	touch $lst
	@ k = 1
	while ($k <= $#irun)
		echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
		@ k++
	end
	conc_4dfp ${lst:r}.conc -l$lst
	if ($status) exit $status
	set format = `cat atlas/${patid}_func_vols.format`
	if ($status) exit $status
	actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
	if ($status) exit $status
	ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
	mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp_mean
	var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
	ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
	mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp_mean
	mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp_mean

	@ i++
end
#exit

FCPROCESS:
####################
# RSFC Processing
####################
foreach patid ($sesnums)
	pushd $patid
	echo $scriptdir/preproc/fcMRI_preproc_180730.csh $patid.params $instruction_file
	$scriptdir/preproc/fcMRI_preproc_180730.csh $patid.params $instruction_file
	if ($status) exit $status
	popd
end
exit

CREATE_SURFACES:
#####################
# Create surfaces
#####################
set freesurfbin = /data/heisenberg/data1/freesurfer5.3/bin/
set freesurfdir = $subdir/fs5.3_native_default/

niftigz_4dfp -n ${basedir}/${subject}/T1/${subject}_mpr_debias_avgT ${basedir}/${subject}/T1/${subject}_mpr_debias_avgT
mkdir ${freesurfdir}
${freesurfbin}/recon-all -all -sd ${freesurfdir} -s ${subject} -i ${basedir}/${subject}/T1/${subject}_mpr_debias_avgT.nii.gz

# Create 32k surfaces
set atlas_name = 'native_fs_LR'
set T1dir = $subdir/T1/
set T2dir = $subdir/T2/
set T1name = ${subject}_mpr_debias_avgT_111_t88

set caret_cmd = /data/heisenberg/data1/mario/FSAVG2FSLR_SCRIPTS/PostFreeSurfer/scripts/FreeSurfer2CaretConvertAndRegisterNonlinear.sh
set CaretAtlasFolder = "/data/heisenberg/data1/mario/FSAVG2FSLR_SCRIPTS/global/templates/standard_mesh_atlases"
set caret5_cmd = /data/cn/data1/linux/bin/caret_command64
set caret7_cmd = /data/heisenberg/data1/mario/wb_dir/workbench/bin_rh_linux64/wb_command
set global_scripts = /data/heisenberg/data1/mario/FSAVG2FSLR_SCRIPTS/global/scripts

if (! -e $atlas_name) mkdir $atlas_name
cp -R $freesurfdir/$subject $atlas_name
cp $T1dir/$T1name.4dfp* $atlas_name
chmod 644 $atlas_name/$T1name*
niftigz_4dfp -n $atlas_name/$T1name $atlas_name/$T1name

fslmaths $atlas_name/$T1name.nii.gz -sub $atlas_name/$T1name.nii.gz $atlas_name/zero.nii.gz
fslmerge -t $atlas_name/zero_.nii.gz $atlas_name/zero.nii.gz $atlas_name/zero.nii.gz $atlas_name/zero.nii.gz
mv -f $atlas_name/zero_.nii.gz $atlas_name/zero.nii.gz

$caret_cmd $basedir $subject $subdir/$atlas_name $subdir/$atlas_name Native $freesurfdir/$subject/ $T1name $subdir/$atlas_name/$T1name.nii.gz $T1name $T1name $CaretAtlasFolder 32000 32 $caret5_cmd $caret7_cmd $subdir/$atlas_name/zero $subdir/$atlas_name/zero $T1name $T1name brainmask_fs $caret_cmd:h $global_scripts

rm -r $atlas_name/$subject

# Resample surfaces to atlas space
set nativedir = $basedir/$subject/native_fs_LR/
set resampledir = $basedir/$subject/7112b_fs_LR/
set workbenchdir = /data/heisenberg/data1/mario/wb_dir/workbench/exe_rh_linux64/

#set T1dir = ${basedir}/${subject}/T1
set T1dir = $basedir/${subject}/T1
set Nativevol = ${resampledir}/${subject}_mpr_debias_avgT.nii.gz
set Atlasvol = ${T1dir}/${subject}_mpr_debias_avgT_111_t88.nii.gz

#Copy native surfaces over to start
mkdir ${resampledir}
pushd ${nativedir}
cp -r * ${resampledir}
rm ${Nativevol}
popd
pushd T1
niftigz_4dfp -n $Atlasvol:r:r $Atlasvol:r:r
cp ${Atlasvol} $resampledir
popd

set t4file = ${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4

#Convert t4 to world matrix
set source_vol_4dfp = ${T1dir}/${subject}_mpr1T_debias.4dfp.ifh
set target_vol_4dfp = /data/petsun43/data1/atlas/TRIO_Y_NDC.4dfp.ifh
set source_vol_nii = ${T1dir}/${subject}_mpr1T_debias.nii
pushd T1
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

	caret_command64 -surface-apply-transformation-matrix ${subject}.L.${surface}.164k_fs_LR.coord.gii ${subject}.L.164k_fs_LR.topo.gii ${subject}.L.${surface}.164k_fs_LR.coord.gii -matrix $matrix

	${workbenchdir}/wb_command -surface-apply-affine ${subject}.R.${surface}.164k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.R.${surface}.164k_fs_LR.surf.gii

	caret_command64 -surface-apply-transformation-matrix ${subject}.R.${surface}.164k_fs_LR.coord.gii ${subject}.R.164k_fs_LR.topo.gii ${subject}.R.${surface}.164k_fs_LR.coord.gii -matrix $matrix
end

#Native surface
pushd Native
foreach surface ( $surfaces )
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.L.${surface}.native.surf.gii ${resampledir}/${world_file} ${subject}.L.${surface}.native.surf.gii

	caret_command64 -surface-apply-transformation-matrix ${subject}.L.${surface}.native.coord.gii ${subject}.L.native.topo.gii ${subject}.L.${surface}.native.coord.gii -matrix $matrix

	${workbenchdir}/wb_command -surface-apply-affine ${subject}.R.${surface}.native.surf.gii ${resampledir}/${world_file} ${subject}.R.${surface}.native.surf.gii

	caret_command64 -surface-apply-transformation-matrix ${subject}.R.${surface}.native.coord.gii ${subject}.R.native.topo.gii ${subject}.R.${surface}.native.coord.gii -matrix $matrix
end
popd

#Fsaverage 32k
pushd fsaverage_LR32k
foreach surface ( $surfaces )
	${workbenchdir}/wb_command -surface-apply-affine ${subject}.L.${surface}.32k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.L.${surface}.32k_fs_LR.surf.gii

	caret_command64 -surface-apply-transformation-matrix ${subject}.L.${surface}.32k_fs_LR.coord.gii ${subject}.L.32k_fs_LR.topo.gii ${subject}.L.${surface}.32k_fs_LR.coord.gii -matrix $matrix

	${workbenchdir}/wb_command -surface-apply-affine ${subject}.R.${surface}.32k_fs_LR.surf.gii ${resampledir}/${world_file} ${subject}.R.${surface}.32k_fs_LR.surf.gii

	caret_command64 -surface-apply-transformation-matrix ${subject}.R.${surface}.32k_fs_LR.coord.gii ${subject}.R.32k_fs_LR.topo.gii ${subject}.R.${surface}.32k_fs_LR.coord.gii -matrix $matrix
end
popd

#Remove unchanges coord files
pushd fsaverage
rm *coord*
popd

exit

SUBCORT_MASK:
##################################
# CREATE SUBCORTICAL MASK & RIBBON
##################################
#$scriptdir/surface_projection/create_subcortical_mask.csh $subject
set maskdir = $subdir/subcortical_mask
set t4file = $subdir/T1/${subject}_mpr1T_to_TRIO_Y_NDC_t4
$scriptdir/surface_projection/create_subcortical_mask_SIC.csh $subject $subdir/fs5.3_native_default/ $t4file $maskdir

$scriptdir/surface_projection/create_ribbon.csh $subject
exit

SUBCORT_MASK_222:
##################################
# CREATE SUBCORTICAL MASK & RIBBON
##################################
set maskdir = $subdir/subcortical_mask_222
set t4file = $subdir/T1/${subject}_mpr1T_to_TRIO_Y_NDC_t4
$scriptdir/surface_projection/create_subcortical_mask_SIC_222.csh $subject $subdir/fs5.3_native_default/ $t4file $maskdir

$scriptdir/surface_projection/create_ribbon_222.csh $subject
exit

GOOD_VOXELS:
##########################
# CREATE GOOD VOXELS MASKS
##########################
$scriptdir/surface_projection/create_goodvoxels_mask.csh ${subject} ${seslist} $basedir/7112b_fs_LR/Ribbon
exit

SURFACE_PROJECTION:
################
# create ciftis
################
# tasks
foreach patid ($sesnums)
	pushd $patid
	source ${patid}.params
	@ k = 2
	while ($k <= $#irun)
		set run = $irun[$k]
		pushd bold$run
		$scriptdir/surface_projection/create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b${run}_faln_dbnd_xr3d_uwrp_atl
		popd
		@ k++
	end
	popd
end
# rest
foreach patid ($sesnums)
	pushd $patid/bold1/
	$scriptdir/surface_projection/create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b1_faln_dbnd_xr3d_uwrp_atl_bpss_resid
	popd
end
exit

UPSAMPLE:
##########################
# Create 222 volume files
##########################
mkdir bold1_222
foreach ses (`cat pre_scans.txt`)
	t4img_4dfp none $ses/bold1/${ses}_b1_faln_dbnd_xr3d_uwrp_atl bold1_222/${ses}_b1_faln_dbnd_xr3d_uwrp_atl_222 -O222
end
exit

SURFACE_PROJECTION_222:
##################################
# Make CIFTIs with 222 subcortical
##################################
foreach patid (`cat temp3.txt`)
	$scriptdir/surface_projection/create_goodvoxels_mask_222.csh ${subject} $patid $subdir/7112b_fs_LR/Ribbon_222 $subdir $subdir/bold1_222 $subdir/$patid/bold1/
	pushd bold1_222
	$scriptdir/surface_projection/create_cifti_goodvoxels_SIC_222.csh $subject $patid ${patid}_b1_faln_dbnd_xr3d_uwrp_atl_bpss_resid_222
	popd
end

FC_GOOD_SURF:
###########################################
# combined good_voxels + surface projection
###########################################
foreach patid ($sesnums)
	pushd $patid
	source ${patid}.params

	# FC proc
	echo $scriptdir/preproc/fcMRI_preproc_180730.csh $patid.params $instruction_file
	$scriptdir/preproc/fcMRI_preproc_180730.csh $patid.params $instruction_file
	if ($status) exit $status

	# good voxels
	echo $patid > temp_patid.txt
	$scriptdir/surface_projection/create_goodvoxels_mask_mod.csh ${subject} `pwd`/temp_patid.txt $subdir/7112b_fs_LR/Ribbon $subdir
	rm temp_patid.txt

	# tasks
	@ k = 2
	while ($k <= $#irun)
		set run = $irun[$k]
		pushd bold$run
		$scriptdir/surface_projection/create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b${run}_faln_dbnd_xr3d_uwrp_atl
		popd
		@ k++
	end

	#rest
	pushd bold1
	$scriptdir/surface_projection/create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b1_faln_dbnd_xr3d_uwrp_atl_bpss_resid
	popd
	popd
end
exit

FC_GOOD_SURF_TASK:
###########################################
# combined good_voxels + surface projection
###########################################
foreach patid ($sesnums)
	pushd $patid
	source ${patid}.params

	# FC proc
	echo $scriptdir/preproc/fcMRI_preproc_180730.csh $patid.params $instruction_file
	#$scriptdir/preproc/fcMRI_preproc_180730.csh $patid.params $instruction_file
	if ($status) exit $status

	# using same good voxels mask as for rest

	# tasks
	@ k = 2
	while ($k <= $#irun)
		set run = $irun[$k]
		pushd bold$run
		$scriptdir/surface_projection/create_cifti_goodvoxels_SIC.csh $subject $patid ${patid}_b${run}_faln_dbnd_xr3d_uwrp_atl_bpss_resid
		popd
		@ k++
	end
exit
