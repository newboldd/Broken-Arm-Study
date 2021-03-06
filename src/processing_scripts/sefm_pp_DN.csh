#!/bin/csh

#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/sefm_pp_AZS.csh,v 1.3 2018/02/11 01:18:19 aaron Exp $
#$Log: sefm_pp_AZS.csh,v $
#Revision 1.3  2018/02/11 01:18:19  aaron
#Automatically determines ped
#
# Revision 1.2  2017/12/26  23:42:34  avi
# automatically correct odd slice-count sefm images
#
# Revision 1.1  2017/07/07  05:50:46  avi
# Initial revision
#
set idstr = '$Id: sefm_pp_AZS.csh,v 1.3 2018/02/11 01:18:19 aaron Exp $'
echo $idstr
set program = $0; set program = $program:t

set mcvbin = /data/nil-bluearc/GMT/MRIConvert/bin

if ( ${#argv} < 1 || ${#argv} > 2 ) then
	echo "Usage:	"$program" <params file> [instructions file]"
	echo " e.g.,	"$program" PSQ0001_s1.params ../PSQ_study.params"
	exit 1
endif
set prmfile = $1
echo "prmfile="$prmfile

if (! -e $prmfile) then
	echo $program": "$prmfile not found
	exit -1
endif
source $prmfile
if (${#argv} > 1) then
	set instructions = $2
	if (! -e $instructions) then
		echo $program": "$instructions not found
		exit -1
	endif
	cat $instructions
	source $instructions
endif

####################
# check input params
####################
if ( ! ${?sefm} ) then
	echo "ERROR: sefm must be specified in one of the params files"
	exit 1
endif

if ( ! ${?patid} ) then
	echo "ERROR: patid must be specified in one of the params files"
	exit 1
endif

if ( ${sorted} == 0 ) then
	echo "ERROR: for now, the DICOMs must be sorted"
	exit 1
endif

#if (-e sefm) /bin/rm -rf sefm
#mkdir sefm
pushd sefm

@ status = 0

#################################################
# make datain file assuming all vols are AP or PA
#################################################
if ( -e datain.txt ) /bin/rm datain.txt
touch datain.txt
set geometry = ""
@ i = 1
while ( $i <= ${#sefm} )
	set readout_time_ms = `jq '.TotalReadoutTime' $patid"_sefm"${i}.json`
	set readout_time_sec = `echo $readout_time_ms | awk '{rt_sec = $1 / 1000; printf "%.6f", rt_sec}'`
	@ num_vols = 3
	
	set geostring = `jq '.EchoTrainLength' $patid"_sefm"${i}.json`
	set geostring = "${geostring}x"
	set geostring = "${geostring}"`jq '.PhaseEncodingSteps' $patid"_sefm"${i}.json`
	set geostring = "${geostring}x"
	set geostring = "${geostring}"`jq '.SliceTiming' $patid"_sefm"${i}.json | wc -l | awk '{print $1-2}'`
	set geometry = ( $geometry $geostring )
	
	@ direction = `jq '.PhaseEncodingDirection' $patid"_sefm"${i}.json | sed 's/j//' | sed 's/"//g'`1
	@ j = 1
	while ( $j <= $num_vols )
		echo "0 "$direction" 0 "$readout_time_sec >> datain.txt
		@ j++
	end
	
	@ i++
end
#################################
# check for inconsistent geometry
#################################
@ i = 1
while ( $i <= ${#sefm} )
	if ( $geometry[1] != $geometry[$i] ) then
		echo $program":" geometry mismatch. resampling sefm$i into sefm1 space
		echo sefm1  geometry = $geometry[1]
		echo sefm$i geometry = $geometry[$i]
		t4_ident ident_t4
		tail -4 ident_t4 >! ident.mat
		flirt -in ${patid}"_sefm"$i.nii -ref ${patid}"_sefm1".nii -out ${patid}"_sefm"$i.nii -applyxfm -init ident.mat
		/bin/rm ident.mat ident_t4
	endif
	@ i++
end
################################
# merge AP and PA into one image
################################
set epi_img = `ls $patid"_sefm"*.nii.gz`
echo    fslmerge -t $patid"_sefm_epi" $epi_img
	    fslmerge -t $patid"_sefm_epi" $epi_img
if ($status) exit $status

####################################################################################
# ensure that ${patid}_sefm_epi matrix is compatible with topup default requirements
####################################################################################
fslinfo ${patid}_sefm_epi >! ${patid}_sefm_epi_info.txt
@ nx = `cat ${patid}_sefm_epi_info.txt | gawk '/^dim1/{print $NF;}'`
if ($nx % 2) @ status = -1
@ ny = `cat ${patid}_sefm_epi_info.txt | gawk '/^dim2/{print $NF;}'`
if ($ny % 2) @ status = -1
if ($status) then
	echo $program":" ${patid}_sefm_epi has topup-incompatible odd slice dimensions
	exit -1
endif
@ nz = `cat ${patid}_sefm_epi_info.txt | gawk '/^dim3/{print $NF;}'`
if ($nz % 2) then
	echo $program":" ${patid}_sefm_epi" has topup-incompatible odd slice count - trimming off top slice..."
	cp ${patid}_sefm_epi.nii.gz   ${patid}_sefm_epi1.nii.gz
	mv ${patid}_sefm_epi_info.txt ${patid}_sefm_epi1_info.txt
	echo $patid $nx $ny $nz | gawk '{printf ("fslroi %s_sefm_epi1 %s_sefm_epi 0 %d 0 %d 0 %d\n", $1, $1, $2, $3, $4 - 1);}' >! fslroi.cmd
	cat    fslroi.cmd
	source fslroi.cmd
	if ($status) exit $status
endif

###############################
# use topup to derive field map
###############################
date
echo    topup --imain=$patid"_sefm_epi" --datain=datain.txt --config=b02b0.cnf --fout=$patid"_sefm_pha" --iout=$patid"_sefm_epi_unwarped" --out=$patid"_sefm_topup"
        topup --imain=$patid"_sefm_epi" --datain=datain.txt --config=b02b0.cnf --fout=$patid"_sefm_pha" --iout=$patid"_sefm_epi_unwarped" --out=$patid"_sefm_topup"
	if ($status) exit $status
date

###########################################
# make mag image and use bet to skull strip
###########################################
fslmaths	$patid"_sefm_epi_unwarped" -Tmean $patid"_sefm_mag"
bet2		$patid"_sefm_mag" $patid"_sefm_mag_brain"

chmod 664 *
popd
exit 0
