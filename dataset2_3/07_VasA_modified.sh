#!/bin/bash

# DESCRIPTION
# ###########
# VasA.sh needs to be run on a participant-by-participant basis. It uses a 4D
# file containing the residuals after first-level GLM estimation to obtain a
# participant-specific vascularization map that is then used to scale this
# participant's contrast images - the resulting rescaled contrast images are to
# be used in group-level analyses.
#
# Based on Kazan et al., 2016, NeuroImage and Kazan et al., 2017, MRM.
# Falk Eippert, FMRIB Centre, University of Oxford. 21/03/2017



project_dir=/data/pt_02661_raw/Heatpain  # general project folder
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located

# add FSL to the path
FSLDIR=${fsl_dir}
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

process_subject() {
  	sub=$1
	printf -v sub "%02d" $subject
	echo "subject: " $sub
	sub_dir=$project_dir/derivatives/sub-sspr${sub}/func
	
	# when some runs were not recorded
	if [[ $subject == 25 || $subject == 27 || $subject == 29 || $subject == 6 ]]; then
		num_runs=2
	elif [[ $subject == 31 ]]; then
		num_runs=3
	else
		num_runs=4
	fi
	
	if [ $subject -gt 16 ]; then
		tr=1.123
	else
		tr=1.12
	fi
			
	for run in $(seq 1 $num_runs); do

		echo "Run: " ${run}
		
		for design in block_2mm_susan onset_2mm_susan block_2mm_susan_even_odd onset_2mm_susan_even_odd block_no_smooth onset_no_smooth; do
			echo $design
			cd $sub_dir
			feat_dir=$sub_dir/run-${run}/feat_${design}.feat
			vasa_dir=$sub_dir/run-${run}/vasa_feat_${design}.feat
			
			cp -r ${feat_dir} ${vasa_dir}
			
			# change to directory containing residuals
			cd ${vasa_dir}/stats

			# linearly detrend (and demean) data
			echo "Detrending 4D residuals..."
			if [ -f tmpReg.txt ]; then
				rm -f tmpReg.txt
			fi
			for t in $(seq 1 1 $(fslval res4d.nii.gz dim4)); do
				echo $t >> tmpReg.txt;
			done
			fsl_glm -i res4d.nii.gz -d tmpReg.txt --out_res=tmpDetrend.nii.gz -m ${vasa_dir}/mask.nii --demean

			# calculate power spectral density image of residuals and take square root
			echo "Calculating PSD and amplitude..."
			fslpspec tmpDetrend.nii.gz tmpPSD.nii.gz
			fslmaths tmpPSD.nii.gz -sqrt tmpSquareRootPSD.nii.gz
			fslmaths tmpSquareRootPSD.nii.gz -div $(echo "$(fslval res4d.nii.gz dim4)/2" | bc -l) tmpNormSquareRootPSD.nii.gz

			# get Nyquist, number of frequencies represented and the frequency bin width
			echo "Selecting frequencies between 0.01 and 0.08 Hz..."
			nyquistFreq=$(echo "1/(2*${tr})" | bc -l)
			numFreqBins=$(fslval tmpNormSquareRootPSD.nii.gz dim4)
			binWidth=$(echo "1/(${numFreqBins}*2)" | bc -l)
			
			# find which images correspond to 0.01Hz (or just above) and to 0.08Hz (or just below)
			freqLo=$binWidth
			indexLo=1
			while [ $(echo "$freqLo < 0.01" | bc -l) -eq 1 ]; do 
				freqLo=$(echo "$freqLo+$binWidth" | bc -l)
				indexLo=$(($indexLo+1))
			done
			freqHi=$binWidth
			indexHi=0
			while [ $(echo "$freqHi < 0.08" | bc -l) -eq 1 ]; do 
				freqHi=$(echo "$freqHi+$binWidth" | bc -l)
				indexHi=$(($indexHi+1))
			done
			
			# select frequencies between 0.01 and 0.08
			fslroi tmpNormSquareRootPSD.nii.gz tmpFiltNormSquareRootPSD.nii.gz $(($indexLo-1)) $((($indexHi-$indexLo)+1))

			# calculate mean amplitude within this frequency band (similar to ALFF)
			echo "Obtaining smoothed and unsmoothed VasA maps..."
			fslmaths tmpFiltNormSquareRootPSD.nii.gz -Tmean tmpVASA.nii.gz

			# output vasa map
			cp tmpVASA.nii.gz vasaSmooth0mm.nii.gz

			# scale contrast images with vasa maps
			echo "Applying VasA to contrast images..."
			
			if [[ $design == block* ]]; then
				if [[ $design == *even_odd ]]; then
					num_contrasts=2
				else
					num_contrasts=1
				fi
			else
				if [[ $design == *even_odd ]]; then
					num_contrasts=5
				else
					num_contrasts=3
				fi
			fi
			for c in $(seq 1 $num_contrasts); do
			    conImg=cope${c}.nii.gz
			    fslmaths ${conImg} -div vasaSmooth0mm.nii.gz $(remove_ext $conImg)_vasa.nii.gz
			    rm -f ${conImg}
			    mv $(remove_ext $conImg)_vasa.nii.gz ${conImg}
			done

			# clean up
			echo "Cleaning up..."
			rm -f tmpReg.txt tmp*.nii.gz
			
			
		done
	done
}

# Loop across subjects for data preparation in parallel
for subject in {1..41}; do
	printf -v sub "%02d" $subject
	process_subject $sub & 
done

# Wait for all background processes to finish
wait
echo "All parallel processes completed."

