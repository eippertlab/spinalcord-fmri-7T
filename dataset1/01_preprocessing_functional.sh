#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
preparation=1
thermnoise=1
quickseg=1
# manual correction
moco=0
csf=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02616/data  # general project folder
code_dir=/data/pt_02616/code  # code folder with subfolder helper_functions
helper_dir=$code_dir/helper_functions/

# needed for denoising:
mppca_dir=/data/u_uhorn_software/mppca_denoise-master/  # MPPCA toolbox downloaded from https://github.com/NYU-DiffusionMRI/mppca_denoise

# general processing software
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located
ants_dir=/data/u_uhorn_software/install/  # where your ANTS installation is located (for zshift only)
afni_dir=/afs/cbs.mpg.de/software/.afni/currentversion.debian-bullseye-amd64/debian-bullseye-amd64/  # where your AFNI installation is (in step csf for 3dFWHMx command only)
# (not working for us, start environment using AFNI command)

# add FSL to the path
FSLDIR=${fsl_dir}
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# add SCT to the path
export PATH="${sct_dir}/bin:${PATH}"

# add ANTS to the path
export PATH="${ants_dir}/bin:${PATH}"

# add ANFI to the path
export PATH="${afni_dir}/bin:${PATH}"

for subject in {1..62}; do
	printf -v sub "%03d" $subject
	echo "subject: " $sub
	
	if [[ "$subject" == 29 || "$subject" == 60 ]]; then
		echo "Subject has no functional data"
		continue
	elif [[ "$subject" == 56 ]]; then
		echo "Subject has severe artifacts"
		continue
	fi

	data_dir=$project_dir/raw_data/sub-${sub}/func
	out_dir=$project_dir/derivatives/sub-${sub}/func
	
	mkdir -p $out_dir
	
	# -------------------------------------------------------------------------------------------
	# 0. As preparation sort the data you need into run folders and give easier names
	# -------------------------------------------------------------------------------------------
	if [ $preparation = 1 ]; then
		cd $out_dir
		echo "Sorting data"
		# find all .nii.gz files
		find "$data_dir" -type f -name "*.nii.gz" | while read -r file; do
			# extract the acquisition
			base_name=$(basename "$file")
			acq=$(echo "$base_name" | sed -n 's/.*_acq-\(.*\)_bold.nii.gz/\1/p')
			mkdir -p "$out_dir/$acq"
			cp $file "$out_dir/$acq/bold.nii.gz"
		done
	fi
	

	# -------------------------------------------------------------------------------------------
	# 1. Thermal denoising using MPPCA
	# -------------------------------------------------------------------------------------------
	if [ $thermnoise = 1 ]; then
		for dir in "$out_dir"/*/; do
			cd "$dir"
			file="bold.nii.gz"
			fname=$(basename "$file" | cut -d. -f1)
			echo "Running thermal denoising on file $fname in folder $dir"
			input="$file"
			folderout="$dir"
			
			# shortcut to run matlab function from command line:
			matlab -nosplash -nodesktop -r "addpath('${helper_dir}'); \
						        try; thermalNoiseRemoval('${mppca_dir}', '${fsl_dir}', '${input}', '${fname}', '${folderout}'); catch ME; disp(ME); end; quit;";
			
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 2. Do a quick first segmentation of the cord
	# -------------------------------------------------------------------------------------------
	if [ $quickseg = 1 ]; then
		for dir in "$out_dir"/*/; do
			cd "$dir"
			echo "Segmenting spinal cord"
			# Quick spinal cord segmentation, does not need to be great, just something in every slice --> check
			sct_deepseg_sc -i bold_denoised_mean.nii.gz -c t2s
			
			# Extract tsnr values within cord
			fslmeants -i bold_denoised_tsnr.nii.gz -m bold_denoised_mean_seg.nii.gz -o bold_denoised_tsnr_cord.txt
			fslmeants -i bold_tsnr.nii.gz -m bold_denoised_mean_seg.nii.gz -o bold_tsnr_cord.txt
		done
	fi
	# -------------------------------------------------------------------------------------------
	# afterwards check quality of segmentation, there just needs to be something in every slice
	# manual correction necessary in:
	# sub-004: almost all runs
	# sub-006: SEsosFLEETp3AP
	# sub-007: almost all runs
	# sub-008: sense1FLEETp3, sense1GREp3, sosFLEETp3, sosFLEETp3AP, sosGREp3
	# sub-009: sense1FLEETp3, SEsosFLEETp3PA, sosFLEETp3PA
	# sub-010: almost all runs
	# sub-011: sosFLEETp3PA
	# sub-012: SEsosFLEETp3AP, SEsosFLEETp3PA, sosFLEETp3PA
	# sub-013: sense1FLEETp3, sense1FLEETp4, SEsosFLEETp3AP
	# sub-014: almost all runs
	# sub-015: SEsosFLEETp3AP, sosFLEETp3PA
	# sub-016: sosFLEETp3PA
	# sub-021: almost all runs
	# sub-022: SEsosFLEETp3AP, SEsosFLEETp3PA, sosFLEETp3, sosGREp3
	# sub-023: almost all runs
	# sub-024: SEsosFLEETp3AP, sosFLEETp3AP, sosGREp3
	# sub-025: almost all runs
	# sub-026: SEsosFLEETp3AP
	# sub-027: almost all runs
	# sub-028: almost all runs
	# sub-030: almost all runs
	# sub-031: sense1FLEETp3
	# sub-032: sense1FLEETp3, sosFLEETp3, sosFLEETp3AP
	# sub-033: sense1FLEETp3, SEsosFLEETp3AP, SEsosFLEETp3PA, sosFLEETp3, sosFLEETp3AP
	# sub-034: almost all runs
	# sub-035: almost all runs
	# ....
	#
	# fsleyes bold_denoised_mean.nii.gz -dr 0 15000 bold_denoised_mean_seg.nii.gz -cm red -a 20
	# fsleyes bold_denoised_mean.nii.gz -dr 0 1500 bold_denoised_mean_seg.nii.gz -cm red -a 20
	#
	# save corrected version under same name
	#
	# if necessary, run the tsnr extraction above again with the manually corrected segmentations
	# (just comment out the sct_deepseg command)
	# -------------------------------------------------------------------------------------------
	
	# -------------------------------------------------------------------------------------------
	# 3. Do motion-correction
	# -------------------------------------------------------------------------------------------
	if [ $moco = 1 ]; then
		for dir in "$out_dir"/*/; do
			cd "$dir"
			
			# Create mask for motion correction using the segmentation
			sct_create_mask -i bold_denoised_mean.nii.gz -p centerline,bold_denoised_mean_seg.nii.gz -size 41mm -o "bold_mean_mask.nii.gz"

			#----- Perform two-step motion correction -----

 			# merge the mean with the time series
			fslmerge -t mean_merged bold_denoised_mean bold_denoised
			
			# MOCO 1 with mean as target volume
 			sct_fmri_moco -i mean_merged.nii.gz -m bold_mean_mask.nii.gz -g 1 -x spline -param sampling=None,iterAvg=0,poly=2,iter=20
 			
 			# remove the mean volume again from the time series
			sct_image -i mean_merged_moco.nii.gz -remove-vol 0 -o first_MOCO.nii.gz
			
 			# create a mean of that motion corrected series
			fslmaths first_MOCO.nii.gz -Tmean first_MOCO_MEAN.nii.gz
			
			# now use this one as the first (target) volume
			fslmerge -t MOCO_mean_merged first_MOCO_MEAN bold_denoised
			
			# delete everything else from first moco
			rm -f moco_params* *_moco.nii.gz
	         	
			# ----------step 2---------
			# MOCO 2 with improved mean as target volume
			sct_fmri_moco -i MOCO_mean_merged.nii.gz -m bold_mean_mask.nii.gz -g 1 -x spline -param sampling=None,iterAvg=0,poly=2,iter=20
			 
			# remove the mean volume again from the time series
			sct_image -i MOCO_mean_merged_moco.nii.gz -remove-vol 0 -o MOCO.nii.gz
			
			# same for the moco params
			sct_image -i moco_params_x.nii.gz -remove-vol 0 -o moco_params_x.nii.gz
			sct_image -i moco_params_y.nii.gz -remove-vol 0 -o moco_params_y.nii.gz
				
			# and create a mean of that motion corrected series
			fslmaths MOCO.nii.gz -Tmean MOCO_MEAN.nii.gz
				
			# delete merged time series
			rm -f *_merged.nii.gz *_moco.nii.gz
				
			# Make a first segmentation of MOCO time series
			sct_deepseg_sc -i MOCO_MEAN.nii.gz -c t2s -o MOCO_MEAN_deepseg.nii.gz

		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# Correct the segmentation manually, save under filename MOCO_MEAN_deepseg_manual.nii.gz
	#
	# fsleyes MOCO_MEAN.nii.gz -dr 0 15000 MOCO_MEAN_deepseg.nii.gz -cm red -a 20
	# fsleyes MOCO_MEAN.nii.gz -dr 0 1500 MOCO_MEAN_deepseg.nii.gz -cm red -a 20
	# -------------------------------------------------------------------------------------------
	# -------------------------------------------------------------------------------------------
	# 4. Small steps after motion-correction & creation of CSF mask
	# manually corrected segmentation is needed!
	# -------------------------------------------------------------------------------------------
	if [ $csf = 1 ]; then
		for dir in "$out_dir"/*/; do
			cd "$dir"
			
			echo "Calculate tSNR of MOCO time series"
			
			# Create temporal mean and std to calculate the tsnr
			fslmaths MOCO.nii.gz -Tstd MOCO_sd.nii.gz
			fslmaths MOCO_MEAN.nii.gz -div MOCO_sd.nii.gz MOCO_tsnr.nii.gz
			# compare to the previous moco step
			if [ $run = 1 ]; then
				fslmaths first_MOCO.nii.gz -Tstd first_MOCO_sd.nii.gz
				fslmaths first_MOCO_MEAN.nii.gz -div first_MOCO_sd.nii.gz first_MOCO_tsnr.nii.gz
				fslmeants -i first_MOCO_tsnr.nii.gz -m MOCO_MEAN_deepseg_manual.nii.gz -o first_moco_tsnr_cord_manual.txt
			fi
			
			# Extract tsnr values in MOCO, original bold and denoised bold using the new mask
			fslmeants -i MOCO_tsnr.nii.gz -m MOCO_MEAN_deepseg_manual.nii.gz -o moco_tsnr_cord_manual.txt
			fslmeants -i bold_tsnr.nii.gz -m MOCO_MEAN_deepseg_manual.nii.gz -o bold_tsnr_cord_manual.txt
			fslmeants -i bold_denoised_tsnr.nii.gz -m MOCO_MEAN_deepseg_manual.nii.gz -o bold_denoised_tsnr_cord_manual.txt
			
			echo "Calculate spatial smoothness"
			# Extract values of spatial smoothness in original bold and denoised bold
			3dFWHMx -acf bold_acf_detrend_cord -input bold.nii.gz  -detrend -mask MOCO_MEAN_deepseg_manual.nii.gz  > bold_acf_detrend_cord_stdout.txt
			3dFWHMx -acf bold_denoised_acf_detrend_cord -input bold_denoised.nii.gz  -detrend -mask MOCO_MEAN_deepseg_manual.nii.gz  > bold_denoised_acf_detrend_cord_stdout.txt
			
		done
	fi
	
done

