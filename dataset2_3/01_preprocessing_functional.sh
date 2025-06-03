#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
preparation=1
thermnoise=1
quickseg=1
# manual correction
moco=0
mean_seg=0
# manual correction
prepare_zshift=0
# manual correction
do_zshift=0
moco_new=0
# manual correction
csf=0
motion=0
esPCA=0
# then switch to python, run: esPCA.py, motion_regr.py, moco_outlier.py, prepare_pnm.py

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder
code_dir=/data/pt_02661_raw/Heatpain/code  # code folder with subfolder helper_functions
helper_dir=$code_dir/helper_functions/

# needed for denoising:
mppca_dir=/data/u_uhorn_software/mppca_denoise-master/  # MPPCA toolbox downloaded from https://github.com/NYU-DiffusionMRI/mppca_denoise

# general processing software
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located
ants_dir=/data/u_uhorn_software/install/  # where your ANTS installation is located (for zshift only)
afni_dir=/afs/cbs.mpg.de/software/.afni/currentversion.debian-bullseye-amd64/debian-bullseye-amd64/  # where your AFNI installation is (in step csf for 3dFWHMx command only)

# needed for zshift:
NLMUpsample_dir=/data/u_uhorn_software/upsampling/demo2/demo/NLMUpsample/  # Upsampling tool from https://sites.google.com/site/pierrickcoupe/softwares/super-resolution/monomodal
nifti_tools_dir=/data/u_uhorn_software/nifti_tools/  # Nifti tools for matlab from https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

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

for subject in {1..41}; do
	printf -v sub "%02d" $subject
	echo "subject: " $sub
	
	if [ $subject -gt 16 ]; then
		data_dir=$project_dir/raw_data/dataset3/sub-sspr${sub}/func
	else
		data_dir=$project_dir/raw_data/dataset2/sub-sspr${sub}/func
	fi
	out_dir=$project_dir/derivatives/sub-sspr${sub}/func
	
	mkdir -p $out_dir
	
	# when some runs were not recorded
	if [[ $subject == 25 || $subject == 27 || $subject == 29 || $subject == 6 ]]; then
		num_runs=2
	elif [[ $subject == 31 ]]; then
		num_runs=3
	else
		num_runs=4
	fi
	
	# -------------------------------------------------------------------------------------------
	# 0. As preparation sort the data you need into run folders and give easier names
	# -------------------------------------------------------------------------------------------
	if [ $preparation = 1 ]; then
		cd $out_dir
		echo "Sorting data"
		for run in  $(seq 1 $num_runs); do
			mkdir -p "run-${run}"
			if [ $subject -gt 16 ]; then
				filename=sub-sspr${sub}_task-longpain_run-${run}_part-mag_bold.nii.gz
			else
				filename=sub-sspr${sub}_task-longpain_run-${run}_bold.nii.gz
			fi
			cp "${data_dir}/${filename}" "${out_dir}/run-${run}/bold.nii.gz"
		done
	fi
	

	# -------------------------------------------------------------------------------------------
	# 1. Thermal denoising using MPPCA
	# -------------------------------------------------------------------------------------------
	if [ $thermnoise = 1 ]; then
		for run in  $(seq 1 $num_runs); do
			cd ${out_dir}/run-${run}
			file="bold.nii.gz"
			fname=$(basename "$file" | cut -d. -f1)
			echo "Running thermal denoising run-${run}"
			echo $file
			input=$file
			folderout=${out_dir}/run-${run}/
			
			# shortcut to run matlab function from command line:
			matlab -nosplash -nodesktop -r "addpath('${helper_dir}'); \
						        try; thermalNoiseRemoval('${mppca_dir}', '${fsl_dir}', '${input}', '${fname}', '${folderout}'); catch ME; disp(ME); end; quit;";
			
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 2. Do a quick first segmentation of the cord
	# -------------------------------------------------------------------------------------------
	if [ $quickseg = 1 ]; then
		for run in  $(seq 1 $num_runs); do
			cd ${out_dir}/run-${run}
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
	# manual correction necessary in sub-sspr03 run-3; sub-sspr07 all runs; sub-sspr12 run-4;
	# sub-sspr18 all runs; sub-sspr29 all runs; sub-sspr30 all runs; sub-sspr32 run-3 and 4
	# sub-sspr34 run-1 and 2
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
		for run in $(seq 1 $num_runs); do
			cd ${out_dir}/run-${run}

			# Create mask for motion correction using the segmentation
			sct_create_mask -i bold_denoised_mean.nii.gz -p centerline,bold_denoised_mean_seg.nii.gz -size 41mm -o "bold_mean_mask.nii.gz"

			#----- Perform two-step motion correction -----
			# for the first run do MOCO on the mean and then concatenate the improved mean with the time series
         		if [ $run == 1 ]; then
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


             		# for the rest of the runs concatenate the mean of the 
	             	# motion-corrected data of the first run and original time series
         		else
				fslmerge -t MOCO_mean_merged ${out_dir}/run-1/first_MOCO_MEAN bold_denoised
	         	fi
	         	
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

		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 4. Create an overall mean and make a segmentation
	# -------------------------------------------------------------------------------------------
	if [ $mean_seg = 1 ]; then
		moco_means=()
		for run in $(seq 1 $num_runs); do
			func_dir=${out_dir}/run-${run}
			cd $func_dir
			moco_means+="$func_dir/MOCO_MEAN.nii.gz "
		done
		seclevel_dir=$out_dir/second_level
		mkdir $seclevel_dir
		cd $seclevel_dir
		fslmerge -t MOCO_MEAN_merge_all_runs.nii.gz ${moco_means[@]}
		fslmaths MOCO_MEAN_merge_all_runs.nii.gz -Tmean MOCO_MEAN_all_runs.nii.gz
		sct_deepseg_sc -i MOCO_MEAN_all_runs.nii.gz -c t2s -o MOCO_MEAN_deepseg.nii.gz

	fi
	
	# -------------------------------------------------------------------------------------------
	# Correct the segmentation manually, save under filename MOCO_MEAN_deepseg_manual.nii.gz
	#
	# fsleyes MOCO_MEAN.nii.gz -dr 0 15000 MOCO_MEAN_deepseg.nii.gz -cm red -a 20
	# fsleyes MOCO_MEAN.nii.gz -dr 0 1500 MOCO_MEAN_deepseg.nii.gz -cm red -a 20
	# -------------------------------------------------------------------------------------------
	
	# -------------------------------------------------------------------------------------------
	# 5. Prepare zshift by using a cord mask and moving it onto the discs
	# -------------------------------------------------------------------------------------------
	if [ $prepare_zshift = 1 ]; then
		seclevel_dir=$out_dir/second_level
		cd ${out_dir}/run-1
		echo "Dilating cord mask"
		sct_maths -i ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -dilate 1  -o MOCO_MEAN_deepseg_manual_dilated_1.nii.gz
		echo "Shift mask towards discs"
		matlab -nosplash -nodesktop -r "addpath('${helper_dir}');shifting_mask_to_discs('${nifti_tools_dir}', '${sub}', '${out_dir}/run-1/MOCO_MEAN_deepseg_manual_dilated_1.nii.gz');quit;";
		if [ $subject -gt 16 ]; then
			echo fsleyes ${out_dir}/run-1/MOCO_MEAN.nii.gz -dr 200 900 ${out_dir}/run-1/MOCO_MEAN_deepseg_manual_dilated_1_shifted.nii.gz -cm red -a 30
		else
			echo fsleyes ${out_dir}/run-1/MOCO_MEAN.nii.gz -dr 3000 10000 ${out_dir}/run-1/MOCO_MEAN_deepseg_manual_dilated_1_shifted.nii.gz -cm red -a 30
		fi
	fi
	
	# -------------------------------------------------------------------------------------------
	# afterwards check how well it matches the discs, 
	# maybe change value in matlab script shifting_mask_to_discs
	# -------------------------------------------------------------------------------------------
	
	# -------------------------------------------------------------------------------------------
	# 6. Calculate zshift to see whether a correction is necessary
	# -------------------------------------------------------------------------------------------
	if [ $do_zshift = 1 ]; then
		for run in $(seq 1 $num_runs); do
			cd ${out_dir}/run-${run}
			
			# factor for upsampling
			scaling=4
			if [ $run == 1 ]; then
				echo "Upsampling and binarizing mask"
				matlab -nosplash -nodesktop -r "addpath('${helper_dir}');upsampling('${NLMUpsample_dir}', '${nifti_tools_dir}', '${out_dir}/run-1/MOCO_MEAN_deepseg_manual_dilated_1_shifted.nii.gz', ${scaling});quit;";
				fslmaths MOCO_MEAN_deepseg_manual_dilated_1_shifted_upsampled.nii.gz -bin disc_mask.nii.gz
			fi
			
			echo "Upsampling MOCO MEAN image"
			matlab -nosplash -nodesktop -r "addpath('${helper_dir}');upsampling('${NLMUpsample_dir}', '${nifti_tools_dir}', '${out_dir}/run-${run}/MOCO_MEAN.nii.gz', ${scaling});quit;";
			
			echo "Restricting image to specific range"
			if [ $subject -gt 16 ]; then
				if [ $subject = 38 ]; then
					fslmaths MOCO_MEAN_upsampled.nii.gz -thr 100 -uthr 900 MOCO_MEAN_upsampled_range.nii.gz
				else
					fslmaths MOCO_MEAN_upsampled.nii.gz -thr 200 -uthr 900 MOCO_MEAN_upsampled_range.nii.gz
				fi
			else
				fslmaths MOCO_MEAN_upsampled.nii.gz -thr 3000 -uthr 10000 MOCO_MEAN_upsampled_range.nii.gz
			fi
			
			echo "Mask the image with disc mask"
    			mask=${out_dir}/run-1/disc_mask.nii.gz
			fslmaths MOCO_MEAN_upsampled_range.nii.gz -mas ${mask} discs.nii.gz
			
			echo "Dilating the discs in x and y direction"
			sct_maths -i discs.nii.gz -dilate 7 -shape disk -dim 2 -o dilated_discs.nii.gz
			
			if [ $run == 1 ]; then
				continue
			else
				echo "Aligning this run with the first one based on the discs"
	    			first=${out_dir}/run-1/dilated_discs.nii.gz
	    			this_run=${out_dir}/run-${run}/dilated_discs.nii.gz
	    			out_name=disc_matching
	    			# register the discs only allowing z-transformation
	    			isct_antsRegistration \
	    			--dimensionality 3 \
	    			--metric CC[${first},${this_run},1,4] \
	    			--transform Translation[1] \
	    			--restrict-deformation 0x0x1 \
	    			--convergence [100,1e-3,10] \
	    			--smoothing-sigmas 2 \
	    			--shrink-factors 1 \
	    			--output ${out_name}
	    			
	    			# and warp this runs' discs by this amount
	    			isct_antsApplyTransforms \
	    			--dimensionality 3 \
	    			--input ${this_run} \
	    			--reference-image ${first} \
	    			--output ${out_name}.nii.gz \
	    			--interpolation BSpline[3] \
	    			--transform ${out_name}0GenericAffine.mat
	    			
	    			# make it human-readable
	    			ConvertTransformFile 3 ${out_name}0GenericAffine.mat ${out_name}0GenericAffine.txt
	    			
	    			# read the file and search for the parameters, the last number is the one you need
	    			parameters_line=$(grep -n "Parameters: " ${out_name}0GenericAffine.txt)
	    			line_number=$(echo $parameters_line | cut -d':' -f1)
	    			numbers=$(sed -n "${line_number}s/Parameters: //p" ${out_name}0GenericAffine.txt)
	    			last_number=$(echo $numbers | awk '{print $NF}')
	    			
	    			# copy the file and replace the number with the correctly scaled number
	    			cp ${out_name}0GenericAffine.txt ${out_name}0GenericAffine_downsampled.txt
	    			zshift=$(echo "${last_number}/${scaling}" | bc -l)
	    			echo "The zshift is: $zshift"
	    			sed -i "s/${last_number}/${zshift}/g" ${out_name}0GenericAffine_downsampled.txt
	    			
	    			# create another transform file for the downsampled data
	    			ConvertTransformFile 3 ${out_name}0GenericAffine_downsampled.txt ${out_name}0GenericAffine_downsampled.mat --convertToAffineType
	    			
	    			
	    		fi
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 7. Do another motion-correction for those runs in which the zshift is strong
	# -------------------------------------------------------------------------------------------
	if [ $moco_new = 1 ]; then
		for run in $(seq 2 $num_runs); do
			cd ${out_dir}/run-${run}
			
			# read downsampled zshift
			parameters_line=$(grep -n "Parameters: " disc_matching0GenericAffine_downsampled.txt)
	    		line_number=$(echo $parameters_line | cut -d':' -f1)
	    		numbers=$(sed -n "${line_number}s/Parameters: //p" disc_matching0GenericAffine_downsampled.txt)
	    		zshift=$(echo $numbers | awk '{print $NF}')
			
			if (( $(echo "$zshift >= 0.75" | bc -l) )) || (( $(echo "$zshift <= -0.75" | bc -l) ));then
				echo "Need to shift subject $sub run $run"
				mkdir -p zshift
				
	    			antsApplyTransforms \
	    			--dimensionality 3 \
	    			--input-image-type 3 \
	    			--input bold_denoised.nii.gz \
	    			--reference-image ${out_dir}/run-1/MOCO_MEAN.nii.gz \
	    			--output bold_denoised_zshifted.nii.gz \
	    			--interpolation BSpline[3] \
	    			--transform disc_matching0GenericAffine_downsampled.mat \
	    			--verbose 1
	    			
	    			# in case the shift is larger than half the voxel thickness (= 3 mm)
	    			# then the most upper or lower slice will be empty --> fill it with the next slice content
	    			if (( $(echo "$zshift >= 1.5" | bc -l) )) || (( $(echo "$zshift <= -1.5" | bc -l) ));then
	    				echo "For this large shift we need to pad the resulting file"
	    				# negative shift = upwards = slice 0 missing
	    				if (( $(echo "$zshift <= -1.5" | bc -l) ));then
	    					fslroi bold_denoised_zshifted.nii.gz second_slice.nii.gz 0 -1 0 -1 1 1 0 -1
	    					fslroi bold_denoised_zshifted.nii.gz second_to_end_slice.nii.gz 0 -1 0 -1 1 -1 0 -1
	    					fslmerge -z bold_denoised_zshifted_padded.nii.gz second_slice.nii.gz second_to_end_slice.nii.gz
	    					fslcpgeom bold_denoised_zshifted.nii.gz bold_denoised_zshifted_padded.nii.gz
	    				# positive shift = downwards = max slice missing
	    				else
	    					dim3=$(fslinfo bold_denoised_zshifted.nii.gz | grep -m 1 dim3 | awk '{print $2}')
	    					printf -v slicesecondmax "%02d" $(($dim3-2))
	    					fslroi bold_denoised_zshifted.nii.gz second_slice.nii.gz 0 -1 0 -1 ${slicesecondmax} 1 0 -1
	    					fslroi bold_denoised_zshifted.nii.gz second_to_end_slice.nii.gz 0 -1 0 -1 0 $((slicesecondmax+1)) 0 -1
	    					fslmerge -z bold_denoised_zshifted_padded.nii.gz second_to_end_slice.nii.gz second_slice.nii.gz
	    				fi
	    				rm -f second*.nii.gz
	    				cp bold_denoised_zshifted_padded.nii.gz zshift/bold_denoised.nii.gz
	    			else
	    				cp bold_denoised_zshifted.nii.gz zshift/bold_denoised.nii.gz
	    			fi
	    			
	    			# take the old segmentation as it is not important how well it matches (just for the mask in moco)
	    			cp bold_denoised_mean_seg.nii.gz zshift/bold_denoised_mean_seg.nii.gz
	    			cd zshift
	    			
	    			fslmaths bold_denoised.nii.gz -Tmean bold_denoised_mean.nii.gz
	    			
	    			# Create mask for motion correction using the segmentation
				sct_create_mask -i bold_denoised_mean.nii.gz -p centerline,bold_denoised_mean_seg.nii.gz -size 41mm -o "bold_mean_mask.nii.gz"

				#----- Perform two-step motion correction -----
				# for the first run we have the moco mean already
		     		# for the rest of the runs concatenate the mean of the 
			     	# motion-corrected data of the first run and original time series

				fslmerge -t MOCO_mean_merged ${out_dir}/run-1/first_MOCO_MEAN bold_denoised
			 	
				# ----------step 2---------
				# MOCO 2 with improved mean as target volume
				sct_fmri_moco -i MOCO_mean_merged.nii.gz -m bold_mean_mask.nii.gz -g 1 -x spline -param sampling=None,iterAvg=0,poly=2,iter=20
				 
				# remove the mean volume again from the time series
				sct_image -i MOCO_mean_merged_moco.nii.gz -remove-vol 0 -o MOCO_zshifted.nii.gz
				
				# same for the moco params
				sct_image -i moco_params_x.nii.gz -remove-vol 0 -o moco_params_x.nii.gz
				sct_image -i moco_params_y.nii.gz -remove-vol 0 -o moco_params_y.nii.gz
					
				# and create a mean of that motion corrected series
				fslmaths MOCO_zshifted.nii.gz -Tmean MOCO_zshifted_MEAN.nii.gz
					
				# delete merged time series
				rm -f *_merged.nii.gz *_moco.nii.gz
				rm -f bold_denoised.nii.gz

			fi
		done
	fi
	# -------------------------------------------------------------------------------------------
	# Check whether the zshifted version looks better, if yes replace the MOCO and 
	# MOCO_MEAN and the moco parameter files with their z-shifted versions
	# run step 4 (mean_seg) again to create a new overall mean and segmentation
	# Correct the segmentation again manually, save under filename MOCO_MEAN_deepseg_manual.nii.gz
	# -------------------------------------------------------------------------------------------
	# -------------------------------------------------------------------------------------------
	# 8. Small steps after motion-correction & creation of CSF mask
	# manually corrected segmentation is needed!
	# -------------------------------------------------------------------------------------------
	if [ $csf = 1 ]; then
		seclevel_dir=$out_dir/second_level
		for run in $(seq 1 $num_runs); do
			cd ${out_dir}/run-${run}
			
			echo "Calculate tSNR of MOCO time series"
			
			# Create temporal mean and std to calculate the tsnr
			fslmaths MOCO.nii.gz -Tstd MOCO_sd.nii.gz
			fslmaths MOCO_MEAN.nii.gz -div MOCO_sd.nii.gz MOCO_tsnr.nii.gz
			 compare to the previous moco step
			if [ $run = 1 ]; then
				fslmaths first_MOCO.nii.gz -Tstd first_MOCO_sd.nii.gz
				fslmaths first_MOCO_MEAN.nii.gz -div first_MOCO_sd.nii.gz first_MOCO_tsnr.nii.gz
				fslmeants -i first_MOCO_tsnr.nii.gz -m ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -o first_moco_tsnr_cord_manual.txt
			fi
			
			# Extract tsnr values in MOCO, original bold and denoised bold using the new mask
			fslmeants -i MOCO_tsnr.nii.gz -m ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -o moco_tsnr_cord_manual.txt
			fslmeants -i bold_tsnr.nii.gz -m ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -o bold_tsnr_cord_manual.txt
			fslmeants -i bold_denoised_tsnr.nii.gz -m ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -o bold_denoised_tsnr_cord_manual.txt
			
			 echo "Calculate spatial smoothness"
			# Extract values of spatial smoothness in original bold and denoised bold
			3dFWHMx -acf bold_acf_detrend_cord -input bold.nii.gz  -detrend -mask ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz  > bold_acf_detrend_cord_stdout.txt
			3dFWHMx -acf bold_denoised_acf_detrend_cord -input bold_denoised.nii.gz  -detrend -mask ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz  > bold_denoised_acf_detrend_cord_stdout.txt

			# Create CSF + cord mask by creating a CSF segmentation, filling the holes and dilating it by 3
			sct_propseg -i MOCO_MEAN.nii.gz -c t2s -CSF
			fslmaths MOCO_MEAN_CSF_seg.nii.gz -fillh MOCO_MEAN_CSF_seg_filled.nii.gz
			sct_maths -i MOCO_MEAN_CSF_seg_filled.nii.gz -dilate 3 -dim 2 -shape disk -o MOCO_MEAN_CSF_seg_filled_dil3.nii.gz 
			fslmaths MOCO_MEAN_CSF_seg_filled_dil3.nii.gz -add ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -bin CSF_and_cord.nii.gz
			fslmaths MOCO_MEAN_CSF_seg_filled_dil3.nii.gz -sub ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz -bin CSF_only.nii.gz
			rm -f MOCO_MEAN_CSF_seg_filled_dil3.nii.gz MOCO_MEAN_CSF_seg_filled.nii.gz
			echo "CSF mask for run-${run} created"
			
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 9. Calculate parameters of abnormally high motion
	# -------------------------------------------------------------------------------------------
	if [ $motion = 1 ]; then
		export PATH="$code_dir:$PATH"
		export CODE_DIR=$code_dir

    		for run in $(seq 1 $num_runs); do
			func_dir=$out_dir/run-${run}
			cd $func_dir
        		
        		moco=MOCO.nii.gz
        		mask=bold_mean_mask.nii.gz
        		echo "Calculating DVARS to determine motion in run-${run}"
        		fsl_motion_outliers -i $moco -m $mask -o DVARS -s DVARS.txt --dvars --nomoco
        		echo "Calculating REFRMS to determine motion in run-${run}"
        		fsl_motion_outliers -i $moco -m $mask -o REFRMS -s REFRMS.txt --refrms --nomoco
        		echo "Calculating a modified version of REFRMS to determine motion in run-${run}"
        		$code_dir/helper_functions/fsl_motion_outliers_edited -i $moco -m $mask -o REFRMS_edited -s REFRMS_edited.txt --refrms --nomoco

            done
        fi


	# -------------------------------------------------------------------------------------------
	# 10. Prepare the extraspinal PCA by creating the mask and extracting the time series
	# -------------------------------------------------------------------------------------------
	if [ $esPCA = 1 ]; then
		for run in $(seq 1 $num_runs); do
			cd ${out_dir}/run-${run}
			
			# Create extra-spinal mask by dilating the CSF+cord mask and then subtracting it
			sct_maths -i CSF_and_cord.nii.gz -dilate 24 -dim 2 -shape disk -o tmp.nii.gz
			fslmaths tmp.nii.gz -sub CSF_and_cord.nii.gz extraspinal_mask.nii.gz
			fslmeants -i MOCO.nii.gz -m extraspinal_mask.nii.gz -o ts_for_espca.txt --showall
			rm -f tmp.nii.gz
			
		done
	fi
done

