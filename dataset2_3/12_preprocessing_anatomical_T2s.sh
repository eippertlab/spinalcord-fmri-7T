#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
preparation=1
combine_images=0
segment=0
register_epi=0
warp_masks=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder where raw and output data is located
code_dir=/data/pt_02661_raw/Heatpain/code  # code folder with subfolder helper_functions
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located

# add FSL to the path
FSLDIR=${fsl_dir}
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# add SCT to the path
export PATH="${sct_dir}/bin:${PATH}"

for subject in {1..41}; do
	printf -v sub "%02d" $subject
	echo "subject: " $sub
	
	# 33 is missing
	if [ $subject -eq 33 ]; then
		continue
	fi
	
	if [ $subject -gt 16 ]; then
		data_dir=$project_dir/raw_data/dataset3/sub-sspr${sub}/anat
	else
		data_dir=$project_dir/raw_data/dataset2/sub-sspr${sub}/anat
	fi
	out_dir=$project_dir/derivatives/sub-sspr${sub}/anat_t2s
	
	mkdir -p $out_dir
	
	# -------------------------------------------------------------------------------------------
	# 0. As preparation search for the anatomy scans and give them easier names
	# -------------------------------------------------------------------------------------------
	if [ $preparation == 1 ]; then
		cd $data_dir
		echo "Getting the data"
		for file in *MEGRE.nii.gz;do
			if [[ "$file" == *"CoilSensitivity"* ]]; then
				echo "skipping files with name coil sensitivity"
				continue
			fi
			# remove subject name from beginning
			new_filename="${file/sub-sspr$sub/t2s}"
			# remove MEGRE from ending
			new_filename="${new_filename%_MEGRE.nii.gz}.nii.gz"
			# copy file to output folder
			cp $file ${out_dir}/$new_filename
		done
	fi

	# -------------------------------------------------------------------------------------------
	# 1. Combine all images
	# -------------------------------------------------------------------------------------------
	if [ $combine_images = 1 ]; then
		cd $out_dir
	    	
	    	# combine echos using root mean square
	    	echo "Combining the echos"
	    	sum_of_squares="sum_of_squares.nii.gz"
	    	for run in {1..3}; do
	    		fslmaths "t2s_run-${run}_echo-1.nii.gz" -mul 0 $sum_of_squares
	    		for echo in {1..5}; do
	    			fslmaths "t2s_run-${run}_echo-${echo}.nii.gz" -sqr -add $sum_of_squares $sum_of_squares
	    		done
	    		fslmaths $sum_of_squares -div 5 $sum_of_squares
	    		rms_image="rms_combined_run-${run}.nii.gz"
	    		fslmaths $sum_of_squares -sqrt $rms_image
	    		rm -f $sum_of_squares
	    	done
	    	
	    	echo "Build average of these"
	    	fslmerge -t rms_combined_merged.nii.gz rms_combined_run-1.nii.gz rms_combined_run-2.nii.gz rms_combined_run-3.nii.gz
	    	grand_avg=rms_combined_mean
	    	fslmaths rms_combined_merged.nii.gz -Tmean ${grand_avg}.nii.gz
	    	
	    	echo "Do motion-correction"
	    	sct_deepseg_sc -i ${grand_avg}.nii.gz -c t2s
	    	sct_create_mask -i ${grand_avg}.nii.gz -size 41mm -p centerline,${grand_avg}_seg.nii.gz
	    	sct_fmri_moco -i rms_combined_merged.nii.gz -m mask_${grand_avg}.nii.gz -g 1 -param sampling=None,iterAvg=0 -x nn
	    	mv rms_combined_merged_moco.nii.gz MEGRE_RMS_MOCO.nii.gz
	    	mv rms_combined_merged_moco_mean.nii.gz MEGRE_RMS_MOCO_MEAN.nii.gz
	    	rm -f moco_params* rms_combined_merged.nii.gz ${grand_avg}
    	fi
    	
    	# -------------------------------------------------------------------------------------------
	# 2. Do a segmentation of cord and gray matter
	# -------------------------------------------------------------------------------------------
    	if [ $segment = 1 ]; then
    		cd $out_dir
    		echo "Starting segmentation"
    		sct_deepseg_sc -i MEGRE_RMS_MOCO_MEAN.nii.gz -c t2s
    		sct_deepseg_gm -i MEGRE_RMS_MOCO_MEAN.nii.gz
    	fi
    	
    	# -------------------------------------------------------------------------------------------
    	# maybe check these segmentations also manually,
    	# e.g. in sub-sspr01 and sub-sspr31 there was a ghost that was also detected in the gm mask
    	# -------------------------------------------------------------------------------------------
    	
    	# -------------------------------------------------------------------------------------------
	# 3. Register MEGRE and EPI
	# -------------------------------------------------------------------------------------------
	if [ $register_epi = 1 ]; then
    		cd $out_dir
    		echo "Registering EPI and MEGRE images"
    		func_dir=$project_dir/derivatives/sub-sspr${sub}/func/second_level
		
		# register once using only translation to not change the overall form of the cord
		# when assessing how distorted the EPI is in comparison --> for Dice coefficient
    		sct_register_multimodal \
    		-i MEGRE_RMS_MOCO_MEAN.nii.gz \
    		-d $func_dir/MOCO_MEAN_all_runs.nii.gz \
    		-iseg MEGRE_RMS_MOCO_MEAN_seg_manual.nii.gz \
    		-dseg $func_dir/MOCO_MEAN_deepseg_manual.nii.gz \
    		-param step=1,type=seg,algo=translation \
    		-x nn \
    		-o MEGRE_RMS_MOCO_MEAN_translation_reg.nii.gz \
    		-owarp warp_MEGRE_RMS_MOCO_MEAN2MOCO_MEAN_all_runs_translation.nii.gz \
    		-owarpinv warp_MOCO_MEAN_all_runs2MEGRE_RMS_MOCO_MEAN_translation.nii.gz
    		
    		# register another time to really warp the gray matter segmentation into the EPI space
    		# here change the cords distortions levels!
    		sct_register_multimodal \
    		-i MEGRE_RMS_MOCO_MEAN.nii.gz \
    		-d $func_dir/MOCO_MEAN_all_runs.nii.gz \
    		-iseg MEGRE_RMS_MOCO_MEAN_seg_manual.nii.gz \
    		-dseg $func_dir/MOCO_MEAN_deepseg_manual.nii.gz \
    		-param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1 \
    		-o MEGRE_RMS_MOCO_MEAN_bspline_reg.nii.gz \
    		-owarp warp_MEGRE_RMS_MOCO_MEAN2MOCO_MEAN_all_runs_bspline.nii.gz \
    		-owarpinv warp_MOCO_MEAN_all_runs2MEGRE_RMS_MOCO_MEAN_bspline.nii.gz

    	fi
    	
    	# -------------------------------------------------------------------------------------------
	# 4. Warp masks from MEGRE onto EPI space
	# -------------------------------------------------------------------------------------------
	if [ $warp_masks = 1 ]; then
    		cd $out_dir
    		echo "Warping some masks into EPI space"
		
		func_dir=$project_dir/derivatives/sub-sspr${sub}/func/second_level
		
		# 1. translation warp for Dice coeff cord mask
		sct_apply_transfo \
		-i MEGRE_RMS_MOCO_MEAN_seg_manual.nii.gz \
		-d $func_dir/MOCO_MEAN_all_runs.nii.gz \
		-w warp_MEGRE_RMS_MOCO_MEAN2MOCO_MEAN_all_runs_translation.nii.gz \
		-x nn \
		-o MEGRE_RMS_MOCO_MEAN_seg_EPIspace_translation.nii.gz
		
		# 2. Bspline warp for gray matter mask
		sct_apply_transfo \
		-i MEGRE_RMS_MOCO_MEAN_gmseg.nii.gz \
		-d $func_dir/MOCO_MEAN_all_runs.nii.gz \
		-w warp*MEGRE*MOCO_MEAN_all_runs_bspline.nii.gz \
		-x nn \
		-o MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline.nii.gz
		
    	fi
        
    	
done
