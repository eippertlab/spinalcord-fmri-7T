#!/bin/bash

# This script requires the usage of ANTs N4BiasFieldCorrection function. We used ANTs version 2.3.1
# Additionally, you will need the helper function Denoising_T1.m provided by us, and Pierrick Coupés Denoising package MRIdenoisingPackage_r01.
# Get the denoising package under https://sites.google.com/site/pierrickcoupe/softwares/denoising/mri-denoising/mri-denoising-software

# --------------------------------
# Which steps do you want to run?
# --------------------------------
preparation=1
first_seg=1
label_vert=1
# visual inspection
# run helper function crop_anat.py
crop=0
bias_denoise=0
proper_seg=0
# visual inspection
proper_label_vert=0
# visual inspection
anat2template=0
merge=0
# visual inspection
create_averages=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder where raw and output data is located
code_dir=/data/pt_02661_raw/Heatpain/code  # code folder with subfolder helper_functions
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located
ants_dir=/afs/cbs.mpg.de/software/ants/2.3.1/ubuntu-bionic-amd64/antsbin/bin # where ANTS bin is located
package_dir=/data/u_uhorn_software/MRIDenoisingPackage_r01_pcode # where the MRI denoising package is located (see description above)
spm_path=/data/u_uhorn_software/spm  # where SPM is installed (the denoising tool uses it)

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

	if [ $subject -gt 16 ]; then
		data_dir=$project_dir/raw_data/dataset3/sub-sspr${sub}/anat
		intensity=800
	else
		data_dir=$project_dir/raw_data/dataset2/sub-sspr${sub}/anat
		intensity=1200
	fi
	out_dir=$project_dir/derivatives/sub-sspr${sub}/anat
	
	mkdir -p $out_dir
	
	# -------------------------------------------------------------------------------------------
	# 0. As preparation search for the anatomy scans and give them easier names
	# -------------------------------------------------------------------------------------------
	if [ $preparation -eq 1 ]; then
		cd $out_dir
		echo "Getting the data"
		
		if [ $subject -gt 16 ]; then
			filename=sub-sspr${sub}_acq-ND_inv-2_MP2RAGE.nii.gz
		else
			filename=sub-sspr${sub}_T1w.nii.gz
		fi
		cp "${data_dir}/${filename}" "${out_dir}/anat.nii.gz"
	fi
	
	# -------------------------------------------------------------------------------------------
	# 1. Do a first rough segmentation
	# -------------------------------------------------------------------------------------------
	if [ $first_seg -eq 1 ]; then
		cd $out_dir
		echo "Doing a first segmentation"
	    
		sct_deepseg_sc -i anat.nii.gz -c t1 -o anat_initial_deepseg.nii.gz
		
		echo "check that segmentation worked, copy this into terminal:"
		echo "fsleyes ${out_dir}/anat.nii.gz -dr 0 ${intensity} ${out_dir}/anat_initial_deepseg.nii.gz --cmap red"
	fi
	
	# -------------------------------------------------------------------------------------------
	# For sub-sspr19 there was some issue with the vertebral labeling that could be fixed by
	# drawing a bit nicer segmentation,
	# I also fixed the upper parts of sub-sspr33
	# -------------------------------------------------------------------------------------------
	
	# -------------------------------------------------------------------------------------------
	# 2. Do a first vertebrae labeling to know where to crop
	# -------------------------------------------------------------------------------------------
	if [ $label_vert -eq 1 ]; then
		cd $out_dir
		echo "Doing a first vertebrae labeling"
	    
		sct_label_vertebrae -i anat.nii.gz -s anat_initial_deepseg.nii.gz -c t1
		
        	# extract labels for registration
        	sct_label_utils -i anat_initial_deepseg_labeled_discs.nii.gz -keep 2,3,4,5,6,7,8,9 -o  disc_labels.nii.gz
		
		echo "check that labeling worked, copy this into terminal:"
		echo "fsleyes ${out_dir}/anat.nii.gz -dr 0 ${intensity} ${out_dir}/anat_initial_deepseg_labeled.nii.gz --cmap brain_colours_x_rain"
		
		# I already checked and noticed the labels are off for these subjects, correct manually
                if [[ ${subject} == 12 || ${subject} == 14 || ${subject} == 17 || ${subject} == 20 || ${subject} == 21 || ${subject} == 24 || ${subject} == 26 || ${subject} == 34 || ${subject} == 35 || ${subject} == 36 || ${subject} == 38 ]]; then
                	sct_label_utils -i anat.nii.gz -create-viewer 3 -o label_c2c3.nii.gz -msg "Click at the posterior tip of C2/C3 inter-vertebral disc"
                	sct_label_vertebrae -i anat.nii.gz -s anat_initial_deepseg.nii.gz -c t1 -initlabel label_c2c3.nii.gz
                	sct_label_utils -i anat_initial_deepseg_labeled_discs.nii.gz -keep 2,3,4,5,6,7,8,9 -o  disc_labels.nii.gz
                fi

	fi
	
	# -------------------------------------------------------------------------------------------
	# at this point run helper function crop_anat.py
	# -------------------------------------------------------------------------------------------
	
	# -------------------------------------------------------------------------------------------
	# 3. Crop image from C1 to just below T1 and in x and y to get rid of aliasing
	# -------------------------------------------------------------------------------------------
	if [ $crop -eq 1 ]; then
		cd $out_dir
		echo "Cropping the anatomical image"
	    	
	    	if [ -e "crop_command.txt" ];then
	    		source "crop_command.txt"
	    		echo fsleyes ${out_dir}/anat.nii.gz -dr 0 ${intensity} ${out_dir}/anat_crop.nii.gz -dr 0 ${intensity} -cm red -a 30
	    	else
	    		echo "Crop command does not exist. Please run python script crop_anat.py"
	    	fi

	fi
	
	# -------------------------------------------------------------------------------------------
	# 4. Bias correction and denoising of the cropped image
	# -------------------------------------------------------------------------------------------
	if [ $bias_denoise -eq 1 ]; then
		cd $out_dir

		echo "Starting bias correction"
	    	
	    	$ants_dir/N4BiasFieldCorrection \
	    	-i anat_crop.nii.gz \
	    	-o anat_crop_biascorr.nii.gz \
	    	-c [200x200x200x200,0.000001]
          	sleep 10
        	echo "Finished bias correction"
        	
        	echo "Starting denoising"
        	
        	# unzip file because toolbox can't deal with .gz files
        	gunzip -k anat_crop_biascorr.nii.gz
    		sleep 3
    		# Run MATLAB function with verbose output
    		matlab -nosplash -nodesktop -r "addpath(genpath('${code_dir}'));addpath('${spm_path}'); try; Denoising_T1('${package_dir}', 1, 2, 1, 1, 3, '_denoise', 0, 'anat_crop_biascorr.nii', '${out_dir}'); catch ME;disp(getReport(ME)); end; quit;"
		
    		# Wait for output file to be created
		while [ ! -f "${out_dir}/anat_crop_biascorr_denoise.nii" ]; do
			sleep 1
		done
		
		# zip resulting file and remove unzipped files
		gzip -f anat_crop_biascorr_denoise.nii
		rm -f anat_crop_biascorr.nii
		
		# wait before executing next call to Matlab; prevents Matlab from randomly crashing on startup
    		sleep 5
    		
    		echo "Finished denoising"
        	
	fi
	
	# -------------------------------------------------------------------------------------------
	# 5. Do a proper segmentation
	# -------------------------------------------------------------------------------------------
	if [ $proper_seg -eq 1 ]; then
		cd $out_dir
		echo "Doing a proper segmentation"
	    	
	    	# segment once (for some subjects this did not work out well, find the centerline manually there)
	    	if [[ ${subject} == 5 || ${subject} == 8 || ${subject} == 26 || ${subject} == 32 || ${subject} == 34 || ${subject} == 35 ]];then
	    		sct_get_centerline -i anat_crop_biascorr_denoise.nii.gz -c t1 -method viewer
	    	else
			sct_deepseg_sc -i anat_crop_biascorr_denoise.nii.gz -c t1 -o anat_deepseg.nii.gz
		fi
		
		# smooth along cord
		if [[ ${subject} == 5 || ${subject} == 8 || ${subject} == 26 || ${subject} == 32 || ${subject} == 34 || ${subject} == 35 ]];then
			sct_smooth_spinalcord -i anat_crop_biascorr_denoise.nii.gz \
			-s anat_crop_biascorr_denoise_centerline.nii.gz -smooth 0,0,6 -v 1 \
			-o anat_crop_biascorr_denoise_smooth.nii.gz

		else
			sct_smooth_spinalcord -i anat_crop_biascorr_denoise.nii.gz \
			-s anat_deepseg.nii.gz -smooth 0,0,6 -v 1 \
			-o  anat_crop_biascorr_denoise_smooth.nii.gz
		fi
        
        	# segment again and save under same name
		sct_deepseg_sc -i anat_crop_biascorr_denoise_smooth.nii.gz -c t1 -o anat_final_deepseg.nii.gz
		
		echo "check that segmentation worked, copy this into terminal:"
		echo "fsleyes ${out_dir}/anat_crop_biascorr_denoise.nii.gz -dr 0 ${intensity} ${out_dir}/anat_final_deepseg.nii.gz --cmap red"
	fi
	# -------------------------------------------------------------------------------------------
	# afterwards check quality of segmentation, manually correct and save under
        # anat_final_deepseg_manual.nii.gz
        # (even if you did not need to correct anything)
        # -------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------
	# 6. Do a proper vertebrae labeling
	# -------------------------------------------------------------------------------------------
	if [ $proper_label_vert -eq 1 ]; then
		cd $out_dir
		echo "Doing a proper labeling of the vertebrae"
		
		if [ -e "anat_final_deepseg_manual.nii.gz" ];then
			sct_label_vertebrae -i anat_crop_biascorr_denoise.nii.gz -s anat_final_deepseg_manual.nii.gz -c t1
			sct_label_utils -i anat_final_deepseg_manual_labeled_discs.nii.gz -keep 3,7 -o  disc_labels_cropped_3_7.nii.gz
        	else
            		echo "You did not manually correct the segmentations"
            	fi
		
		echo "check that the labeling worked, copy this into terminal:"
		echo fsleyes ${out_dir}/anat_crop_biascorr_denoise.nii.gz -dr 0 ${intensity} \
		${out_dir}/anat_final_deepseg_manual_labeled.nii.gz --cmap brain_colours_x_rain
		
		# here I checked and it did not find the correct labels, so correct them manually
		if [[ ${subject} == 2 || ${subject} == 12 || ${subject} == 14 || ${subject} == 20 || ${subject} == 21 || ${subject} == 24 || ${subject} == 26 || ${subject} == 34 || ${subject} == 35 || ${subject} == 36 || ${subject} == 38 ]];then
			sct_label_utils -i anat_crop_biascorr_denoise.nii.gz -create-viewer 3 -o label_c2c3.nii.gz
			sct_label_vertebrae -i anat_crop_biascorr_denoise.nii.gz -s anat_final_deepseg_manual.nii.gz -c t1 -initlabel label_c2c3.nii.gz
			sct_label_utils -i anat_final_deepseg_manual_labeled_discs.nii.gz -keep 3,7 -o disc_labels_cropped_3_7.nii.gz
		# here this did not even help, so we need to label all manually
		elif [[ ${subject} == 6 ]]; then
			sct_label_utils -i anat_crop_biascorr_denoise.nii.gz -create-viewer 2,3,4,5,6,7,8 \
			-o disc_labels_cropped.nii.gz -msg "Place labels at the posterior tip of each inter-vertebral disc. E.g. Label 3: C2/C3, Label 4: C3/C4, etc."
			sct_label_utils -i disc_labels_cropped.nii.gz  -keep 3,7 -o disc_labels_cropped_3_7.nii.gz
		fi
	fi
	
	# -------------------------------------------------------------------------------------------
	# 7. Register anatomy to template
	# -------------------------------------------------------------------------------------------
	if [ $anat2template -eq 1 ]; then
		cd $out_dir
		mkdir -p "seg_centermass_seg_bsplinesyn"
		
		echo "Registering the anatomy to the PAM50 template centermass seg-based bspline"
		sct_register_to_template -i anat_crop_biascorr_denoise.nii.gz \
		-s anat_final_deepseg_manual.nii.gz \
		-c t1 -ldisc disc_labels_cropped_3_7.nii.gz \
		-param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn \
		-ofolder ${out_dir}/seg_centermass_seg_bsplinesyn
		
		mkdir -p "seg_centermass_im_bsplinesyn"
		
		echo "Registering the anatomy to the PAM50 template centermass image-based bspline"
		sct_register_to_template -i anat_crop_biascorr_denoise.nii.gz \
		-s anat_final_deepseg_manual.nii.gz \
		-c t1 -ldisc disc_labels_cropped_3_7.nii.gz \
		-param step=1,type=seg,algo=centermass:step=2,type=im,algo=bsplinesyn,metric=MI \
		-ofolder ${out_dir}/seg_centermass_im_bsplinesyn
		
	fi

done

# -------------------------------------------------------------------------------------------
# 8. Create merged images of all anatomical images in template space
# -------------------------------------------------------------------------------------------
if [ $merge -eq 1 ]; then
	echo "Collating all template space anatomical images - seg-based"
	group_dir=$project_dir/derivatives/group_analysis
	mkdir -p ${group_dir}
	cd ${group_dir}
	sub_anat=()
	for subject in {1..41}; do
		printf -v sub "%02d" $subject
		sub_anat+="$project_dir/derivatives/sub-sspr${sub}/anat/seg_centermass_seg_bsplinesyn/anat2template.nii.gz "
	done
	fslmerge -t anat2template_segcentermass_segbsplinesyn_merged.nii.gz ${sub_anat[@]}
	
	echo "Collating all template space anatomical images - image-based"
	group_dir=$project_dir/derivatives/group_analysis
	mkdir -p ${group_dir}
	cd ${group_dir}
	sub_anat=()
	for subject in {1..41}; do
		printf -v sub "%02d" $subject
		sub_anat+="$project_dir/derivatives/sub-sspr${sub}/anat/seg_centermass_im_bsplinesyn/anat2template.nii.gz "
	done
	fslmerge -t anat2template_segcentermass_imbsplinesyn_merged.nii.gz ${sub_anat[@]}

fi

# -------------------------------------------------------------------------------------------
# after that we selected manually which subjects would be best warped with seg+seg-based and
# which would be best warped with seg+image-based normalization (see step below)
# -------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# 9. Create averages of all anatomical images in template space
# -------------------------------------------------------------------------------------------
if [ $create_averages -eq 1 ]; then
	# copy the best normalization into the main folder
	echo "Collecting optimal images - dataset 1"
	sub_anat=()
	for subject in {1..16}; do
		printf -v sub "%02d" $subject
		out_dir=$project_dir/derivatives/sub-sspr${sub}/anat
		if [[ ${subject} == 1 || ${subject} == 2 || ${subject} == 11 ]]; then
			norm_path=$out_dir/seg_centermass_im_bsplinesyn
		else
			norm_path=$out_dir/seg_centermass_seg_bsplinesyn
		fi
		cp -r ${norm_path}/* ${out_dir}/
		sub_anat+="${norm_path}/anat2template.nii.gz "
	done
	echo "Creating an average of all template space anatomical images in dataset 1"
	group_dir=$project_dir/derivatives/group_analysis
	mkdir -p ${group_dir}
	cd ${group_dir}
	fslmerge -t anat2template_optimal_merged_dataset1.nii.gz ${sub_anat[@]}
	fslmaths anat2template_optimal_merged_dataset1.nii.gz -Tmean anat2template_optimal_avg_dataset1.nii.gz
	
	# copy the best normalization into the main folder
	echo "Collecting optimal images - dataset 2"
	sub_anat=()
	for subject in {17..41}; do
		printf -v sub "%02d" $subject
		out_dir=$project_dir/derivatives/sub-sspr${sub}/anat
		if [[ ${subject} == 18 || ${subject} == 19 || ${subject} == 22 || ${subject} == 33 || ${subject} == 35 || ${subject} == 38 ]]; then
			norm_path=$out_dir/seg_centermass_im_bsplinesyn
		else
			norm_path=$out_dir/seg_centermass_seg_bsplinesyn
		fi
		cp -r ${norm_path}/* ${out_dir}/
		sub_anat+="${norm_path}/anat2template.nii.gz "
	done
	echo "Creating an average of all template space anatomical images in dataset 2"
	group_dir=$project_dir/derivatives/group_analysis
	mkdir -p ${group_dir}
	cd ${group_dir}
	fslmerge -t anat2template_optimal_merged_dataset2.nii.gz ${sub_anat[@]}
	fslmaths anat2template_optimal_merged_dataset2.nii.gz -Tmean anat2template_optimal_avg_dataset2.nii.gz
fi
