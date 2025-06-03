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
normalize_func=0
tsnr2template=0
create_average_tsnr=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02616/data  # general project folder where raw and output data is located
code_dir=/data/pt_02616/code  # code folder with subfolder helper_functions
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

for subject in {1..62}; do
	printf -v sub "%03d" $subject
	echo "subject: " $sub

	data_dir=$project_dir/raw_data/sub-${sub}/anat
	intensity=1200

	out_dir=$project_dir/derivatives/sub-${sub}/anat
	
	mkdir -p $out_dir
	
	# -------------------------------------------------------------------------------------------
	# 0. As preparation search for the anatomy scans and give them easier names
	# -------------------------------------------------------------------------------------------
	if [ $preparation -eq 1 ]; then
		cd $out_dir
		echo "Getting the data"
		
		filename=sub-${sub}_T1w.nii.gz

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
	# a lot missing in sub-050 --> drew this again
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
                if [[ (${subject} = 57) ]]; then
                	sct_label_utils -i anat.nii.gz -create-viewer 3 -o label_c2c3.nii.gz -msg "Click at the posterior tip of C2/C3 inter-vertebral disc"
                	sct_label_vertebrae -i anat.nii.gz -s anat_initial_deepseg.nii.gz -c t1 -initlabel label_c2c3.nii.gz
                	sct_label_utils -i anat_initial_deepseg_labeled_discs.nii.gz -keep 2,3,4,5,6,7,8,9 -o  disc_labels.nii.gz
                # here this did not even help, so we need to label all manually
		elif [[ (${subject} = 29) || (${subject} = 37) ]]; then
			sct_label_utils -i anat.nii.gz -create-viewer 1,2,3,4,5,6,7,8,9,10 \
			-o disc_labels.nii.gz -msg "Place labels at the posterior tip of each inter-vertebral disc. E.g. Label 2: C1/C2, Label 3: C2/C3, Label 4: C3/C4, etc."
			sct_label_vertebrae -i anat.nii.gz -s anat_initial_deepseg.nii.gz -c t1 -discfile disc_labels.nii.gz
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
		sct_deepseg_sc -i anat_crop_biascorr_denoise.nii.gz -c t1 -o anat_deepseg.nii.gz
		
		# smooth along cord
		sct_smooth_spinalcord -i anat_crop_biascorr_denoise.nii.gz \
		-s anat_deepseg.nii.gz -smooth 0,0,6 -v 1 \
		-o  anat_crop_biascorr_denoise_smooth.nii.gz
        
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
		if [[ (${subject} -eq 55) || (${subject} -eq 57) ]]; then
                	sct_label_utils -i anat_crop_biascorr_denoise.nii.gz -create-viewer 3 -o label_c2c3.nii.gz
			sct_label_vertebrae -i anat_crop_biascorr_denoise.nii.gz -s anat_final_deepseg_manual.nii.gz -c t1 -initlabel label_c2c3.nii.gz
			sct_label_utils -i anat_final_deepseg_manual_labeled_discs.nii.gz -keep 3,7 -o disc_labels_cropped_3_7.nii.gz
                elif [[ (${subject} -eq 11) || (${subject} -eq 17) || (${subject} -eq 24) ]]; then
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
		-param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1 \
		-ofolder ${out_dir}/seg_centermass_seg_bsplinesyn_slicewise
		
		mkdir -p "seg_centermass_im_bsplinesyn"
		
		echo "Registering the anatomy to the PAM50 template centermass image-based bspline"
		sct_register_to_template -i anat_crop_biascorr_denoise.nii.gz \
		-s anat_final_deepseg_manual.nii.gz \
		-c t1 -ldisc disc_labels_cropped_3_7.nii.gz \
		-param step=1,type=seg,algo=centermass:step=2,type=im,algo=bsplinesyn,metric=MI,slicewise=1 \
		-ofolder ${out_dir}/seg_centermass_im_bsplinesyn_slicewise
		
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
	for subject in {1..62}; do
		printf -v sub "%03d" $subject
		sub_anat+="$project_dir/derivatives/sub-${sub}/anat/seg_centermass_seg_bsplinesyn_slicewise/anat2template.nii.gz "
	done
	fslmerge -t anat2template_segcentermass_segbsplinesyn_slicewise_merged.nii.gz ${sub_anat[@]}
	
	echo "Collating all template space anatomical images - image-based"
	group_dir=$project_dir/derivatives/group_analysis
	mkdir -p ${group_dir}
	cd ${group_dir}
	sub_anat=()
	for subject in {1..62}; do
		printf -v sub "%03d" $subject
		sub_anat+="$project_dir/derivatives/sub-${sub}/anat/seg_centermass_im_bsplinesyn_slicewise/anat2template.nii.gz "
	done
	fslmerge -t anat2template_segcentermass_imbsplinesyn_slicewise_merged.nii.gz ${sub_anat[@]}
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
	echo "Collecting optimal images"
	sub_anat=()
	for subject in {1..62}; do
		printf -v sub "%03d" $subject
		out_dir=$project_dir/derivatives/sub-${sub}/anat
		if [[ ${subject} == 1 || ${subject} == 2 || ${subject} == 3 || ${subject} == 4 || ${subject} == 6 || ${subject} == 7 || ${subject} == 10 || ${subject} == 11 || ${subject} == 13 || ${subject} == 14 || ${subject} == 16 || ${subject} == 17 || ${subject} == 18 || ${subject} == 19 || ${subject} == 20 || ${subject} == 21 || ${subject} == 25 || ${subject} == 29 || ${subject} == 34 || ${subject} == 37 || ${subject} == 40 || ${subject} == 43 || ${subject} == 44 || ${subject} == 46 || ${subject} == 47 || ${subject} == 48 || ${subject} == 49 || ${subject} == 50 || ${subject} == 51 || ${subject} == 53 || ${subject} == 58 || ${subject} == 59 || ${subject} == 60 || ${subject} == 61 ]]; then
			norm_path=$out_dir/seg_centermass_im_bsplinesyn
		else
			norm_path=$out_dir/seg_centermass_seg_bsplinesyn
		fi
		cp -r ${norm_path}/* ${out_dir}/
		sub_anat+="${norm_path}/anat2template.nii.gz "
	done
	echo "Creating an average of all template space anatomical images in dataset"
	group_dir=$project_dir/derivatives/group_analysis
	mkdir -p ${group_dir}
	cd ${group_dir}
	fslmerge -t anat2template_optimal_merged.nii.gz ${sub_anat[@]}
	fslmaths anat2template_optimal_merged.nii.gz -Tmean anat2template_optimal_avg.nii.gz
fi

# -------------------------------------------------------------------------------------------
# 10. Register functional image to template space
# -------------------------------------------------------------------------------------------
if [ $normalize_func -eq 1 ]; then
	for subject in {1..62}; do
		printf -v sub "%03d" $subject
		
		if [[ "$subject" == 29 || "$subject" == 60 ]]; then
			echo "Subject has no functional data"
			continue
		elif [[ "$subject" == 56 ]]; then
			echo "Subject has severe artifacts"
			continue
		fi
		anat_dir=$project_dir/derivatives/sub-${sub}/anat
		func_dir=$project_dir/derivatives/sub-${sub}/func/sosGREp3
		
		func_anat_dir=$project_dir/derivatives/sub-${sub}/func_anat/sosGREp3
		mkdir -p $func_anat_dir
	
		cd $func_anat_dir
		
		# register the moco mean to the template by using the warping fields from anat2template transformation
		# in the anat folder there is already the optimal of the two warping fields
		
		mkdir -p "optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn"
		
		sct_register_multimodal \
		-i  ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
		-d ${func_dir}/MOCO_MEAN.nii.gz \
		-iseg ${sct_dir}/data/PAM50/template/PAM50_cord.nii.gz \
		-dseg ${func_dir}/MOCO_MEAN_deepseg_manual.nii.gz \
		-param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,smooth=1,slicewise=1,iter=3 \
		-initwarp ${anat_dir}/warp_template2anat.nii.gz \
		-initwarpinv ${anat_dir}/warp_anat2template.nii.gz \
		-x spline -ofolder ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn
	done	
fi
	
# -------------------------------------------------------------------------------------------
# 10. Warp tSNR images to template space
# -------------------------------------------------------------------------------------------
if [ $tsnr2template -eq 1 ]; then
	for subject in {1..62}; do
		printf -v sub "%03d" $subject
		if [[ "$subject" == 29 || "$subject" == 60 ]]; then
			echo "Subject has no functional data"
			continue
		elif [[ "$subject" == 56 ]]; then
			echo "Subject has severe artifacts"
			continue
		fi	
		func_dir=${project_dir}/derivatives/sub-${sub}/func/sosGREp3
		func_anat_dir=$project_dir/derivatives/sub-${sub}/func_anat/sosGREp3
		
		echo "Warping tSNR image into template space"
		sct_apply_transfo -i ${func_dir}/bold_tsnr.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN2PAM50_t2s.nii.gz \
				-o ${func_dir}/bold_tsnr_templatespace.nii.gz
		sct_apply_transfo -i ${func_dir}/bold_denoised_tsnr.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN2PAM50_t2s.nii.gz \
				-o ${func_dir}/bold_denoised_tsnr_templatespace.nii.gz
	done
fi


# -------------------------------------------------------------------------------------------
# 5. Create an average tSNR image in template space
# -------------------------------------------------------------------------------------------
if [ $create_average_tsnr -eq 1 ]; then
	group_dir=${project_dir}/derivatives/group_analysis
	mkdir -p $group_dir
	cd ${group_dir}
	
	echo "Creating average tSNR images"
	sub_raw_tsnr=()
	sub_denoised_tsnr=()
	for subject in {1..62}; do
		printf -v sub "%03d" $subject
		if [[ "$subject" == 29 || "$subject" == 60 ]]; then
			echo "Subject has no functional data"
			continue
		elif [[ "$subject" == 56 ]]; then
			echo "Subject has severe artifacts"
			continue
		fi	
		func_dir=${project_dir}/derivatives/sub-${sub}/func/sosGREp3
		sub_raw_tsnr+="${func_dir}/bold_tsnr_templatespace.nii.gz "
		sub_denoised_tsnr+="${func_dir}/bold_denoised_tsnr_templatespace.nii.gz "

	done
	fslmerge -t Raw_tsnr_all_subs_sosGREp3.nii.gz ${sub_raw_tsnr[@]}
	fslmaths Raw_tsnr_all_subs_sosGREp3.nii.gz -Tmean Raw_tsnr_all_subs_sosGREp3_mean.nii.gz
	fslmerge -t Denoised_tsnr_all_subs_sosGREp3.nii.gz ${sub_denoised_tsnr[@]}
	fslmaths Denoised_tsnr_all_subs_sosGREp3.nii.gz -Tmean Denoised_tsnr_all_subs_sosGREp3_mean.nii.gz

fi
