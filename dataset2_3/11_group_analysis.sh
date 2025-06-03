#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
normalize_func=1
copes2template=0
tsnr2template=0
tsnr2template_100=0
res4d2template=0
create_average_epi=0
create_average_tsnr=0
create_average_tsnr_100=0
signal_limits=0
collect_copes=0
create_masks=0
cut_masks=0
randomize=0
randomize_atlas=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder where raw and output data is located
code_dir=/data/pt_02661_raw/Heatpain/code  # code folder with subfolder helper_functions
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located
mask_dir=$project_dir/derivatives/masks

# add FSL to the path
FSLDIR=${fsl_dir}
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# add SCT to the path
export PATH="${sct_dir}/bin:${PATH}"

#process_subject() {
#  	sub=$1
for subject in {1..41}; do

	printf -v sub "%02d" $subject
	echo "subject: " $sub	
	
	func_dir=$project_dir/derivatives/sub-sspr${sub}/func
	anat_dir=$project_dir/derivatives/sub-sspr${sub}/anat
	
	func_anat_dir=$project_dir/derivatives/sub-sspr${sub}/func_anat
	mkdir -p $func_anat_dir


	# -------------------------------------------------------------------------------------------
	# 1. Register functional image to template space
	# -------------------------------------------------------------------------------------------
	if [ $normalize_func -eq 1 ]; then
		seclevel_dir=${func_dir}/second_level
		cd $func_anat_dir
		
		# register the moco mean to the template by using the warping fields from anat2template transformation
		# in the anat folder there is already the optimal of the two warping fields
		
		mkdir -p "optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn"
		
		sct_register_multimodal \
		-i  ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
		-d ${seclevel_dir}/MOCO_MEAN_all_runs.nii.gz \
		-iseg ${sct_dir}/data/PAM50/template/PAM50_cord.nii.gz \
		-dseg ${seclevel_dir}/MOCO_MEAN_deepseg_manual.nii.gz \
		-param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,smooth=1,slicewise=1,iter=3 \
		-initwarp ${anat_dir}/warp_template2anat.nii.gz \
		-initwarpinv ${anat_dir}/warp_anat2template.nii.gz \
		-x spline -ofolder ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn
			
	fi
	
    	# -------------------------------------------------------------------------------------------
	# 2. Transform copes into template space
	# -------------------------------------------------------------------------------------------
	if [ $copes2template -eq 1 ]; then
		for design in block_design_2mm_susan onset_design_2mm_susan block_design_2mm_susan_vasa onset_design_2mm_susan_vasa block_design_2mm_susan_vasa_even_odd onset_design_2mm_susan_vasa_even_odd block_design_no_smooth_vasa onset_design_no_smooth_vasa; do

			echo "Warping cope image into template space ${design}"
			out_dir=${project_dir}/derivatives/sub-sspr${sub}/copes_templatespace/${design}
			mkdir -p $out_dir
			cd $out_dir
			feat_dir=${func_dir}/second_level/${design}.gfeat
			if [[ $design == block* ]]; then
				if [[ $design == *even_odd ]]; then
					num_copes=2
				else
					num_copes=1
				fi
			else
				if [[ $design == *even_odd ]]; then
					num_copes=5
				else
					num_copes=3
				fi
			fi
			for cope in $(seq 1 $num_copes); do
				sct_apply_transfo -i ${feat_dir}/cope${cope}.feat/stats/cope1.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN_all_runs2PAM50_t2s.nii.gz \
				-o cope${cope}.nii.gz
			done
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 3. Transform tSNR images into template space
	# -------------------------------------------------------------------------------------------
	if [ $tsnr2template -eq 1 ]; then
		echo "Building a mean of all runs' tSNR images"
		
		# when some runs were not recorded
		if [[ $subject == 25 || $subject == 27 || $subject == 29 || $subject == 6 ]]; then
			num_runs=2
		elif [[ $subject == 31 ]]; then
			num_runs=3
		else
			num_runs=4
		fi
		
		func_dir=${project_dir}/derivatives/sub-sspr${sub}/func
		raw_tsnr=()
		denoised_tsnr=()
		#moco_tsnr=()
		for run in  $(seq 1 $num_runs); do
			run_dir=${func_dir}/run-${run}
			raw_tsnr+="${run_dir}/bold_tsnr.nii.gz "
			denoised_tsnr+="${run_dir}/bold_denoised_tsnr.nii.gz "
		done
		sec_level_dir=${func_dir}/second_level
		cd ${sec_level_dir}
		fslmerge -t Raw_tsnr_all_runs.nii.gz ${raw_tsnr[@]}
		fslmaths Raw_tsnr_all_runs.nii.gz -Tmean Raw_tsnr_all_runs_mean.nii.gz
		fslmerge -t Denoised_tsnr_all_runs.nii.gz ${denoised_tsnr[@]}
		fslmaths Denoised_tsnr_all_runs.nii.gz -Tmean Denoised_tsnr_all_runs_mean.nii.gz
		#fslmerge -t MOCO_tsnr_all_runs.nii.gz ${moco_tsnr[@]}
		#fslmaths MOCO_tsnr_all_runs.nii.gz -Tmean MOCO_tsnr_all_runs_mean.nii.gz
		
		echo "Warping tSNR image into template space"
		sct_apply_transfo -i Raw_tsnr_all_runs_mean.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN_all_runs2PAM50_t2s.nii.gz \
				-o Raw_tsnr_all_runs_mean_templatespace.nii.gz
		sct_apply_transfo -i Denoised_tsnr_all_runs_mean.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN_all_runs2PAM50_t2s.nii.gz \
				-o Denoised_tsnr_all_runs_mean_templatespace.nii.gz
	
	fi
	
	# -------------------------------------------------------------------------------------------
	# 3a. Transform tSNR images into template space --> this time only the first 100 volumes
	# to compare it to screening dataset
	# -------------------------------------------------------------------------------------------
	if [ $tsnr2template_100 -eq 1 ]; then
		echo "Taking the first 100 volumes"
		
		func_dir=${project_dir}/derivatives/sub-sspr${sub}/func
		run_dir=${func_dir}/run-1
		cd ${run_dir}
		fslroi ${run_dir}/bold.nii.gz ${run_dir}/bold_100.nii.gz 0 100
		fslroi ${run_dir}/bold_denoised.nii.gz ${run_dir}/bold_denoised_100.nii.gz 0 100
		
		echo "Calculating tsnr again"
		fslmaths bold_100.nii.gz -Tmean bold_100_mean.nii.gz
		fslmaths bold_100.nii.gz -Tstd bold_100_sd.nii.gz
		fslmaths bold_100_mean.nii.gz -div bold_100_sd.nii.gz bold_100_tsnr.nii.gz
		
		fslmaths bold_denoised_100.nii.gz -Tmean bold_denoised_100_mean.nii.gz
		fslmaths bold_denoised_100.nii.gz -Tstd bold_denoised_100_sd.nii.gz
		fslmaths bold_denoised_100_mean.nii.gz -div bold_denoised_100_sd.nii.gz bold_denoised_100_tsnr.nii.gz

		echo "Warping tSNR image into template space"
		sct_apply_transfo -i ${run_dir}/bold_100_tsnr.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN_all_runs2PAM50_t2s.nii.gz \
				-o ${run_dir}/bold_100_tsnr_templatespace.nii.gz
		sct_apply_transfo -i ${run_dir}/bold_denoised_100_tsnr.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN_all_runs2PAM50_t2s.nii.gz \
				-o ${run_dir}/bold_denoised_100_tsnr_templatespace.nii.gz
	
	fi
	
	# -------------------------------------------------------------------------------------------
	# 4. Transform dummy GLM res4d images into template space
	# -------------------------------------------------------------------------------------------
	if [ $res4d2template -eq 1 ]; then
		echo "Warp the res4d images of the dummy GLM into template space"
		
		# when some runs were not recorded
		if [[ $subject == 25 || $subject == 27 || $subject == 29 || $subject == 6 ]]; then
			num_runs=2
		elif [[ $subject == 31 ]]; then
			num_runs=3
		else
			num_runs=4
		fi
		
		func_dir=${project_dir}/derivatives/sub-sspr${sub}/func
		for run in  $(seq 1 $num_runs); do
			feat_dir=${func_dir}/run-${run}/feat_dummy_2mm_susan.feat/stats
			cd $feat_dir
			sct_apply_transfo -i ${feat_dir}/res4d.nii.gz  \
				-d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
				-w ${func_anat_dir}/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/warp_MOCO_MEAN_all_runs2PAM50_t2s.nii.gz \
				-o ${feat_dir}/res4d_templatespace.nii.gz
		done
	
	fi
	
#}
done

# Loop across subjects for data preparation in parallel
#for subject in {1..41}; do
#	printf -v sub "%02d" $subject
#	process_subject $sub & 
#done

# Wait for all background processes to finish
#wait
#echo "All parallel processes completed."


# -------------------------------------------------------------------------------------------
# 4. Create an average EPI in template space
# -------------------------------------------------------------------------------------------
if [ $create_average_epi -eq 1 ]; then
	group_dir=${project_dir}/derivatives/group_analysis
	mkdir -p $group_dir
	cd ${group_dir}
	
	echo "Creating average EPI in dataset 2"
	sub_epi=()
	for subject in {1..16}; do
		printf -v sub "%02d" $subject
		main_dir=${project_dir}/derivatives/sub-sspr${sub}
		sub_epi+="$main_dir/func_anat/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/MOCO_MEAN_all_runs_reg.nii.gz "
	done
	fslmerge -t EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_merged_dataset2.nii.gz ${sub_epi[@]}
	fslmaths EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_merged_dataset2.nii.gz -Tmean EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_mean_dataset2.nii.gz
	fslmaths EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_merged_dataset2.nii.gz -Tstd EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_std_dataset2.nii.gz
	fslmaths EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_mean_dataset2.nii.gz -div EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_std_dataset2.nii.gz EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_tsnr_dataset2.nii.gz
	
	echo "Creating average EPI in dataset 3"
	sub_epi=()
	for subject in {17..41}; do
		printf -v sub "%02d" $subject
		main_dir=${project_dir}/derivatives/sub-sspr${sub}
		sub_epi+="$main_dir/func_anat/optimal_anat_norm_EPI_seg_centermass_seg_bsplinesyn/MOCO_MEAN_all_runs_reg.nii.gz "
	done
	fslmerge -t EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_merged_dataset3.nii.gz ${sub_epi[@]}
	fslmaths EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_merged_dataset3.nii.gz -Tmean EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_mean_dataset3.nii.gz
	fslmaths EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_merged_dataset3.nii.gz -Tstd EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_std_dataset3.nii.gz
	fslmaths EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_mean_dataset3.nii.gz -div EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_std_dataset3.nii.gz EPI_MOCO_optimal_norm_segcentermass_segbsplinesyn_tsnr_dataset3.nii.gz

fi

# -------------------------------------------------------------------------------------------
# 5. Create an average tSNR image in template space
# -------------------------------------------------------------------------------------------
if [ $create_average_tsnr -eq 1 ]; then
	group_dir=${project_dir}/derivatives/group_analysis
	mkdir -p $group_dir
	cd ${group_dir}
	
	echo "Creating average tSNR images in dataset 2"
	sub_raw_tsnr=()
	sub_denoised_tsnr=()
	sub_moco_tsnr=()
	for subject in {1..16}; do
		printf -v sub "%02d" $subject
		main_dir=${project_dir}/derivatives/sub-sspr${sub}
		sub_raw_tsnr+="$main_dir/func/second_level/Raw_tsnr_all_runs_mean_templatespace.nii.gz "
		sub_denoised_tsnr+="$main_dir/func/second_level/Denoised_tsnr_all_runs_mean_templatespace.nii.gz "
		sub_moco_tsnr+="$main_dir/func/second_level/MOCO_tsnr_all_runs_mean_templatespace.nii.gz "
	done
	fslmerge -t Raw_tsnr_all_subs_dataset2.nii.gz ${sub_raw_tsnr[@]}
	fslmaths Raw_tsnr_all_subs_dataset2.nii.gz -Tmean Raw_tsnr_all_subs_dataset2_mean.nii.gz
	fslmerge -t Denoised_tsnr_all_subs_dataset2.nii.gz ${sub_denoised_tsnr[@]}
	fslmaths Denoised_tsnr_all_subs_dataset2.nii.gz -Tmean Denoised_tsnr_all_subs_dataset2_mean.nii.gz
	fslmerge -t MOCO_tsnr_all_subs_dataset2.nii.gz ${sub_moco_tsnr[@]}
	fslmaths MOCO_tsnr_all_subs_dataset2.nii.gz -Tmean MOCO_tsnr_all_subs_dataset2_mean.nii.gz
	
	echo "Creating average tSNR images in dataset 3"
	sub_raw_tsnr=()
	sub_denoised_tsnr=()
	sub_moco_tsnr=()
	for subject in {17..41}; do
		printf -v sub "%02d" $subject
		main_dir=${project_dir}/derivatives/sub-sspr${sub}
		sub_raw_tsnr+="$main_dir/func/second_level/Raw_tsnr_all_runs_mean_templatespace.nii.gz "
		sub_denoised_tsnr+="$main_dir/func/second_level/Denoised_tsnr_all_runs_mean_templatespace.nii.gz "
		sub_moco_tsnr+="$main_dir/func/second_level/MOCO_tsnr_all_runs_mean_templatespace.nii.gz "
	done
	fslmerge -t Raw_tsnr_all_subs_dataset3.nii.gz ${sub_raw_tsnr[@]}
	fslmaths Raw_tsnr_all_subs_dataset3.nii.gz -Tmean Raw_tsnr_all_subs_dataset3_mean.nii.gz
	fslmerge -t Denoised_tsnr_all_subs_dataset3.nii.gz ${sub_denoised_tsnr[@]}
	fslmaths Denoised_tsnr_all_subs_dataset3.nii.gz -Tmean Denoised_tsnr_all_subs_dataset3_mean.nii.gz
	fslmerge -t MOCO_tsnr_all_subs_dataset3.nii.gz ${sub_moco_tsnr[@]}
	fslmaths MOCO_tsnr_all_subs_dataset3.nii.gz -Tmean MOCO_tsnr_all_subs_dataset3_mean.nii.gz
	
	echo "Creating average tSNR images in both datasets combined"
	sub_raw_tsnr=()
	sub_denoised_tsnr=()
	sub_moco_tsnr=()
	for subject in {1..41}; do
		if [[ $subject -eq 1 || $subject -eq 12 ]]; then
			continue
		else
			printf -v sub "%02d" $subject
			main_dir=${project_dir}/derivatives/sub-sspr${sub}
			sub_raw_tsnr+="$main_dir/func/second_level/Raw_tsnr_all_runs_mean_templatespace.nii.gz "
			sub_denoised_tsnr+="$main_dir/func/second_level/Denoised_tsnr_all_runs_mean_templatespace.nii.gz "
			sub_moco_tsnr+="$main_dir/func/second_level/MOCO_tsnr_all_runs_mean_templatespace.nii.gz "
		fi
	done
	fslmerge -t Raw_tsnr_all_subs_both.nii.gz ${sub_raw_tsnr[@]}
	fslmaths Raw_tsnr_all_subs_both.nii.gz -Tmean Raw_tsnr_all_subs_both_mean.nii.gz
	fslmerge -t Denoised_tsnr_all_subs_both.nii.gz ${sub_denoised_tsnr[@]}
	fslmaths Denoised_tsnr_all_subs_both.nii.gz -Tmean Denoised_tsnr_all_subs_both_mean.nii.gz
	fslmerge -t MOCO_tsnr_all_subs_both.nii.gz ${sub_moco_tsnr[@]}
	fslmaths MOCO_tsnr_all_subs_both.nii.gz -Tmean MOCO_tsnr_all_subs_both_mean.nii.gz
fi

# -------------------------------------------------------------------------------------------
# 5a. Create an average tSNR image in template space --> only first 100 volumes
# -------------------------------------------------------------------------------------------
if [ $create_average_tsnr_100 -eq 1 ]; then
	group_dir=${project_dir}/derivatives/group_analysis
	mkdir -p $group_dir
	cd ${group_dir}
	
	echo "Creating average tSNR images in dataset 2"
	sub_raw_tsnr=()
	sub_denoised_tsnr=()
	for subject in {1..16}; do
		printf -v sub "%02d" $subject
		main_dir=${project_dir}/derivatives/sub-sspr${sub}
		sub_raw_tsnr+="$main_dir/func/run-1/bold_100_tsnr_templatespace.nii.gz "
		sub_denoised_tsnr+="$main_dir/func/run-1/bold_denoised_100_tsnr_templatespace.nii.gz "

	done
	fslmerge -t Raw_tsnr_100_all_subs_dataset2.nii.gz ${sub_raw_tsnr[@]}
	fslmaths Raw_tsnr_100_all_subs_dataset2.nii.gz -Tmean Raw_tsnr_100_all_subs_dataset2_mean.nii.gz
	fslmerge -t Denoised_tsnr_100_all_subs_dataset2.nii.gz ${sub_denoised_tsnr[@]}
	fslmaths Denoised_tsnr_100_all_subs_dataset2.nii.gz -Tmean Denoised_tsnr_100_all_subs_dataset2_mean.nii.gz

fi

# -------------------------------------------------------------------------------------------
# 6. create masks to check where all the subjects have signal so that we can limit randomize to that range
# -------------------------------------------------------------------------------------------
if [ $signal_limits -eq 1 ]; then

	echo "Checking where the datasets have signal"
	for dataset in dataset2 dataset3 both; do
		group_dir=${project_dir}/derivatives/group_analysis/${dataset}
		mkdir -p $group_dir
		
		sub_cope_mask=()
		if [ ${dataset} = "dataset2" ]; then
			start=1
			end=16
			total_num=16
		elif [ ${dataset} = "dataset3" ]; then
			start=17
			end=41
			total_num=25
		elif [ ${dataset} = "both" ]; then
			start=1
			end=41
			total_num=39
		fi
		for subject in $(seq $start $end); do
			echo "Subject ${subject}"
			printf -v sub "%02d" $subject
			
			# exclude the subjects that have been measured twice
			if [[ (${dataset} = "both") && ($subject -eq 1 || $subject -eq 12) ]];then
				continue
			else
				feat_dir=${project_dir}/derivatives/sub-sspr${sub}/copes_templatespace/block_design_2mm_susan
				cd $feat_dir
				
				# as the feats can be positive and negative, binarize in both directions
				fslmaths cope1.nii.gz -mul -1 -bin cope1_neg_bin.nii.gz
				fslmaths cope1.nii.gz -bin cope1_pos_bin.nii.gz
				fslmaths cope1_neg_bin.nii.gz -add cope1_pos_bin.nii.gz cope1_mask.nii.gz

				# delete intermediate results
				rm -f cope1_pos_bin.nii.gz
		   		rm -f cope1_neg_bin.nii.gz
			   	
			   	if [ ${subject} -eq ${end} ];then
			   		sub_cope_mask+="$feat_dir/cope1_mask.nii.gz "
			   	else
			   		sub_cope_mask+="$feat_dir/cope1_mask.nii.gz -add "
			   	fi
			fi
		done
		
		# add all copes
		cd $group_dir
		fslmaths ${sub_cope_mask[@]} all_copes1_sum.nii.gz

		# now divide by number of subjects so that 1 always means all subjects have data there, less than 1 means there are gaps
		fslmaths all_copes1_sum.nii.gz -div ${total_num} all_copes1_normed_sum.nii.gz
		
		echo "fsleyes ${group_dir}/all_copes1_normed_sum.nii.gz ${sct_dir}/data/PAM50/template/PAM50_cord.nii.gz"

	done
fi
# -------------------------------------------------------------------------------------------
# now you can lay the created normed sum image and the cord mask from PAM50 together and check where the borders are
# (the acceptable range of slices needs to have a value of 1 everywhere in the cord)
# fsleyes all_copes1_normed_sum.nii.gz ${sct_dir}/data/PAM50/template/PAM50_cord.nii.gz
# dataset 2: 	785 - 855
# dataset 3: 	790 - 859
# both: 	790 - 855
# -------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# 7. Collect cope images
# -------------------------------------------------------------------------------------------
if [ $collect_copes = 1 ]; then
	for design in block_design_2mm_susan onset_design_2mm_susan block_design_2mm_susan_vasa onset_design_2mm_susan_vasa block_design_2mm_susan_vasa_even_odd onset_design_2mm_susan_vasa_even_odd block_design_no_smooth_vasa onset_design_no_smooth_vasa; do

		echo "${design}"
		if [[ $design == block* ]]; then
			if [[ $design == *even_odd ]]; then
				num_copes=2
			else
				num_copes=1
			fi
		else
			if [[ $design == *even_odd ]]; then
				num_copes=5
			else
				num_copes=3
			fi
		fi
		
		for dataset in dataset2 dataset3 both; do
			echo "${dataset}"
			group_dir=${project_dir}/derivatives/group_analysis/${dataset}/${design}
			mkdir -p $group_dir
			
			# then collect all copes
			for cope in $(seq 1 $num_copes); do
				echo "Cope ${cope}"
				
				sub_cope=()
				if [ ${dataset} = "dataset2" ]; then
					start=1
					end=16
				elif [ ${dataset} = "dataset3" ]; then
					start=17
					end=41
				elif [ ${dataset} = "both" ]; then
					start=1
					end=41
				fi
				for subject in $(seq $start $end); do
					echo "Subject ${subject}"
					printf -v sub "%02d" $subject
					# exclude the subjects that have been measured twice
					if [[ (${dataset} = "both") && ($subject -eq 1 || $subject -eq 12) ]];then
						continue
					else
						feat_dir=${project_dir}/derivatives/sub-sspr${sub}/copes_templatespace/${design}
					   	sub_cope+="$feat_dir/cope${cope}.nii.gz "
					fi
				done
				cd $group_dir
				fslmerge -t cope${cope}_merged.nii.gz ${sub_cope[@]}
			done
		done
	done
fi


# -------------------------------------------------------------------------------------------
# 8. Create template space PAM50 masks
# -------------------------------------------------------------------------------------------
if [ $create_masks -eq 1 ]; then
	
	mask_dir=${project_dir}/derivatives/masks/
	echo "building horn masks"

	# copy and binarize masks for horns (VHL, VHR, DHL, DHR)
	declare -A atlas_files=(["vh_left"]=30 ["vh_right"]=31 ["dh_left"]=34 ["dh_right"]=35)
	for region in "${!atlas_files[@]}"; do
		cp ${sct_dir}/data/PAM50/atlas/PAM50_atlas_${atlas_files[$region]}.nii.gz $mask_dir/${region}.nii.gz
		fslmaths $mask_dir/${region}.nii.gz -thr 0.5 -bin $mask_dir/${region}.nii.gz
	done
	# copy and binarize spinal cord levels
	cp ${sct_dir}/data/PAM50/template/PAM50_spinal_levels.nii.gz $mask_dir/PAM50_spinal_levels.nii.gz
	for level in {5..8};do
		fslmaths $mask_dir/PAM50_spinal_levels.nii.gz -thr ${level} -uthr ${level} -bin $mask_dir/C${level}.nii.gz
	done
	
	# then use the levels to mask the horns so that they are limited to specific levels
	for region in "${!atlas_files[@]}"; do
		for level in {5..8};do
			fslmaths $mask_dir/${region}.nii.gz -mas $mask_dir/C${level}.nii.gz $mask_dir/${region}_C${level}.nii.gz
		done
	done
fi

# -------------------------------------------------------------------------------------------
# 9. Cut masks depending on the dataset used
# -------------------------------------------------------------------------------------------
if [ $cut_masks -eq 1 ]; then
	
	mask_dir=${project_dir}/derivatives/masks/
	
	for dataset in dataset2 dataset3 both; do
		echo "Cutting masks for ${dataset}"
		
		out_path=${project_dir}/derivatives/group_analysis/${dataset}
		cd ${out_path}
		mkdir -p masks
		
		if [[ ${dataset} = "dataset2" ]]; then
			cutting="785 70"
		elif [ ${dataset} = "dataset3" ]; then
			cutting="790 69"
		elif [ ${dataset} = "both" ]; then
			cutting="790 65"
		fi
			
		declare -A atlas_files=(["vh_left"]=30 ["vh_right"]=31 ["dh_left"]=34 ["dh_right"]=35)
		for region in "${!atlas_files[@]}"; do
			fslroi ${mask_dir}/${region}.nii.gz ${out_path}/masks/${region}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
		done

		for level in {5..8};do
			fslroi ${mask_dir}/C${level}.nii.gz ${out_path}/masks/C${level}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
		done
		
		for region in "${!atlas_files[@]}"; do
			for level in {5..8};do
				fslroi ${mask_dir}/${region}_C${level}.nii.gz ${out_path}/masks/${region}_C${level}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
			done
		done
		
		fslroi ${mask_dir}/PAM50_gm_bin.nii.gz ${out_path}/masks/PAM50_gm_bin_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
		
		for mask in atlas_C6_L1L2_left_PAM50 atlas_C6_L3L4_left_PAM50 atlas_C6_L5L6_left_PAM50 atlas_C6_L7ToL10_left_PAM50; do
			fslroi ${mask_dir}/${mask}.nii.gz ${out_path}/masks/${mask}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
		done
		for mask in atlas_C6_L1L2_right_PAM50 atlas_C6_L3L4_right_PAM50 atlas_C6_L5L6_right_PAM50 atlas_C6_L7ToL10_right_PAM50; do
			fslroi ${mask_dir}/${mask}.nii.gz ${out_path}/masks/${mask}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
		done
	done
fi

# -------------------------------------------------------------------------------------------
# 9. Run randomize
# -------------------------------------------------------------------------------------------
if [ $randomize -eq 1 ]; then
	
	mask_dir=${project_dir}/derivatives/masks/
	
	for design in block_design_2mm_susan onset_design_2mm_susan block_design_2mm_susan_vasa onset_design_2mm_susan_vasa block_design_2mm_susan_vasa_even_odd onset_design_2mm_susan_vasa_even_odd block_design_unsmoothed_vasa onset_design_unsmoothed_vasa; do
	
		echo "${design}"
		if [[ $design == block* ]]; then
			if [[ $design == *even_odd ]]; then
				num_copes=2
			else
				num_copes=1
			fi
		else
			if [[ $design == *even_odd ]]; then
				num_copes=5
			else
				num_copes=3
			fi
		fi

		for dataset in dataset2 dataset3 both; do
			echo "${dataset}"
			if [[ ${dataset} = "dataset2" ]]; then
				cutting="785 70"
			elif [ ${dataset} = "dataset3" ]; then
				cutting="790 69"
			elif [ ${dataset} = "both" ]; then
				cutting="790 65"
			fi

			group_dir=${project_dir}/derivatives/group_analysis/${dataset}/${design}
			cd $group_dir
			
			mkdir -p PAM50_cord
			cd PAM50_cord
			# prepare a cord mask as a mask, cut into the window that is needed so that all subjects have data there
			fslroi $sct_dir/data/PAM50/template/PAM50_cord.nii.gz PAM50_cord_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
			
			for cope in $(seq 1 $num_copes); do
				echo "Cope ${cope}" 
				cd $group_dir
				# cut the copes into the window that is needed so that all subjects have data there
				fslroi cope${cope}_merged.nii.gz cope${cope}_merged_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
				
				cd PAM50_cord
				randomise -i ${group_dir}/cope${cope}_merged_cut.nii.gz -m PAM50_cord_cut.nii.gz -v 2 -1 --uncorrp -T -x -c 2.3 -C 2.3 -o cope${cope}
				fslmaths cope${cope}_tfce_p_tstat1.nii.gz -thr 0.99 -bin uncorr_p_mask_cope${cope}_99.nii.gz
				fslmaths cope${cope}_tstat1.nii.gz -mas uncorr_p_mask_cope${cope}_99.nii.gz uncorr_t_map_cope${cope}_99.nii.gz
				cluster -i cope${cope}_tfce_corrp_tstat1 -t 0.95 -c cope${cope}_tstat1 --scalarname="1-p" > cluster_cope${cope}_tfce.csv
				fslmaths cope${cope}_tfce_corrp_tstat1.nii.gz -thr 0.95 -bin corr_p_mask_cope${cope}.nii.gz
				fslmaths cope${cope}_tstat1.nii.gz -mas corr_p_mask_cope${cope}.nii.gz corr_t_map_cope${cope}.nii.gz
				
				# and then test target ROI
				mask=dh_left_C6
				cd $group_dir
				mkdir -p ${mask}
				cp ${mask_dir}/${mask}.nii.gz ${mask}/${mask}.nii.gz
				cd ${mask}
				fslroi ${mask}.nii.gz ${mask}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
				randomise -i ${group_dir}/cope${cope}_merged_cut.nii.gz -m ${mask}_cut.nii.gz -v 2 -1 --uncorrp -T -x -c 2.3 -C 2.3 -o cope${cope}
				fslmaths cope${cope}_tfce_p_tstat1.nii.gz -thr 0.99 -bin uncorr_p_mask_cope${cope}_99.nii.gz
				fslmaths cope${cope}_tstat1.nii.gz -mas uncorr_p_mask_cope${cope}_99.nii.gz uncorr_t_map_cope${cope}_99.nii.gz
				cluster -i cope${cope}_tfce_corrp_tstat1 -t 0.95 -c cope${cope}_tstat1 --scalarname="1-p" > cluster_cope${cope}_tfce.csv
				fslmaths cope${cope}_tfce_corrp_tstat1.nii.gz -thr 0.95 -bin corr_p_mask_cope${cope}.nii.gz
				fslmaths cope${cope}_tstat1.nii.gz -mas corr_p_mask_cope${cope}.nii.gz corr_t_map_cope${cope}.nii.gz
			done
		done
	done
fi

# -------------------------------------------------------------------------------------------
# 10. Run randomize within the new atlas
# -------------------------------------------------------------------------------------------
if [ $randomize_atlas -eq 1 ]; then
	
	mask_dir=${project_dir}/derivatives/atlas/
	
	for design in block_design_2mm_susan_vasa onset_design_2mm_susan_vasa block_design_no_smooth_vasa onset_design_no_smooth_vasa; do

		echo "${design}"
		if [[ $design == block* ]]; then
			if [[ $design == *even_odd ]]; then
				num_copes=2
			else
				num_copes=1
			fi
		else
			if [[ $design == *even_odd ]]; then
				num_copes=5
			else
				num_copes=3
			fi
		fi
		
		for dataset in dataset2 dataset3 both; do
			echo "${dataset}"
			if [[ ${dataset} = "dataset2" ]]; then
				cutting="785 70"
			elif [ ${dataset} = "dataset3" ]; then
				cutting="790 69"
			elif [ ${dataset} = "both" ]; then
				cutting="790 65"
			fi

			group_dir=${project_dir}/derivatives/group_analysis/${dataset}/${design}
			cd $group_dir
			
			for cope in $(seq 1 $num_copes); do
				echo "Cope ${cope}" 
				
				# and then test several masks
				for mask in atlas_C6_L1L2_left_PAM50 atlas_C6_L3L4_left_PAM50 atlas_C6_L5L6_left_PAM50 atlas_C6_L7ToL10_left_PAM50; do
					cd $group_dir
					mkdir -p ${mask}
					
					cp ${mask_dir}/${mask}.nii.gz ${mask}/${mask}.nii.gz
					cd ${mask}
					fslroi ${mask}.nii.gz ${mask}_cut.nii.gz 0 -1 0 -1 ${cutting} 0 -1
					randomise -i ${group_dir}/cope${cope}_merged_cut.nii.gz -m ${mask}_cut.nii.gz -v 2 -1 --uncorrp -T -x -c 2.3 -C 2.3 -o cope${cope}
					fslmaths cope${cope}_tfce_p_tstat1.nii.gz -thr 0.95 -bin uncorr_p_mask_cope${cope}.nii.gz
					fslmaths cope${cope}_tstat1.nii.gz -mas uncorr_p_mask_cope${cope}.nii.gz uncorr_t_map_cope${cope}.nii.gz
					fslmaths cope${cope}_tfce_p_tstat1.nii.gz -thr 0.99 -bin uncorr_p_mask_cope${cope}_99.nii.gz
					fslmaths cope${cope}_tstat1.nii.gz -mas uncorr_p_mask_cope${cope}_99.nii.gz uncorr_t_map_cope${cope}_99.nii.gz
					cluster -i cope${cope}_tfce_corrp_tstat1 -t 0.95 -c cope${cope}_tstat1 --scalarname="1-p" > cluster_cope${cope}_tfce.csv
					fslmaths cope${cope}_tfce_corrp_tstat1.nii.gz -thr 0.95 -bin corr_p_mask_cope${cope}.nii.gz
					fslmaths cope${cope}_tstat1.nii.gz -mas corr_p_mask_cope${cope}.nii.gz corr_t_map_cope${cope}.nii.gz
				done
			done
		done
	done
fi





