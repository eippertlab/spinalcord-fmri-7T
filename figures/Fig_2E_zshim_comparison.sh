#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
reconstruct_zshim=1
segmentation=1
warp_to_template=1
extract_signal_intensity=1
average_images=1

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder
group_dir=$project_dir/derivatives/group_analysis/zshim
code_dir=${project_dir}/final_pipeline  # code folder with subfolder helper_functions
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located

# add FSL to the path
FSLDIR=${fsl_dir}
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# add SCT to the path
export PATH="${sct_dir}/bin:${PATH}"


for subject in {17..41}; do
	printf -v sub "%02d" $subject
	echo "subject: " $sub
	
	raw_sub_dir=$project_dir/raw_data/dataset3/sub-sspr${sub}
	out_sub_dir=$project_dir/derivatives/sub-sspr${sub}
	zshim_dir=$out_sub_dir/func/zshim
	anat_dir=$out_sub_dir/anat
	
	mkdir -p ${zshim_dir}
	
	# -------------------------------------------------------------------------------------------
	# 1. Get the data without z-shim and with z-shim for the specific phase encoding direction
	# -------------------------------------------------------------------------------------------
	if [ $reconstruct_zshim = 1 ]; then
		echo "Checking phase encoding direction"
		row_name=sub-sspr${sub}
		column_name=EPI_encoding_dir
		file=$project_dir/raw_data/dataset3/participants.tsv
		column_index=$(awk -v col="$column_name" 'NR==1{for (i=1; i<=NF; i++) if ($i == col) { print i; exit }}' "$file")
		PE=$(awk -v row="$row_name" -v col_index="$column_index" 'BEGIN{FS="\t"; found=0} {if ($1 == row) { found=1; print $(col_index); exit; } } END{ if (found==0) print "Value not found" }' "$file")
		echo PE: ${PE}
		cd $zshim_dir
		ref_scan=${raw_sub_dir}/func/sub-sspr${sub}_task-rest_acq-zref_dir-${PE}_part-mag_bold.nii.gz
		echo "Getting data without zshim"
		fslroi $ref_scan ${zshim_dir}/no_zshim.nii.gz 15 1
		
		echo "Reconstructing data with zshim"
		# open the file with the selected zshims per slice
		zshim_values=${raw_sub_dir}/ZShim_AutoRef_Thr100_${PE}.csv
		mapfile -t zshims < <(awk '{
			# Check if the number is scientific notation
			if ($1 ~ /[eE]/) {
				printf "%.0f\n", $1
			} else {
				print $1
			}
			}' "${zshim_values}")
		# split images (zhims in time and slices)
		mkdir -p tmp
		cp $ref_scan tmp/
		cd tmp
		fslsplit $ref_scan "zref" -t
		for j in zref0*; do
			fslsplit $j $(remove_ext $j) -z
		done
		# create file name matrix based on zshim moment and slice number
		filenames=""
		slice=0
		for val in ${zshims[@]};do
			filenames+="zref$(printf '%04d' $val)$(printf '%04d' $slice).nii.gz "
			slice=$((slice+1))
		done
		# create reconstructed zshim image and delete tmp folder
		fslmerge -z zshim_recon ${filenames}
		cd $zshim_dir
		cp tmp/zshim_recon.nii.gz $zshim_dir/
		rm -Rf tmp
	fi
	
	# -------------------------------------------------------------------------------------------
	# 2. Do a cord segmentation in the reconstructed and no z-shim data as well as in the z-ref
	# -------------------------------------------------------------------------------------------
	if [ $segmentation = 1 ]; then
		echo "Checking phase encoding direction"
		row_name=sub-sspr${sub}
		column_name=EPI_encoding_dir
		file=$project_dir/raw_data/dataset3/participants.tsv
		column_index=$(awk -v col="$column_name" 'NR==1{for (i=1; i<=NF; i++) if ($i == col) { print i; exit }}' "$file")
		PE=$(awk -v row="$row_name" -v col_index="$column_index" 'BEGIN{FS="\t"; found=0} {if ($1 == row) { found=1; print $(col_index); exit; } } END{ if (found==0) print "Value not found" }' "$file")
		echo PE: ${PE}
		
		echo "Building a mean of the zshim ref scan"
		cd ${zshim_dir}
		ref_scan=${raw_sub_dir}/func/sub-sspr${sub}_task-rest_acq-zref_dir-${PE}_part-mag_bold.nii.gz
		
    		#fslmaths ${ref_scan} -Tmean zref_mean.nii.gz
		#sct_deepseg_sc -i zref_mean.nii.gz -c t2s
		
		echo "Extracting signal in this cord mask from zshim and no zshim"
		for mode in zshim_recon no_zshim; do
		    	#sct_deepseg_sc -i ${mode}.nii.gz -c t2s
		    	fslmeants -i ${mode}.nii.gz -m zref_mean_seg.nii.gz --showall -o ${mode}_signal_all.txt
		done
		
	fi
	
	# -------------------------------------------------------------------------------------------
	# for sub-sspr30 and sub-sspr31 I needed to correct the segmentations as they did not have anything in some slices
	# -------------------------------------------------------------------------------------------
	
	# -------------------------------------------------------------------------------------------
	# 3. Warp the z-ref mean to template space, then use this warp for the zshim/no zshim scans
	# -------------------------------------------------------------------------------------------
	if [ $warp_to_template = 1 ]; then
		cd ${zshim_dir}
		mkdir -p template_space
		cd template_space
		
		echo "Warping the z-shim reference mean into template space"
		# use the existing anat 2 template warp as a starting point
		sct_register_multimodal \
		-i ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz \
		-d ${zshim_dir}/zref_mean.nii.gz \
		-iseg ${sct_dir}/data/PAM50/template/PAM50_cord.nii.gz \
		-dseg ${zshim_dir}/zref_mean_seg.nii.gz \
		-param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,smooth=1,slicewise=1,iter=3 \
		-initwarp ${anat_dir}/warp_template2anat.nii.gz \
		-initwarpinv ${anat_dir}/warp_anat2template.nii.gz \
		-x spline -ofolder ${zshim_dir}/template_space
		
		for mode in zshim_recon no_zshim; do
			echo "Warping the ${mode} EPI into template space"
			sct_apply_transfo -i ${zshim_dir}/${mode}.nii.gz -d ${sct_dir}/data/PAM50/template/PAM50_t2s.nii.gz -w warp_zref_mean2PAM50_t2s.nii.gz -o ${mode}_template_space.nii.gz
		done
	fi
	
done

# ---------------------------------------------------------------------------------------------------------
# 4. Extract the signal in the PAM50 cord mask for both the z-shim and no z-shim template space image
# ---------------------------------------------------------------------------------------------------------
if [ $extract_signal_intensity = 1 ]; then
	# limit to area where all subjects have data, dataset 3 has a range in template space 790 - 859
	mkdir -p ${group_dir}
	cd ${group_dir}
	fslroi $sct_dir/data/PAM50/template/PAM50_cord.nii.gz PAM50_cord_cut.nii.gz 0 -1 0 -1 790 69 0 -1
	
	for subject in {17..41}; do
		printf -v sub "%02d" $subject
		echo "subject: " $sub
		out_sub_dir=$project_dir/derivatives/sub-sspr${sub}
		zshim_dir=$out_sub_dir/func/zshim
		cd ${zshim_dir}/template_space
		for mode in zshim_recon no_zshim; do
			fslroi ${mode}_template_space.nii.gz ${mode}_template_space_cut.nii.gz 0 -1 0 -1 790 69 0 -1
			fslmeants -i ${mode}_template_space_cut.nii.gz -m ${group_dir}/PAM50_cord_cut.nii.gz --showall -o ${mode}_template_space_PAM50_cord_cut_all.txt
		done
	done
fi

# ---------------------------------------------------------------------------------------------------------
# 5. Average both the z-shim and no z-shim template space images for all subjects
# ---------------------------------------------------------------------------------------------------------
if [ $average_images = 1 ]; then

	mkdir -p ${group_dir}
	cd ${group_dir}
	for mode in zshim_recon no_zshim; do
		echo "Gathering all subjects' ${mode} data"
		all_imgs=()
		for subject in {17..41}; do
			printf -v sub "%02d" $subject
			zshim_dir=$project_dir/derivatives/sub-sspr${sub}/func/zshim/template_space
			all_imgs+="$zshim_dir/${mode}_template_space.nii.gz "
		done
		fslmerge -t ${mode}_all_subs.nii.gz ${all_imgs[@]}
		fslmaths ${mode}_all_subs.nii.gz -Tmean ${mode}_mean.nii.gz
	done
fi

