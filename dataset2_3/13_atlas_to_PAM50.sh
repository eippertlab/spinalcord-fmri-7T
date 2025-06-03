#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
preparation=1
upsampling_for_visuals=0
registration_atlas_ref=0
warp_labels=0
register_neighbors=0
merge=0
downsampling=0
binarize=0
pad=0
cut_middle=0
combine_dorsal=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder where raw and output data is located
code_dir=/data/pt_02661_raw/Heatpain/code  # code folder with subfolder helper_functions
fsl_dir=/afs/cbs.mpg.de/software/fsl/6.0.3/ubuntu-bionic-amd64  # where your FSL is located
sct_dir=/data/u_uhorn_software/sct_6.1  # where your Spinal cord toolbox is located
ants_dir=/data/u_uhorn_software/install/  # where your ANTS installation is located (for zshift only)

# add FSL to the path
FSLDIR=${fsl_dir}
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# add SCT to the path
export PATH="${sct_dir}/bin:${PATH}"

# add ANTS to the path
export PATH="${ants_dir}/bin:${PATH}"

# my input data = atlas images
atlas_path=${project_dir}/derivatives/atlas
atlas_gm=nifti/C6_GM_binary.nii

# the template gm
template_gm=${sct_dir}/data/PAM50/template/PAM50_gm.nii.gz

# where C6 is
middle=834
start=819
end=849

ref_slice_gm=PAM50_gm_upsampled_bin_refslice.nii.gz

cd ${atlas_path}

if [ $preparation = 1 ]; then

	echo "Upsampling and binarizing the PAM50 gray matter"
	sct_resample -i ${template_gm} -mm 0.1x0.1x0.5 -o PAM50_gm_upsampled.nii.gz
	sct_maths -i PAM50_gm_upsampled.nii.gz -bin 0.5 -o PAM50_gm_upsampled_bin.nii.gz

	# Taking a mid C6 slice
	fslroi PAM50_gm_upsampled_bin.nii.gz ${ref_slice_gm} 0 -1 0 -1 ${middle} 1
fi

# upsample some more masks that can be used in figures
if [ $upsampling_for_visuals = 1 ]; then
	template_cord=${sct_dir}/data/PAM50/template/PAM50_cord.nii.gz
	sct_resample -i ${template_cord} -mm 0.1x0.1x0.5 -o PAM50_cord_upsampled.nii.gz
	sct_maths -i PAM50_cord_upsampled.nii.gz -bin 0.5 -o PAM50_cord_upsampled_bin.nii.gz
	
	template_csf=${sct_dir}/data/PAM50/template/PAM50_csf.nii.gz
	sct_resample -i ${template_csf} -mm 0.1x0.1x0.5 -o PAM50_csf_upsampled.nii.gz
	sct_maths -i PAM50_csf_upsampled.nii.gz -bin 0.5 -o PAM50_csf_upsampled_bin.nii.gz
fi


if [ $registration_atlas_ref = 1 ]; then

	echo "Registering the reference slice and the atlas image"

	isct_antsRegistration \
		--dimensionality 2 \
		--output [affine_syn_transform,affine_syn_warped_image.nii.gz] \
		--initial-moving-transform [${ref_slice_gm},${atlas_gm},1] \
		--transform Affine[1] \
		--metric MeanSquares[${ref_slice_gm},${atlas_gm},1,0] \
		--use-histogram-matching 0 \
		--convergence [100x10,1e-6,10] \
		--smoothing-sigmas 1x0mm \
		--shrink-factors 2x1 \
		--transform SyN[0.1,3,0] \
		--metric MeanSquares[${ref_slice_gm},${atlas_gm},1,0] \
		--convergence [100x10,1e-6,10] \
		--smoothing-sigmas 0x0mm \
		--shrink-factors 4x1

	fslcpgeom ${ref_slice_gm} affine_syn_warped_image.nii.gz
fi

if [ $warp_labels = 1 ]; then
	
	echo "Warping all atlas labels to reference slice"
	for label in C6_L1L2_manual C6_L3L4_manual C6_L5L6_manual C6_L7ToL10_manual; do
		echo "${label}"
		
		isct_antsApplyTransforms \
			--dimensionality 2 \
			--input nifti/${label}.nii.gz \
			--reference-image ${ref_slice_gm} \
			--output ${label}_slice${middle}_warped.nii.gz \
			--transform affine_syn_transform1Warp.nii.gz \
			--transform affine_syn_transform0GenericAffine.mat \
			--verbose 1

		fslcpgeom ${ref_slice_gm} ${label}_slice${middle}_warped.nii.gz
	done
fi


if [ $register_neighbors = 1 ]; then

	echo "Registering each slice with the reference slice and warping the labels"
	
	for slice in  $(seq $start $end); do
		
		echo "Slice ${slice}"
		
		if [ $slice -eq $middle ]; then
			continue
		fi
	
		fslroi PAM50_gm_upsampled_bin.nii.gz PAM50_gm_upsampled_bin_${slice}.nii.gz 0 -1 0 -1 ${slice} 1
		
		isct_antsRegistration \
			--dimensionality 2 \
			--output [bsplinesyn_transform_slice${slice}] \
			--transform BSplineSyN[0.1,3,0] \
			--metric MeanSquares[PAM50_gm_upsampled_bin_${slice}.nii.gz,${ref_slice_gm},1,0] \
			--convergence [100x10,1e-6,10] \
			--smoothing-sigmas 0x0mm \
			--shrink-factors 4x1 \
			--verbose 1
		
		for label in C6_L1L2_manual C6_L3L4_manual C6_L5L6_manual C6_L7ToL10_manual; do
			echo "${label}"
			
			isct_antsApplyTransforms \
				--dimensionality 2 \
				--input nifti/${label}.nii.gz \
				--reference-image PAM50_gm_upsampled_bin_${slice}.nii.gz \
				--output ${label}_slice${slice}_warped.nii.gz \
				--transform bsplinesyn_transform_slice${slice}0Warp.nii.gz \
				--transform affine_syn_transform1Warp.nii.gz \
				--transform affine_syn_transform0GenericAffine.mat \
				--verbose 1

			fslcpgeom PAM50_gm_upsampled_bin_${slice} ${label}_slice${slice}_warped.nii.gz
		done

	done
	
fi

if [ $merge = 1 ]; then

	echo "Merge all slices into composite images"
	for label in C6_L1L2_manual C6_L3L4_manual C6_L5L6_manual C6_L7ToL10_manual; do
		echo "${label}"
		fslmerge -z ${label}_Highres.nii.gz $(/bin/ls -1 ${label}_slice*_warped.nii.gz)
	done
	
fi

if [ $downsampling = 1 ]; then

	echo "Downsampling the resulting images"
	for label in C6_L1L2_manual C6_L3L4_manual C6_L5L6_manual C6_L7ToL10_manual; do
		echo "${label}"
		sct_resample -i ${label}_Highres.nii.gz -mm 0.5x0.5x0.5 -o ${label}_PAM50.nii.gz
	done
	
fi

if [ $binarize = 1 ]; then

	echo "Binarizing the resulting images"
	for label in C6_L1L2_manual C6_L3L4_manual C6_L5L6_manual C6_L7ToL10_manual; do
		echo "${label}"
		sct_maths -i ${label}_PAM50.nii.gz -bin 0.5 -o ${label}_PAM50_bin.nii.gz
	done
fi

if [ $pad = 1 ]; then
	
	echo "Creating a blank image"
	
	# Extract header information using fslhd
	header_info=$(fslhd -x $template_gm)
	echo "$header_info" > "header_info_PAM50.xml"
	
	# create empty image
	blank="blank.nii.gz"
	fslcreatehd "header_info_PAM50.xml" $blank
	
	# cut from 0 to start
	zero_below="zero_below.nii.gz"
	fslroi ${blank} ${zero_below} 0 -1 0 -1 0 ${start}
	
	# cut from end to dim3
	zero_above="zero_above.nii.gz"
	dims=$(fslinfo $template_gm | grep -E 'dim[1-3]' | awk '{print $2}')
	dim3=$(echo $dims | awk '{print $3}')
	fslroi ${blank} ${zero_above} 0 -1 0 -1 $(($end+1)) $(($dim3-$end-1))
	
	for label in C6_L1L2 C6_L3L4 C6_L5L6 C6_L7ToL10; do
		echo "${label}"
		
		label_img=${label}_manual_PAM50_bin.nii.gz
		
		fslmerge -z atlas_${label}_PAM50.nii.gz $zero_below $label_img $zero_above
		
	done
	
fi


if [ $cut_middle = 1 ]; then
	echo "Cutting all masks in the middle to have left and right masks"
	for label in C6_L1L2 C6_L3L4 C6_L5L6 C6_L7ToL10; do
		echo "${label}"
		fslmaths atlas_${label}_PAM50.nii.gz -roi 0 70 0 -1 0 -1 0 1 atlas_${label}_right_PAM50.nii.gz
		fslmaths atlas_${label}_PAM50.nii.gz -roi 71 70 0 -1 0 -1 0 1 atlas_${label}_left_PAM50.nii.gz
	done
fi


if [ $combine_dorsal = 1 ]; then
	echo "Combining laminae I to VI to a dorsal mask"
	fslmaths atlas_C6_L1L2_PAM50.nii.gz -add atlas_C6_L3L4_PAM50.nii.gz -add atlas_C6_L5L6_PAM50.nii.gz atlas_C6_dorsal_PAM50.nii.gz
	
	echo "Cutting in the middle to have a left and right dorsal mask"
	fslmaths atlas_C6_dorsal_PAM50.nii.gz -roi 0 69 0 -1 0 -1 0 1 atlas_C6_dorsal_right_PAM50.nii.gz
	fslmaths atlas_C6_dorsal_PAM50.nii.gz -roi 72 69 0 -1 0 -1 0 1 atlas_C6_dorsal_left_PAM50.nii.gz
	
fi

