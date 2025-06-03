#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
higherlevel_subject=1
higherlevel_subject_no_smooth=0
higherlevel_subject_even_odd=0
higherlevel_subject_vasa=0
higherlevel_subject_no_smooth_vasa=0
higherlevel_subject_even_odd_vasa=0

# --------------------------------
# Define all directories
# --------------------------------
project_dir=/data/pt_02661_raw/Heatpain  # general project folder
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

# either you use this loop
# for subject in {1..41}; do

# or you use this function to process each subject and parallelize in code at the bottom
process_subject() {
  	sub=$1

	printf -v sub "%02d" $subject
	echo "subject: " $sub
	
	out_dir=$project_dir/derivatives/sub-sspr${sub}/func
	
	# when some runs were not recorded
	if [[ $subject == 25 || $subject == 27 || $subject == 29 || $subject == 6 ]]; then
		num_runs=2
	elif [[ $subject == 31 ]]; then
		num_runs=3
	else
		num_runs=4
	fi
	
	# -------------------------------------------------------------------------------------------
	# 1. Run second level feat in subject space (2mm susan smoothing)
	# -------------------------------------------------------------------------------------------
	if [ $higherlevel_subject = 1 ]; then
		echo "Running subject level analysis"
		
		seclevel_dir=$out_dir/second_level
		cd $seclevel_dir
		
		for design in block onset; do
			# Replace the feat-generated 'reg' subdirectory with the dummy identity one
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/feat_${design}_2mm_susan.feat
				cd $feat_dir
				cp $FSLDIR/etc/flirtsch/ident.mat reg/example_func2standard.mat
				cp $seclevel_dir/MOCO_MEAN_all_runs.nii.gz reg/standard.nii.gz
			done
			
			# copy the template feat into the folder, fill it with specific info and run it
			cd $seclevel_dir
			cp ${code_dir}/templates/${design}_design_second_level.fsf ${design}_design_second_level.fsf
			# to give path names you need another delimiter than /
			sed -i "s#MY_OUTPUTDIR#${seclevel_dir}#g" ${design}_design_second_level.fsf
			sed -i "s#MY_FSLDIR#${fsl_dir}#g" ${design}_design_second_level.fsf
			# as some subjects have less runs, we need to adapt the file for this
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/feat_${design}_2mm_susan.feat
				escaped_feat_dir=$(echo "$feat_dir" | sed 's/\//\\\//g')
				text_to_insert="# 4D AVW data or FEAT directory (${run})\nset feat_files(${run}) \"${escaped_feat_dir}\"\n"
				sed -i "/MY_FEATDIRS/ i\\
				${text_to_insert}" ${design}_design_second_level.fsf
				
				text_to_insert="# Higher-level EV value for EV 1 and input ${run}\nset fmri(evg${run}.1) 1\n"
				sed -i "/MY_EV_VALUES/ i\\
				${text_to_insert}" ${design}_design_second_level.fsf
				
				text_to_insert="# Group membership for input ${run}\nset fmri(groupmem.${run}) 1\n"
				sed -i "/MY_MEMBERSHIP/ i\\
				${text_to_insert}" ${design}_design_second_level.fsf
				
				sed -i "s/set fmri(npts) 4/set fmri(npts) ${num_runs}/g" ${design}_design_second_level.fsf
				sed -i "s/set fmri(multiple) 4/set fmri(multiple) ${num_runs}/g" ${design}_design_second_level.fsf
				
			done
			sed -i "/MY_FEATDIRS/d" ${design}_design_second_level.fsf
			sed -i "/MY_EV_VALUES/d" ${design}_design_second_level.fsf
			sed -i "/MY_MEMBERSHIP/d" ${design}_design_second_level.fsf
			
			# replace subject name
			sed -i "s/sub-sspr01/sub-sspr${sub}/g" ${design}_design_second_level.fsf
			
			# run it!
			feat ${design}_design_second_level.fsf
		done

	fi
	
	# -------------------------------------------------------------------------------------------
	# 1. Run second level feat in subject space (no smoothing)
	# -------------------------------------------------------------------------------------------
	if [ $higherlevel_subject_no_smooth = 1 ]; then
		echo "Running subject level analysis"
		
		seclevel_dir=$out_dir/second_level
		cd $seclevel_dir
		
		for design in block onset; do
			# Replace the feat-generated 'reg' subdirectory with the dummy identity one
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/feat_${design}_no_smooth.feat
				cd $feat_dir
				cp $FSLDIR/etc/flirtsch/ident.mat reg/example_func2standard.mat
				cp $seclevel_dir/MOCO_MEAN_all_runs.nii.gz reg/standard.nii.gz
			done
			
			# copy the template feat into the folder, fill it with specific info and run it
			cd $seclevel_dir
			cp ${code_dir}/templates/${design}_design_second_level.fsf ${design}_design_second_level_no_smooth.fsf
			
			# change output name
			sed -i "s/${design}_design_2mm_susan/${design}_design_no_smooth/g" ${design}_design_no_smooth.fsf
			
			# to give path names you need another delimiter than /
			sed -i "s#MY_OUTPUTDIR#${seclevel_dir}#g" ${design}_design_second_level_no_smooth.fsf
			sed -i "s#MY_FSLDIR#${fsl_dir}#g" ${design}_design_second_level_no_smooth.fsf
			# as some subjects have less runs, we need to adapt the file for this
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/feat_${design}_no_smooth.feat
				escaped_feat_dir=$(echo "$feat_dir" | sed 's/\//\\\//g')
				text_to_insert="# 4D AVW data or FEAT directory (${run})\nset feat_files(${run}) \"${escaped_feat_dir}\"\n"
				sed -i "/MY_FEATDIRS/ i\\
				${text_to_insert}" ${design}_design_second_level_no_smooth.fsf
				
				text_to_insert="# Higher-level EV value for EV 1 and input ${run}\nset fmri(evg${run}.1) 1\n"
				sed -i "/MY_EV_VALUES/ i\\
				${text_to_insert}" ${design}_design_second_level_no_smooth.fsf
				
				text_to_insert="# Group membership for input ${run}\nset fmri(groupmem.${run}) 1\n"
				sed -i "/MY_MEMBERSHIP/ i\\
				${text_to_insert}" ${design}_design_second_level_no_smooth.fsf
				
				sed -i "s/set fmri(npts) 4/set fmri(npts) ${num_runs}/g" ${design}_design_second_level_no_smooth.fsf
				sed -i "s/set fmri(multiple) 4/set fmri(multiple) ${num_runs}/g" ${design}_design_second_level_no_smooth.fsf
				
			done
			sed -i "/MY_FEATDIRS/d" ${design}_design_second_level_no_smooth.fsf
			sed -i "/MY_EV_VALUES/d" ${design}_design_second_level_no_smooth.fsf
			sed -i "/MY_MEMBERSHIP/d" ${design}_design_second_level_no_smooth.fsf
			
			# replace subject name
			sed -i "s/sub-sspr01/sub-sspr${sub}/g" ${design}_design_second_level_no_smooth.fsf
			
			# run it!
			feat ${design}_design_second_level_no_smooth.fsf
		done

	fi
	
	# -------------------------------------------------------------------------------------------
	# 2. Run second level feat in subject space (2 mm susan, even & odd)
	# -------------------------------------------------------------------------------------------
	if [ $higherlevel_subject_even_odd = 1 ]; then
		echo "Running subject level analysis"
		
		seclevel_dir=$out_dir/second_level
		cd $seclevel_dir
		
		for design in block onset; do
			# Replace the feat-generated 'reg' subdirectory with the dummy identity one
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/feat_${design}_2mm_susan_even_odd.feat
				cd $feat_dir
				cp $FSLDIR/etc/flirtsch/ident.mat reg/example_func2standard.mat
				cp $seclevel_dir/MOCO_MEAN_all_runs.nii.gz reg/standard.nii.gz
			done
			
			# copy the template feat into the folder, fill it with specific info and run it
			cd $seclevel_dir
			cp ${code_dir}/templates/${design}_design_even_odd_second_level.fsf ${design}_design_even_odd_second_level.fsf
			# to give path names you need another delimiter than /
			sed -i "s#MY_OUTPUTDIR#${seclevel_dir}#g" ${design}_design_even_odd_second_level.fsf
			sed -i "s#MY_FSLDIR#${fsl_dir}#g" ${design}_design_even_odd_second_level.fsf
			# as some subjects have less runs, we need to adapt the file for this
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/feat_${design}_2mm_susan_even_odd.feat
				escaped_feat_dir=$(echo "$feat_dir" | sed 's/\//\\\//g')
				text_to_insert="# 4D AVW data or FEAT directory (${run})\nset feat_files(${run}) \"${escaped_feat_dir}\"\n"
				sed -i "/MY_FEATDIRS/ i\\
				${text_to_insert}" ${design}_design_even_odd_second_level.fsf
				
				text_to_insert="# Higher-level EV value for EV 1 and input ${run}\nset fmri(evg${run}.1) 1\n"
				sed -i "/MY_EV_VALUES/ i\\
				${text_to_insert}" ${design}_design_even_odd_second_level.fsf
				
				text_to_insert="# Group membership for input ${run}\nset fmri(groupmem.${run}) 1\n"
				sed -i "/MY_MEMBERSHIP/ i\\
				${text_to_insert}" ${design}_design_even_odd_second_level.fsf
				
				sed -i "s/set fmri(npts) 4/set fmri(npts) ${num_runs}/g" ${design}_design_even_odd_second_level.fsf
				sed -i "s/set fmri(multiple) 4/set fmri(multiple) ${num_runs}/g" ${design}_design_even_odd_second_level.fsf
				
			done
			sed -i "/MY_FEATDIRS/d" ${design}_design_even_odd_second_level.fsf
			sed -i "/MY_EV_VALUES/d" ${design}_design_even_odd_second_level.fsf
			sed -i "/MY_MEMBERSHIP/d" ${design}_design_even_odd_second_level.fsf
			
			# replace subject name
			sed -i "s/sub-sspr01/sub-sspr${sub}/g" ${design}_design_even_odd_second_level.fsf
			
			# run it!
			feat ${design}_design_even_odd_second_level.fsf
		done

	fi
	
	# -------------------------------------------------------------------------------------------
	# 3. Run second level feat in subject space (2 mm susan, VasA)
	# -------------------------------------------------------------------------------------------
	if [ $higherlevel_subject_vasa = 1 ]; then
		echo "Running subject level analysis"
		
		seclevel_dir=$out_dir/second_level
		cd $seclevel_dir
		
		for design in block onset; do
			# Replace the feat-generated 'reg' subdirectory with the dummy identity one
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/vasa_feat_${design}_2mm_susan.feat
				cd $feat_dir
				cp $FSLDIR/etc/flirtsch/ident.mat reg/example_func2standard.mat
				cp $seclevel_dir/MOCO_MEAN_all_runs.nii.gz reg/standard.nii.gz
			done
			
			# copy the template feat into the folder, fill it with specific info and run it
			cd $seclevel_dir
			cp ${code_dir}/templates/vasa_${design}_design_second_level.fsf vasa_${design}_design_second_level.fsf
			
			# to give path names you need another delimiter than /
			sed -i "s#MY_OUTPUTDIR#${seclevel_dir}#g" vasa_${design}_design_second_level.fsf
			sed -i "s#MY_FSLDIR#${fsl_dir}#g" vasa_${design}_design_second_level.fsf
			# as some subjects have less runs, we need to adapt the file for this
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/vasa_feat_${design}_2mm_susan.feat
				escaped_feat_dir=$(echo "$feat_dir" | sed 's/\//\\\//g')
				text_to_insert="# 4D AVW data or FEAT directory (${run})\nset feat_files(${run}) \"${escaped_feat_dir}\"\n"
				sed -i "/MY_FEATDIRS/ i\\
				${text_to_insert}" vasa_${design}_design_second_level.fsf
				
				text_to_insert="# Higher-level EV value for EV 1 and input ${run}\nset fmri(evg${run}.1) 1\n"
				sed -i "/MY_EV_VALUES/ i\\
				${text_to_insert}" vasa_${design}_design_second_level.fsf
				
				text_to_insert="# Group membership for input ${run}\nset fmri(groupmem.${run}) 1\n"
				sed -i "/MY_MEMBERSHIP/ i\\
				${text_to_insert}" vasa_${design}_design_second_level.fsf
				
				sed -i "s/set fmri(npts) 4/set fmri(npts) ${num_runs}/g" vasa_${design}_design_second_level.fsf
				sed -i "s/set fmri(multiple) 4/set fmri(multiple) ${num_runs}/g" vasa_${design}_design_second_level.fsf
				
			done
			sed -i "/MY_FEATDIRS/d" vasa_${design}_design_second_level.fsf
			sed -i "/MY_EV_VALUES/d" vasa_${design}_design_second_level.fsf
			sed -i "/MY_MEMBERSHIP/d" vasa_${design}_design_second_level.fsf
			
			# replace subject name
			sed -i "s/sub-sspr01/sub-sspr${sub}/g" vasa_${design}_design_second_level.fsf
			
			# run it!
			feat vasa_${design}_design_second_level.fsf
		done

	fi
	
	# -------------------------------------------------------------------------------------------
	# 3. Run second level feat in subject space (no smooth, VasA)
	# -------------------------------------------------------------------------------------------
	if [ $higherlevel_subject_no_smooth_vasa = 1 ]; then
		echo "Running subject level analysis"
		
		seclevel_dir=$out_dir/second_level
		cd $seclevel_dir
		
		for design in block; do
			# Replace the feat-generated 'reg' subdirectory with the dummy identity one
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/vasa_feat_${design}_no_smooth.feat
				cd $feat_dir
				cp $FSLDIR/etc/flirtsch/ident.mat reg/example_func2standard.mat
				cp $seclevel_dir/MOCO_MEAN_all_runs.nii.gz reg/standard.nii.gz
			done
			
			# copy the template feat into the folder, fill it with specific info and run it
			cd $seclevel_dir
			cp ${code_dir}/templates/vasa_${design}_design_second_level.fsf vasa_${design}_design_second_level_no_smooth.fsf
			
			# change output name
			sed -i "s/${design}_design_2mm_susan_vasa/${design}_design_no_smooth_vasa/g" vasa_${design}_design_second_level_no_smooth.fsf
			
			# to give path names you need another delimiter than /
			sed -i "s#MY_OUTPUTDIR#${seclevel_dir}#g" vasa_${design}_design_second_level_no_smooth.fsf
			sed -i "s#MY_FSLDIR#${fsl_dir}#g" vasa_${design}_design_second_level_no_smooth.fsf
			# as some subjects have less runs, we need to adapt the file for this
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/vasa_feat_${design}_no_smooth.feat
				escaped_feat_dir=$(echo "$feat_dir" | sed 's/\//\\\//g')
				text_to_insert="# 4D AVW data or FEAT directory (${run})\nset feat_files(${run}) \"${escaped_feat_dir}\"\n"
				sed -i "/MY_FEATDIRS/ i\\
				${text_to_insert}" vasa_${design}_design_second_level_no_smooth.fsf
				
				text_to_insert="# Higher-level EV value for EV 1 and input ${run}\nset fmri(evg${run}.1) 1\n"
				sed -i "/MY_EV_VALUES/ i\\
				${text_to_insert}" vasa_${design}_design_second_level_no_smooth.fsf
				
				text_to_insert="# Group membership for input ${run}\nset fmri(groupmem.${run}) 1\n"
				sed -i "/MY_MEMBERSHIP/ i\\
				${text_to_insert}" vasa_${design}_design_second_level_no_smooth.fsf
				
				sed -i "s/set fmri(npts) 4/set fmri(npts) ${num_runs}/g" vasa_${design}_design_second_level_no_smooth.fsf
				sed -i "s/set fmri(multiple) 4/set fmri(multiple) ${num_runs}/g" vasa_${design}_design_second_level_no_smooth.fsf
				
			done
			sed -i "/MY_FEATDIRS/d" vasa_${design}_design_second_level_no_smooth.fsf
			sed -i "/MY_EV_VALUES/d" vasa_${design}_design_second_level_no_smooth.fsf
			sed -i "/MY_MEMBERSHIP/d" vasa_${design}_design_second_level_no_smooth.fsf
			
			# replace subject name
			sed -i "s/sub-sspr01/sub-sspr${sub}/g" vasa_${design}_design_second_level_no_smooth.fsf
			
			# run it!
			feat vasa_${design}_design_second_level_no_smooth.fsf
		done

	fi
	
	
	# -------------------------------------------------------------------------------------------
	# 4. Run second level feat in subject space (2 mm susan, VasA, even odd)
	# -------------------------------------------------------------------------------------------
	if [ $higherlevel_subject_even_odd_vasa = 1 ]; then
		echo "Running subject level analysis"
		
		seclevel_dir=$out_dir/second_level
		cd $seclevel_dir
		
		for design in block onset; do
			# Replace the feat-generated 'reg' subdirectory with the dummy identity one
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/vasa_feat_${design}_2mm_susan_even_odd.feat
				cd $feat_dir
				cp $FSLDIR/etc/flirtsch/ident.mat reg/example_func2standard.mat
				cp $seclevel_dir/MOCO_MEAN_all_runs.nii.gz reg/standard.nii.gz
			done
			
			# copy the template feat into the folder, fill it with specific info and run it
			cd $seclevel_dir
			cp ${code_dir}/templates/vasa_${design}_design_second_level_even_odd.fsf vasa_${design}_design_second_level_even_odd.fsf
			# to give path names you need another delimiter than /
			sed -i "s#MY_OUTPUTDIR#${seclevel_dir}#g" vasa_${design}_design_second_level_even_odd.fsf
			sed -i "s#MY_FSLDIR#${fsl_dir}#g" vasa_${design}_design_second_level_even_odd.fsf
			# as some subjects have less runs, we need to adapt the file for this
			for run in $(seq 1 $num_runs); do
				feat_dir=$out_dir/run-${run}/vasa_feat_${design}_2mm_susan_even_odd.feat
				escaped_feat_dir=$(echo "$feat_dir" | sed 's/\//\\\//g')
				text_to_insert="# 4D AVW data or FEAT directory (${run})\nset feat_files(${run}) \"${escaped_feat_dir}\"\n"
				sed -i "/MY_FEATDIRS/ i\\
				${text_to_insert}" vasa_${design}_design_second_level_even_odd.fsf
				
				text_to_insert="# Higher-level EV value for EV 1 and input ${run}\nset fmri(evg${run}.1) 1\n"
				sed -i "/MY_EV_VALUES/ i\\
				${text_to_insert}" vasa_${design}_design_second_level_even_odd.fsf
				
				text_to_insert="# Group membership for input ${run}\nset fmri(groupmem.${run}) 1\n"
				sed -i "/MY_MEMBERSHIP/ i\\
				${text_to_insert}" vasa_${design}_design_second_level_even_odd.fsf
				
				sed -i "s/set fmri(npts) 4/set fmri(npts) ${num_runs}/g" vasa_${design}_design_second_level_even_odd.fsf
				sed -i "s/set fmri(multiple) 4/set fmri(multiple) ${num_runs}/g" vasa_${design}_design_second_level_even_odd.fsf
				
			done
			sed -i "/MY_FEATDIRS/d" vasa_${design}_design_second_level_even_odd.fsf
			sed -i "/MY_EV_VALUES/d" vasa_${design}_design_second_level_even_odd.fsf
			sed -i "/MY_MEMBERSHIP/d" vasa_${design}_design_second_level_even_odd.fsf
			
			# replace subject name
			sed -i "s/sub-sspr01/sub-sspr${sub}/g" vasa_${design}_design_second_level_even_odd.fsf
			
			# run it!
			feat vasa_${design}_design_second_level_even_odd.fsf
		done

	fi
	
# use done to end normal loop
# done
# or this bracket to end process
}

# Loop across subjects for data preparation in parallel
for subject in {1..41}; do
	printf -v sub "%02d" $subject
	process_subject $sub & 
done

# Wait for all background processes to finish
wait
echo "All parallel processes completed."

