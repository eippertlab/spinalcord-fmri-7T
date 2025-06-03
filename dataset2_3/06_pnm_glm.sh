#!/bin/bash

# --------------------------------
# Which steps do you want to run?
# --------------------------------
PNM=1
create_onsets=0
create_onsets_even_odd=0
first_level=0
first_level_wo_smooth=0
first_level_even_odd=0

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
	
	if [ $subject -gt 16 ]; then
		data_dir=$project_dir/raw_data/dataset3/sub-sspr${sub}/func
	else
		data_dir=$project_dir/raw_data/dataset2/sub-sspr${sub}/func
	fi
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
	# 1. Put everything in a PNM
	# -------------------------------------------------------------------------------------------
	if [ $PNM -eq 1 ]; then
		echo "Putting together PNM"
		for run in $(seq 1 $num_runs); do
			func_dir=$out_dir/run-${run}
			cd $func_dir
			
			# pnm + csf + moco params + esPCA
			rm -rf pnm
			mkdir pnm
			fslFixText physio.txt pnm/pnm_input.txt
			sr=$(<physio_sr.txt)
			echo "Using a physio file with sampling rate " $sr
			if [ $subject -gt 16 ]; then
				tr=1.123
			else
				tr=1.12
			fi
			popp -i pnm/pnm_input.txt -o pnm/pnm -s ${sr} --tr=${tr} --smoothcard=0.1 --smoothresp=0.1 --resp=3 --cardiac=2 --trigger=1 --pulseox_trigger
			obase=${func_dir}/pnm/pnm
			pnm_evs -i MOCO.nii.gz -c pnm/pnm_card.txt -r pnm/pnm_resp.txt -o pnm/pnm --tr=${tr} --oc=4 --or=4 --multc=2 --multr=2 --csfmask=CSF_only.nii.gz --sliceorder=down --slicedir=z
			ls -1 `imglob -extensions ${obase}ev0*` > pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_x.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_y.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_x_diff.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_y_diff.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_x_squared.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_y_squared.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_x_squared_diff.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/moco_params_y_squared_diff.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/espca_component-1.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/espca_component-2.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/espca_component-3.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/espca_component-4.nii.gz" >> pnm/pnm_evlist.txt
			echo "${func_dir}/espca_component-5.nii.gz" >> pnm/pnm_evlist.txt
			
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 2. Create onset files based on events.tsv
	# -------------------------------------------------------------------------------------------
	if [ $create_onsets -eq 1 ]; then
		echo "Create onsets for GLM"
		for run in $(seq 1 $num_runs); do
			func_dir=$out_dir/run-${run}
			cd $func_dir
			event_file=$data_dir/sub-sspr${sub}_task-longpain_run-${run}_events.tsv
			# read column 1 = onsets without the header
			onsets=$(awk -F'\t' -v col=1 'NR > 1 {print $col}' "$event_file")
			# create output file
			output_file="design_block_regr.txt"
			duration=30
			name=1
			results=()
			results_plus_3=()
			while read -r onset; do
				results+=("$onset")
				tmp=$(echo "$onset + 3.0" | bc -l)
				results_plus_3+=($tmp)
			done <<< $onsets
			paste -d ' ' <(printf "%s\n" "${results[@]}") <(yes "$duration" | head -n "${#results[@]}") <(yes "$name" | head -n "${#results[@]}") > "$output_file"
			# same for onset regressor just shorter duration
			output_file="design_onset_regr.txt"
			duration=3
			paste -d ' ' <(printf "%s\n" "${results[@]}") <(yes "$duration" | head -n "${#results[@]}") <(yes "$name" | head -n "${#results[@]}") > "$output_file"
			# and the regressor you need to model in the onset design for the rest of the 30 s
			# (starting at 3s but with duration 27s)
			output_file="design_onset_regr_rest.txt"
			duration=27
			paste -d ' ' <(printf "%s\n" "${results_plus_3[@]}") <(yes "$duration" | head -n "${#results_plus_3[@]}") <(yes "$name" | head -n "${#results_plus_3[@]}") > "$output_file"
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 3. Create even/odd onset files based on events.tsv
	# -------------------------------------------------------------------------------------------
	if [ $create_onsets_even_odd -eq 1 ]; then
		echo "Create onsets for GLM"
		for run in $(seq 1 $num_runs); do
			func_dir=$out_dir/run-${run}
			cd $func_dir
			event_file=$data_dir/sub-sspr${sub}_task-longpain_run-${run}_events.tsv
			# read column 1 = onsets without the header
			onsets=$(awk -F'\t' -v col=1 'NR > 1 {print $col}' "$event_file")
			# collect all onsets and store new values
			name=1
			even=()
			odd=()
			results_plus_3_even=()
			results_plus_3_odd=()
			counter=1
			
			# as the very first painful stimulus might be always way stronger, we switch even/odd for runs 3 and 4
			if [[ $run -gt 2 ]]; then
				while read -r onset; do
					if [ $((counter % 2)) -eq 0 ]; then
						odd+=("$onset")
						tmp=$(echo "$onset + 3.0" | bc -l)
						results_plus_3_odd+=($tmp)
					else
						even+=("$onset")
						tmp=$(echo "$onset + 3.0" | bc -l)
						results_plus_3_even+=($tmp)
					fi
					counter=$((counter+1))
				done <<< $onsets
			else
				while read -r onset; do
					if [ $((counter % 2)) -eq 0 ]; then
						even+=("$onset")
						tmp=$(echo "$onset + 3.0" | bc -l)
						results_plus_3_even+=($tmp)
					else
						odd+=("$onset")
						tmp=$(echo "$onset + 3.0" | bc -l)
						results_plus_3_odd+=($tmp)
					fi
					counter=$((counter+1))
				done <<< $onsets
			fi
			# block design
			duration=30
			output_file="design_block_regr_even.txt"
			paste -d ' ' <(printf "%s\n" "${even[@]}") <(yes "$duration" | head -n "${#even[@]}") <(yes "$name" | head -n "${#even[@]}") > "$output_file"
			output_file="design_block_regr_odd.txt"
			paste -d ' ' <(printf "%s\n" "${odd[@]}") <(yes "$duration" | head -n "${#odd[@]}") <(yes "$name" | head -n "${#odd[@]}") > "$output_file"
			# onset design (onset)
			duration=3
			output_file="design_onset_regr_even.txt"
			paste -d ' ' <(printf "%s\n" "${even[@]}") <(yes "$duration" | head -n "${#even[@]}") <(yes "$name" | head -n "${#even[@]}") > "$output_file"
			output_file="design_onset_regr_odd.txt"
			paste -d ' ' <(printf "%s\n" "${odd[@]}") <(yes "$duration" | head -n "${#odd[@]}") <(yes "$name" | head -n "${#odd[@]}") > "$output_file"
			# onset design (rest)
			duration=27
			output_file="design_onset_regr_rest_even.txt"
			paste -d ' ' <(printf "%s\n" "${results_plus_3_even[@]}") <(yes "$duration" | head -n "${#results_plus_3_even[@]}") <(yes "$name" | head -n "${#results_plus_3_even[@]}") > "$output_file"
			output_file="design_onset_regr_rest_odd.txt"
			paste -d ' ' <(printf "%s\n" "${results_plus_3_odd[@]}") <(yes "$duration" | head -n "${#results_plus_3_odd[@]}") <(yes "$name" | head -n "${#results_plus_3_odd[@]}") > "$output_file"
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 4. Run first level feat with 2mm susan smoothing
	# -------------------------------------------------------------------------------------------
	if [ $first_level -eq 1 ]; then
		echo "Running feat"
		
		for design in block onset; do
			for run in $(seq 1 $num_runs); do
				echo "Run $run"
				func_dir=$out_dir/run-${run}
				cd $func_dir
				
				# remove existing feat directory
				rm -rf feat_${design}_2mm_susan.feat
				
				cp ${code_dir}/templates/${design}_design.fsf ${design}_design.fsf
				
				# if MOCO outliers were found, use them in first level, if not delete these lines
				if [ ! -f "OUTLIERS.txt" ]; then
					sed -i "s/set fmri(confoundevs) 1/set fmri(confoundevs) 0/g" ${design}_design.fsf
					sed -i "s/# Confound EVs text file for analysis 1//g" ${design}_design.fsf
					sed -i 's/set confoundev_files(1) "MY_OUTPUTDIR\/OUTLIERS.txt"//g' ${design}_design.fsf
				fi
				
				# template used run-1 of subject 01 --> replace with current subject
				sed -i "s/sub-sspr01/sub-sspr${sub}/g" ${design}_design.fsf
				sed -i "s/run-1/run-${run}/g" ${design}_design.fsf
				
				# replace paths
				sed -i "s#MY_OUTPUTDIR#${func_dir}#g" ${design}_design.fsf
				
				# for this subject the number of volumes is different
				if [[ ${subject} == 21 && ${run} == 1 ]];then
					sed -i "s/set fmri(npts) 460/set fmri(npts) 459/g" ${design}_design.fsf
				fi
				
				# sometimes FILM crashed (either turn off FILM or use another FSL version)
				if [[ (${subject} == 13 && ${run} == 1 && ${design} = "block") || (${subject} == 19 && ${run} == 2 && ${design} = "block") || (${subject} == 19 && ${run} == 3 && ${design} = "block") || (${subject} == 23 && ${run} == 4 && ${design} = "block") || (${subject} == 1 && ${run} == 3 && ${design} = "onset") || (${subject} == 2 && ${run} == 4 && ${design} = "onset") || (${subject} == 3 && ${run} == 4 && ${design} = "onset") || (${subject} == 11 && ${run} == 1 && ${design} = "onset") || (${subject} == 19 && ${run} == 4 && ${design} = "onset") || (${subject} == 26 && ${run} == 4 && ${design} = "onset") || (${subject} == 31 && ${run} == 3 && ${design} = "onset") ]];then

					sed -i "s/set fmri(prewhiten_yn) 1/set fmri(prewhiten_yn) 0/g" ${design}_design.fsf

				fi
				
				# run it!
				feat ${design}_design.fsf

			done
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 5. Run first level feat without smoothing
	# -------------------------------------------------------------------------------------------
	if [ $first_level_wo_smooth -eq 1 ]; then
		echo "Running feat"
		
		for design in block onset; do
			for run in $(seq 1 $num_runs); do
				echo "Run $run"
				func_dir=$out_dir/run-${run}
				cd $func_dir
				
				# remove existing feat directory
				rm -rf feat_${design}_no_smooth.feat
				
				cp ${code_dir}/templates/${design}_design.fsf ${design}_design_no_smooth.fsf
				
				# remove the smoothing and change output name
				sed -i "s/set fmri(smooth) 2/set fmri(smooth) 0.0/g" ${design}_design_no_smooth.fsf
				sed -i "s/feat_${design}_2mm_susan/feat_${design}_no_smooth/g" ${design}_design_no_smooth.fsf
				
				# if MOCO outliers were found, use them in first level, if not delete these lines
				if [ ! -f "OUTLIERS.txt" ]; then
					sed -i "s/set fmri(confoundevs) 1/set fmri(confoundevs) 0/g" ${design}_design_no_smooth.fsf
					sed -i "s/# Confound EVs text file for analysis 1//g" ${design}_design_no_smooth.fsf
					sed -i 's/set confoundev_files(1) "MY_OUTPUTDIR\/OUTLIERS.txt"//g' ${design}_design_no_smooth.fsf
				fi
				
				# template used run-1 of subject 01 --> replace with current subject
				sed -i "s/sub-sspr01/sub-sspr${sub}/g" ${design}_design_no_smooth.fsf
				sed -i "s/run-1/run-${run}/g" ${design}_design_no_smooth.fsf
				
				# replace paths
				
				sed -i "s#MY_OUTPUTDIR#${func_dir}#g" ${design}_design_no_smooth.fsf
				
				# for this subject the number of volumes is different
				if [[ ${subject} == 21 && ${run} == 1 ]];then
					sed -i "s/set fmri(npts) 460/set fmri(npts) 459/g" ${design}_design_no_smooth.fsf
				fi
				
				# sometimes FILM crashed (either turn off FILM or use another FSL version)
				if [[ (${subject} == 5 && ${run} == 4 && ${design} = "block") || (${subject} == 23 && ${run} == 1) || (${subject} == 23 && ${run} == 4 && ${design} = "block") || (${subject} == 11 && ${run} == 3 && ${design} = "onset") || (${subject} == 12 && ${run} == 4 && ${design} = "onset") || (${subject} == 39 && ${run} == 3 && ${design} = "onset") ]]; then
					sed -i "s/set fmri(prewhiten_yn) 1/set fmri(prewhiten_yn) 0/g" ${design}_design_no_smooth.fsf
				fi
				
				# run it!
				feat ${design}_design_no_smooth.fsf

			done
		done
	fi
	
	# -------------------------------------------------------------------------------------------
	# 5. Run first level feats with even/odd trials with 2mm susan smoothing
	# -------------------------------------------------------------------------------------------
	if [ $first_level_even_odd -eq 1 ]; then
		for design in block onset; do
			for run in $(seq 1 $num_runs); do
				echo "Running feat run-${run}"
				
				func_dir=$out_dir/run-${run}
				cd $func_dir
				
				# remove existing feat directory
				rm -rf feat_${design}_2mm_susan_even_odd.feat

				cp ${code_dir}/templates/${design}_design_even_odd.fsf ${design}_design_even_odd.fsf
				
				# if MOCO outliers were found use them in first level, if not delete these lines
				if [ ! -f "OUTLIERS.txt" ]; then
					sed -i "s/set fmri(confoundevs) 1/set fmri(confoundevs) 0/g" ${design}_design_even_odd.fsf
					sed -i "s/# Confound EVs text file for analysis 1//g" ${design}_design_even_odd.fsf
					sed -i 's/set confoundev_files(1) "MY_OUTPUTDIR\/OUTLIERS.txt"//g' ${design}_design_even_odd.fsf
				fi
				
				# replace paths
				sed -i "s#MY_OUTPUTDIR#${func_dir}#g" ${design}_design_even_odd.fsf
				sed -i "s#feat_block_and_onset_2mm_susan#feat_block_and_onset_2mm_susan_${type}#g" ${design}_design_even_odd.fsf
				
				# for this subject the number of volumes is different
				if [[ ${subject} == 21 && ${run} == 1 ]];then
					sed -i "s/set fmri(npts) 460/set fmri(npts) 459/g" ${design}_design_even_odd.fsf
				fi
				
				if [[ (${subject} == 17 && ${run} == 2 && ${design} == "block") || (${subject} == 19 && ${run} == 2 && ${design} == "block") || (${subject} == 23 && ${run} == 1) || (${subject} == 31 && ${run} == 3 && ${design} == "block") || (${subject} == 4 && ${run} == 3 && ${design} == "onset") || (${subject} == 11 && ${run} == 3 && ${design} == "onset") || (${subject} == 12 && ${run} == 3 && ${design} == "onset") || (${subject} == 16 && ${run} == 1 && ${design} == "onset") || (${subject} == 23 && ${run} == 4 && ${design} == "onset") || (${subject} == 34 && ${run} == 3 && ${design} == "onset") ]]; then
				
					sed -i "s/set fmri(prewhiten_yn) 1/set fmri(prewhiten_yn) 0/g" ${design}_design_even_odd.fsf
				fi
				
				# run it!
				feat ${design}_design_even_odd.fsf
			done
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

