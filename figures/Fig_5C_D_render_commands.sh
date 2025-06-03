#!/bin/bash

result_dir=/data/pt_02661_raw/Heatpain/derivatives/results/screenshots

mkdir -p $result_dir

design=onset

echo $design
cd $result_dir
mkdir -p $design
cd $design

for dataset in dataset2 dataset3 both; do

	echo $dataset
	cd $result_dir/$design
	mkdir -p $dataset
	cd $dataset
	
	if [[ (${dataset} = "dataset2") ]]; then
		x_coord=76
		y_coord=67
		z_coords=(835 \
		854 \
		791)
	fi
	
	if [[ (${dataset} = "dataset3") ]]; then
		x_coord=76
		y_coord=67
		z_coords=(832 \
		843 \
		795)
	fi
	
	if [[ (${dataset} = "both") ]]; then
		x_coord=76
		y_coord=67
		z_coords=(833 \
		843 \
		791)
	fi
	
	for (( z_coord=0; z_coord<${#z_coords[@]}; z_coord++ ));do 
		# coronal view
		fsleyes render -of coronal_${x_coord}_${y_coord}_${z_coords[$z_coord]}.png -yh -zh -hc -hl -vl ${x_coord} ${y_coord} ${z_coords[$z_coord]} --xcentre  0.03423  0.655 --xzoom 3100 /data/u_uhorn_software/sct_6.1/data/PAM50/template/PAM50_t2.nii.gz -dr 0 4500 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/PAM50_cord/uncorr_t_map_cope1_99.nii.gz --cmap red-yellow -dr 0 4
		
		# coronal view with cursor
		fsleyes render -of coronal_${x_coord}_${y_coord}_${z_coords[$z_coord]}_cursor.png -yh -zh -hl -vl ${x_coord} ${y_coord} ${z_coords[$z_coord]} --xcentre  0.03423  0.655 --xzoom 3100 /data/u_uhorn_software/sct_6.1/data/PAM50/template/PAM50_t2.nii.gz -dr 0 4500 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/PAM50_cord/uncorr_t_map_cope1_99.nii.gz --cmap red-yellow -dr 0 4
		
		# coronal view with C5/C6/C7
		fsleyes render -of coronal_${x_coord}_${y_coord}_${z_coords[$z_coord]}_levels.png -yh -zh -hl -vl ${x_coord} ${y_coord} ${z_coords[$z_coord]} --xcentre  0.03423  0.655 --xzoom 3100 /data/u_uhorn_software/sct_6.1/data/PAM50/template/PAM50_t2.nii.gz -dr 0 4500 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/PAM50_cord/uncorr_t_map_cope1_99.nii.gz --cmap red-yellow -dr 0 4 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/both/C5_cut.nii.gz --overlayType label --lut random --outline --outlineWidth 1 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/both/C6_cut.nii.gz --overlayType label --lut random --outline --outlineWidth 1 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/both/C7_cut.nii.gz --overlayType label --lut random --outline --outlineWidth 1
		
		# axial view
		fsleyes render -of axial_${x_coord}_${y_coord}_${z_coords[$z_coord]}.png -xh -yh -hc -hl -vl ${x_coord} ${y_coord} ${z_coords[$z_coord]} --zcentre -0.00322  0.03382 --zzoom 2000 /data/u_uhorn_software/sct_6.1/data/PAM50/template/PAM50_t2.nii.gz -dr 0 4500 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/PAM50_cord/uncorr_t_map_cope1_99.nii.gz --cmap red-yellow -dr 0 4
		
	done
done


echo $design
cd $result_dir
mkdir -p $design
cd $design

x_coord=76
y_coord=67
z_coords=(832 \
833 \
843)

for dataset in dataset2 dataset3 both; do
	echo $dataset
	cd $result_dir/$design
	mkdir -p $dataset
	cd $dataset
	
	for (( z_coord=0; z_coord<${#z_coords[@]}; z_coord++ ));do 
		# coronal view with dorsal horn left
		fsleyes render -of coronal_${x_coord}_${y_coord}_${z_coords[$z_coord]}_DHL_C6.png -yh -zh -hc -hl -vl ${x_coord} ${y_coord} ${z_coords[$z_coord]} --xcentre  0.03423  0.655 --xzoom 3100 /data/u_uhorn_software/sct_6.1/data/PAM50/template/PAM50_t2.nii.gz -dr 0 4500 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/dh_left_C6/corr_t_map_cope1.nii.gz --cmap red-yellow -dr 0 4 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/dh_left_C6/dh_left_C6_cut.nii.gz --overlayType label --lut harvard-oxford-cortical --outline --outlineWidth 3
		
		# axial view with dorsal horn left
		fsleyes render -of axial_${x_coord}_${y_coord}_${z_coords[$z_coord]}_DHL_C6.png -xh -yh -hc -hl -vl ${x_coord} ${y_coord} ${z_coords[$z_coord]} --zcentre -0.00322  0.03382 --zzoom 2000 /data/u_uhorn_software/sct_6.1/data/PAM50/template/PAM50_t2.nii.gz -dr 0 4500 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/dh_left_C6/corr_t_map_cope1.nii.gz --cmap red-yellow -dr 0 4 /data/pt_02661_raw/Heatpain/derivatives/group_analysis/${dataset}/${design}_design_2mm_susan_vasa/dh_left_C6/dh_left_C6_cut.nii.gz --overlayType label --lut harvard-oxford-cortical --outline --outlineWidth 3

	done
done

