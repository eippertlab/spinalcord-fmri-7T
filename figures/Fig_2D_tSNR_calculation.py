#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 13:31:17 2025

@author: uhorn
"""

import os
import numpy as np
import nibabel as nib

dataset1_path = '/data/pt_02616/data/derivatives/group_analysis/'
dataset2_path = '/data/pt_02661_raw/Heatpain/derivatives/group_analysis'
sct_path = '/data/u_uhorn_software/sct_6.1'

raw_tsnr1_file = os.path.join(dataset1_path, 
                         'Raw_tsnr_all_subs_sosGREp3_mean.nii.gz')

raw_tsnr2_file = os.path.join(dataset2_path, 
                         'Raw_tsnr_100_all_subs_dataset1_mean.nii.gz')

raw_tsnr1 = nib.load(raw_tsnr1_file)
raw_tsnr1_img = raw_tsnr1.get_fdata()

raw_tsnr2 = nib.load(raw_tsnr2_file)
raw_tsnr2_img = raw_tsnr2.get_fdata()

cord_mask = os.path.join(sct_path, 'data', 'PAM50', 
                         'template', 'PAM50_cord.nii.gz')
mask = nib.load(cord_mask)
mask_img = mask.get_fdata()

# cutting dataset1
mask_img[:, :, 855:] = 0
mask_img[:, :, :785] = 0

raw_mean1 = np.mean(raw_tsnr1_img[mask_img == 1])
raw_mean2 = np.mean(raw_tsnr2_img[mask_img == 1])

print("The tSNR in the cord in dataset 1 is {}".format(round(raw_mean1, 3)))
print("The tSNR in the cord in dataset 2 is {}".format(round(raw_mean2, 3)))
