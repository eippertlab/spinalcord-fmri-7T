# define where to crop anatomical images
# Ulrike Horn
# uhorn@cbs.mpg.de
# 30th Jan 2024

import numpy as np
import pandas as pd
import os
import nibabel as nib


output_path = '/data/pt_02661_raw/Heatpain/derivatives/'

subjects = np.arange(1, 42).tolist()

for isub in subjects:
    sub = 'sub-sspr' + str("{:02d}".format(isub))
    print('Starting with subject: ' + sub)

    os.chdir(output_path + sub + os.sep + 'anat')

    label_file = 'anat_initial_deepseg_labeled.nii.gz'

    # read labels file and cut everything above C1 and below T1
    labels = nib.load(label_file)
    label_img = labels.get_fdata()
    dims = label_img.shape
    [C1_x, C1_y, C1_z] = np.where(label_img == 1)
    max_C1 = max(C1_z)
    [T1_x, T1_y, T1_z] = np.where(label_img == 8)
    min_T1 = min(T1_z)

    # read cord segmentation
    seg_file = 'anat_initial_deepseg.nii.gz'
    seg = nib.load(seg_file)
    seg_img = seg.get_fdata()

    # find within that z-range those dimensions in x and y where the segmentation is still 1
    [seg_x, seg_y, seg_z] = np.where(seg_img == 1)
    to_be_cropped = np.where((seg_z > max_C1) | (seg_z < min_T1))
    seg_x_crop = np.delete(seg_x, to_be_cropped)
    seg_y_crop = np.delete(seg_y, to_be_cropped)
    seg_z_crop = np.delete(seg_z, to_be_cropped)
    min_x = min(seg_x_crop)
    max_x = max(seg_x_crop)
    min_y = min(seg_y_crop)
    max_y = max(seg_y_crop)

    # find middle of image
    mid_x = round(dims[0]/2)
    mid_y = round(dims[1]/2)

    # define where it should be cropped
    crop_x_min = mid_x - 50
    crop_x_max = mid_x + 50
    if isub == 26:
        crop_y_min = mid_y - 50
        crop_y_max = mid_y + 40
    else:
        crop_y_min = mid_y - 60
        crop_y_max = mid_y + 60

    # and check whether the segmentation is always still in there
    assert crop_x_min <= min_x, "Segmentation is too large in x direction"
    assert crop_x_max >= max_x, "Segmentation is too large in x direction"
    assert crop_y_min <= min_y, "Segmentation is too large in y direction"
    assert crop_y_max >= max_y, "Segmentation is too large in y direction"

    # and save the command to be run
    command = 'fslroi anat.nii.gz anat_crop.nii.gz ' +\
              str(crop_x_min) + ' ' + str(crop_x_max - crop_x_min) + ' ' +\
              str(crop_y_min) + ' ' + str(crop_y_max - crop_y_min) + ' ' +\
              str(min_T1) + ' ' + str(max_C1 - min_T1) +\
              '\nfslroi anat_initial_deepseg.nii.gz anat_initial_deepseg_crop.nii.gz ' + \
              str(crop_x_min) + ' ' + str(crop_x_max - crop_x_min) + ' ' + \
              str(crop_y_min) + ' ' + str(crop_y_max - crop_y_min) + ' ' + \
              str(min_T1) + ' ' + str(max_C1 - min_T1)

    with open('crop_command.txt', 'w') as file:
        file.write(command)
