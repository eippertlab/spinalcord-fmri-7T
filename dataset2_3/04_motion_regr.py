# extract higher order motion parameters
# Ulrike Horn
# uhorn@cbs.mpg.de
# 15th Sep 2023

import numpy as np
import pandas as pd
import os
import nibabel as nib

output_path = '/data/pt_02661_raw/Heatpain/derivatives/'

subjects = np.arange(1, 42).tolist()

for isub in subjects:
    sub = 'sub-sspr' + str("{:02d}".format(isub))
    print('Starting with subject: ' + sub)

    # Some subjects have fewer runs
    if (isub == 25) | (isub == 27) | (isub == 29) | (isub == 6):
        num_runs = 2
    elif isub == 31:
        num_runs = 3
    else:
        num_runs = 4

    func_path = os.path.join(output_path, sub, 'func')
    for irun in range(num_runs):
        print('run ' + str(irun+1))
        run_path = func_path + os.sep + 'run-' + str(irun+1)

        # Process x parameters
        x_nifti = nib.load(run_path + os.sep + 'moco_params_x.nii.gz')
        x_data = x_nifti.get_fdata()
        dims = x_data.shape
        # Calculate diff_img, sq_img, and sq_diff_img
        diff_img = np.zeros(dims)
        sq_img = np.zeros(dims)
        sq_diff_img = np.zeros(dims)
        for t in range(1, dims[3]):
            diff_img[:, :, :, t] = x_data[:, :, :, t] - x_data[:, :, :, t - 1]
            sq_img[:, :, :, t] = x_data[:, :, :, t] ** 2
            sq_diff_img[:, :, :, t] = diff_img[:, :, :, t] ** 2
        # Save the modified data to new NIfTI files
        nib.save(nib.Nifti1Image(diff_img, x_nifti.affine), run_path + os.sep + 'moco_params_x_diff.nii.gz')
        nib.save(nib.Nifti1Image(sq_img, x_nifti.affine), run_path + os.sep + 'moco_params_x_squared.nii.gz')
        nib.save(nib.Nifti1Image(sq_diff_img, x_nifti.affine), run_path + os.sep + 'moco_params_x_squared_diff.nii.gz')

        # Process y parameters
        y_nifti = nib.load(run_path + os.sep + 'moco_params_y.nii.gz')
        y_data = y_nifti.get_fdata()
        dims = y_data.shape
        # Calculate diff_img, sq_img, and sq_diff_img
        diff_img = np.zeros(dims)
        sq_img = np.zeros(dims)
        sq_diff_img = np.zeros(dims)
        for t in range(1, dims[3]):
            diff_img[:, :, :, t] = y_data[:, :, :, t] - y_data[:, :, :, t - 1]
            sq_img[:, :, :, t] = y_data[:, :, :, t] ** 2
            sq_diff_img[:, :, :, t] = diff_img[:, :, :, t] ** 2
        # Save the modified data to new NIfTI files
        nib.save(nib.Nifti1Image(diff_img, y_nifti.affine), run_path + os.sep + 'moco_params_y_diff.nii.gz')
        nib.save(nib.Nifti1Image(sq_img, y_nifti.affine), run_path + os.sep + 'moco_params_y_squared.nii.gz')
        nib.save(nib.Nifti1Image(sq_diff_img, y_nifti.affine), run_path + os.sep + 'moco_params_y_squared_diff.nii.gz')
