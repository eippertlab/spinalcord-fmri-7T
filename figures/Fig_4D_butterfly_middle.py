import numpy as np
import os
import nibabel as nib
import matplotlib.pyplot as plt

subjects = np.arange(1, 42).tolist()
# subject 33 does not have an MEGRE scan
subjects.remove(33)

task = 'longpain'

output_path = '/data/pt_02661_raw/Heatpain/derivatives/'

for isub in subjects:
    sub = 'sub-sspr' + str("{:02d}".format(isub))
    print('Starting with subject: ' + sub)

    os.chdir(os.path.join(output_path, sub, 'anat_t2s'))

    # read mask file
    mask_file = 'MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline.nii.gz'
    mask = nib.load(mask_file)
    mask_img = mask.get_fdata()

    upper_right_all = mask_img.copy()
    upper_left_all = mask_img.copy()
    lower_right_all = mask_img.copy()
    lower_left_all = mask_img.copy()
    # Loop through the slices and extract 2D image
    dims = mask_img.shape
    for slice_index in range(mask_img.shape[2]):
        slice_image = mask_img[:, :, slice_index]

        # find indices where mask is 1
        rows, cols = np.where(slice_image == 1)
        min_row, max_row = rows.min(), rows.max()
        min_col, max_col = cols.min(), cols.max()

        # determine middle in left/right direction
        middle_lr = int(min_row + (max_row - min_row) / 2)

        try:
            # then search within this line starting from min_col for the first value to be 1
            first_value = np.where(slice_image[middle_lr, min_col:max_col])[0][0]

            # this is the point where the bridge between the cord's "wings" is = the middle between upper and lower
            middle_ul = min_col + first_value
        # it is possible that this fails if the left and right parts are not connected
        except:
            # then go to the left and the right and track where the ones are
            try:
                left_value = np.where(slice_image[middle_lr - 1, min_col:max_col])[0][0]
            except:
                left_value = np.nan
            try:
                right_value = np.where(slice_image[middle_lr + 1, min_col:max_col])[0][0]
            except:
                right_value = np.nan

            # unlikely but possible: three voxels in a row are 0 --> go even further left and right
            if np.isnan(left_value) and np.isnan(right_value):
                try:
                    left_value = np.where(slice_image[middle_lr - 2, min_col:max_col])[0][0]
                except:
                    left_value = np.nan
                try:
                    right_value = np.where(slice_image[middle_lr + 2, min_col:max_col])[0][0]
                except:
                    right_value = np.nan
            # then take the minimum of both (more dorsal) to include enough in the ventral part
            first_value = int(np.nanmin([left_value, right_value]))
            middle_ul = min_col + first_value

        # now cut the image into 4 parts leaving the middle line between left and right out
        upper_right = slice_image.copy()
        upper_right[min_row:middle_lr + 1, :] = 0
        upper_right[:, min_col:middle_ul] = 0
        upper_right_all[:, :, slice_index] = upper_right

        upper_left = slice_image.copy()
        upper_left[middle_lr:max_row + 1, :] = 0
        upper_left[:, min_col:middle_ul] = 0
        upper_left_all[:, :, slice_index] = upper_left

        lower_right = slice_image.copy()
        lower_right[min_row:middle_lr + 1, :] = 0
        lower_right[:, middle_ul:max_col + 1] = 0
        lower_right_all[:, :, slice_index] = lower_right

        lower_left = slice_image.copy()
        lower_left[middle_lr:max_row + 1, :] = 0
        lower_left[:, middle_ul:max_col + 1] = 0
        lower_left_all[:, :, slice_index] = lower_left

        # display the slice with middle points
        # background_color = 'none'
        # color_upper_right = [1, 0, 0, 0.5]  # red
        # color_upper_left = [0, 0, 1, 0.5]  # blue
        # color_lower_right = [0, 1, 0, 0.5]  # green
        # color_lower_left = [1, 1, 0, 0.5]  # yellow
        # cmap_upper_right = plt.cm.colors.ListedColormap([background_color, color_upper_right])
        # cmap_upper_left = plt.cm.colors.ListedColormap([background_color, color_upper_left])
        # cmap_lower_right = plt.cm.colors.ListedColormap([background_color, color_lower_right])
        # cmap_lower_left = plt.cm.colors.ListedColormap([background_color, color_lower_left])
        # plt.imshow(slice_image, cmap='gray')
        # plt.imshow(upper_right, cmap=cmap_upper_right)
        # plt.imshow(upper_left, cmap=cmap_upper_left)
        # plt.imshow(lower_right, cmap=cmap_lower_right)
        # plt.imshow(lower_left, cmap=cmap_lower_left)
        # plt.show()

    # Save in nifti
    # here I needed to save it with left and right reversed as the images were displayed in radiological convention
    img = nib.Nifti1Image(upper_left_all, mask.affine, mask.header)
    img.to_filename(os.path.join(output_path, sub, 'anat_t2s',
                                 'MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline_right_ventral.nii.gz'))
    img = nib.Nifti1Image(lower_left_all, mask.affine, mask.header)
    img.to_filename(os.path.join(output_path, sub, 'anat_t2s',
                                 'MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline_right_dorsal.nii.gz'))
    img = nib.Nifti1Image(upper_right_all, mask.affine, mask.header)
    img.to_filename(os.path.join(output_path, sub, 'anat_t2s',
                                 'MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline_left_ventral.nii.gz'))
    img = nib.Nifti1Image(lower_right_all, mask.affine, mask.header)
    img.to_filename(os.path.join(output_path, sub, 'anat_t2s',
                                 'MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline_left_dorsal.nii.gz'))
    print('Saved all as niftis')
