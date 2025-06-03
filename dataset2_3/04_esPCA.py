# slice-wise extra-spinal PCA
# Ulrike Horn
# uhorn@cbs.mpg.de
# with code snippets from Alice Dabbagh
# 6th Sep 2023

import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA
import nibabel as nib

output_path = '/data/pt_02661_raw/Heatpain/derivatives/'

subjects = np.arange(1, 42).tolist()

cut_off = 5

for isub in subjects:

    sub = 'sub-sspr' + str("{:02d}".format(isub))
    print('Starting with subject: ' + sub)

    func_path = os.path.join(output_path, sub, 'func')

    # Some subjects have fewer runs
    if (isub == 25) | (isub == 27) | (isub == 29) | (isub == 6):
        num_runs = 2
    elif isub == 31:
        num_runs = 3
    else:
        num_runs = 4

    for irun in range(num_runs):
        print('run ' + str(irun + 1))
        run_path = func_path + os.sep + 'run-' + str(irun + 1)
        ts = pd.read_csv(run_path + os.sep + 'ts_for_espca.txt', header=None, delim_whitespace=True).transpose()
        ts.rename(columns={0: "x", 1: "y", 2: "z"}, inplace=True)

        imgref = nib.load(run_path + os.sep + 'MOCO.nii.gz')
        [ref_xdim, ref_ydim, num_slices, num_vol] = imgref.shape

        # for each slice make PCA
        slices = np.unique(ts['z'])
        array_list = []
        var_list = []
        for z in slices:
            slice_data = ts.loc[ts['z'] == z].copy()
            slice_data.drop(['x', 'y', 'z'], axis=1, inplace=True)
            pca = PCA(n_components=cut_off)
            principal_comps = pca.fit_transform(slice_data.transpose())
            explained_var = pca.explained_variance_ratio_
            array_list.append(principal_comps)
            var_list.append(explained_var)

        # and save the result as well as explained variances
        result_3d = np.stack(array_list, axis=2)  # has shape n_vols x n_comps x n_slices
        explained_vars = pd.DataFrame(var_list)
        explained_vars.to_csv(run_path + os.sep + 'espca_explained_vars.csv', header=None, index=None)

        # create images that can be used in PNM
        for comp in range(cut_off):
            result = np.ones((ref_xdim, ref_ydim, num_slices, num_vol), dtype=np.float32)
            for sl in range(num_slices):
                for tr in range(num_vol):
                    result[:, :, sl, tr] = result[:, :, sl, tr] * result_3d[tr, comp, sl]

            img = nib.Nifti1Image(result, imgref.affine, imgref.header)
            img.to_filename(run_path + os.sep + "espca_component-" + str(comp+1) + ".nii.gz")
            print('Saved as nifti')
