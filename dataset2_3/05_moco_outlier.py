# find outlier volumes with excessive motion
# Ulrike Horn
# uhorn@cbs.mpg.de
# 18th Sep 2023

import numpy as np
import pandas as pd
import os
import nibabel as nib
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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

        # load raw metric data
        metricDVARS = np.loadtxt(run_path + os.sep + 'DVARS.txt')
        metricREFRMS = np.loadtxt(run_path + os.sep + 'REFRMS_edited.txt')
        metricREFRMS = metricREFRMS[1:]

        # calculate own thresholds
        tmpThreshDeDVARS = np.percentile(metricDVARS, 75) + np.subtract(*np.percentile(metricDVARS, [75, 25]))
        tmpThreshZ2DVARS = np.mean(metricDVARS) + 2 * np.std(metricDVARS)
        tmpThreshZ3DVARS = np.mean(metricDVARS) + 3 * np.std(metricDVARS)

        tmpThreshDeREFRMS = np.percentile(metricREFRMS, 75) + np.subtract(*np.percentile(metricREFRMS, [75, 25]))
        tmpThreshZ2REFRMS = np.mean(metricREFRMS) + 2 * np.std(metricREFRMS)
        tmpThreshZ3REFRMS = np.mean(metricREFRMS) + 3 * np.std(metricREFRMS)

        # get number of outliers according to own thresholds
        threshZ2DVARS = np.where(metricDVARS >= tmpThreshZ2DVARS)[0]
        threshZ3DVARS = np.where(metricDVARS >= tmpThreshZ3DVARS)[0]

        threshZ2REFRMS = np.where(metricREFRMS >= tmpThreshZ2REFRMS)[0]
        threshZ3REFRMS = np.where(metricREFRMS >= tmpThreshZ3REFRMS)[0]

        # plot data with thresholds
        fig, axs = plt.subplots(2, 1, figsize=(16, 12))
        axs[0].set_title('Subject {} DVARS run {}'.format(sub, irun+1))
        axs[0].plot(metricDVARS, 'k', label='metric DVARS')
        axs[0].plot(np.arange(len(metricDVARS)), [tmpThreshZ2DVARS] * len(metricDVARS), 'b--', label='Z2 Threshold')
        axs[0].plot(np.arange(len(metricDVARS)), [tmpThreshZ3DVARS] * len(metricDVARS), 'b', label='Z3 Threshold')
        axs[0].plot(np.arange(len(metricDVARS)), [tmpThreshDeDVARS] * len(metricDVARS), 'g', label='De Threshold')
        axs[0].legend(loc='lower right')
        axs[1].set_title('Subject {} REFRMS run {}'.format(sub, irun + 1))
        axs[1].plot(metricREFRMS, 'k', label='metric REFRMS (edited)')
        axs[1].plot(np.arange(len(metricREFRMS)), [tmpThreshZ2REFRMS] * len(metricREFRMS), 'b--', label='Z2 Threshold')
        axs[1].plot(np.arange(len(metricREFRMS)), [tmpThreshZ3REFRMS] * len(metricREFRMS), 'b', label='Z3 Threshold')
        axs[1].plot(np.arange(len(metricREFRMS)), [tmpThreshDeREFRMS] * len(metricREFRMS), 'g', label='De Threshold')
        axs[1].legend(loc='lower right')
        plt.savefig(run_path + os.sep + 'MOCO_outliers.png')
        plt.show()

        # get number of outliers according to default threshold
        numDVARSOutliers = len(threshZ2DVARS)
        indicesDVARSOutliers = threshZ2DVARS
        numREFRMSOutliers = len(threshZ2REFRMS)
        indicesREFRMSOutliers = threshZ2REFRMS

        # Combine the indices of outliers and calculate total and unique outliers
        all_indices = np.concatenate([threshZ2DVARS, threshZ2REFRMS])
        numTotalOutliers = len(all_indices)
        unique_indices = np.unique(all_indices)
        numUniqueOutliers = len(unique_indices)

        # If there are outliers save in file
        if numTotalOutliers == 0 and numUniqueOutliers == 0:
            print('No outliers for subject {} run {}'.format(sub, irun+1))
        else:
            np.savez(run_path + os.sep + 'OUTLIERS.npz',
                     numDVARSOutliers=numDVARSOutliers,
                     indicesDVARSOutliers=indicesDVARSOutliers,
                     numREFRMSOutliers=numREFRMSOutliers,
                     indicesREFRMSOutliers=indicesREFRMSOutliers,
                     indicesOutliers=unique_indices)
            # Create a binary matrix for outliers
            if numTotalOutliers > 0:
                outlier_matrix = np.zeros((len(metricDVARS), len(unique_indices)))
                for i, idx in enumerate(unique_indices):
                    outlier_matrix[idx, i] = 1

                # save the binary matrix as a text file with a space delimiter
                np.savetxt(run_path + os.sep + 'OUTLIERS.txt', outlier_matrix, delimiter=' ')

            del unique_indices
