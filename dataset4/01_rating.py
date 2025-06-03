# Analysis of continuous rating data

import numpy as np
import pandas as pd
import os
import glob
import mne
from scipy.interpolate import interp1d
from scipy import signal, stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import seaborn as sns
import matplotlib.patches as patches

main_path = '/data/pt_02661_raw/Heatpain_behav/data'
raw_path = os.path.join(main_path, 'raw_data')
output_path = os.path.join(main_path, 'derivatives')
result_path = os.path.join(output_path, 'results')

color1 = [0/255, 90/255, 50/255]

if not os.path.exists(output_path):
    os.makedirs(output_path)

if not os.path.exists(result_path):
    os.makedirs(result_path)

interval = [-5, 35]  # interval for cutting and plotting the data
task = 'continrating'
desired_sr = 36  # 36 Hz

subjects = np.arange(1, 21)
excluded = []

preprocess = False
plot = False
group_plot = True

new_rc_params = {"font.family": 'Arial', "font.size": 12, "font.serif": [],
                 "svg.fonttype": 'none'}
mpl.rcParams.update(new_rc_params)

for isub in subjects:
    sub_id = 'sub-cr' + '{:02d}'.format(isub)
    print(sub_id)
    physio_files = glob.glob(os.path.join(raw_path, sub_id, 'beh', sub_id + '_task-' + task + '*rating.csv'))
    physio_files = np.array(sorted(physio_files))

    sub_out_path = os.path.join(output_path, sub_id, 'rating')
    if not os.path.exists(sub_out_path):
        os.makedirs(sub_out_path)

    if preprocess:
        samples_before_all = []
        samples_after_all = []
        for i, p in enumerate(physio_files):
            df = pd.read_csv(p, header=0, index_col=False)

            # the problem is there is no fixed sampling rate,
            # but I want the epochs to be the same size
            # find onsets and track how long the epochs should be
            # maybe we even need to interpolate? --> think not: SR is consistently between 36 and 37
            for trial in np.unique(df['Trial']):
                tmp = df[df['Trial'] == trial].reset_index()
                onset = tmp.Time[tmp['Heat'].values.argmax()]
                start_epoch = onset + interval[0]
                samples_before = ((tmp['Time'] >= start_epoch) & (tmp['Time'] <= onset)).sum()
                samples_before_all.append(samples_before)
                end_epoch = onset + interval[1]
                samples_after = ((tmp['Time'] >= onset) & (tmp['Time'] <= end_epoch)).sum()
                samples_after_all.append(samples_after)
        print('You need min {} and max {} samples before the onset.'.format(np.min(samples_before_all),
                                                                            np.max(samples_before_all)))

        print('The sampling rate goes from {} to {} Hz.'.format(np.min(samples_before_all)/(-interval[0]),
                                                                np.max(samples_before_all)/(-interval[0])))
        print('You need min {} and max {} samples after the onset.'.format(np.min(samples_after_all),
                                                                           np.max(samples_after_all)))
        print('The sampling rate goes from {} to {} Hz.'.format(np.min(samples_after_all) / interval[1],
                                                                np.max(samples_after_all) / interval[1]))

        sampling_rate = int(np.round(np.mean(samples_after_all)/interval[1]))
        print('Sampling rate: {}'.format(sampling_rate))

        

        epochs_all = []
        trial_infos_all = pd.DataFrame()
        for i, p in enumerate(physio_files):
            df = pd.read_csv(p, header=0, index_col=False)

            for trial in np.unique(df['Trial']):
                tmp = df[df['Trial'] == trial].reset_index()
                interpolator = interp1d(tmp['Time'], tmp['Rating'], kind='linear',
                                        bounds_error=False, fill_value='extrapolate')
                new_index = np.arange(0.0, 60.0, 1 / desired_sr)
                interpolated_values = interpolator(new_index)
                
                # this is just to confirm things:
                interpolator = interp1d(tmp['Time'], tmp['Heat'], kind='linear',
                                        bounds_error=False, fill_value='extrapolate')
                interpolated_heat = interpolator(new_index)
                onset = interpolated_heat.argmax()
                #print(onset)
                # the onset in the original is at 15s
                epoch = np.array(interpolated_values[15*desired_sr + interval[0] * desired_sr:15*desired_sr + interval[1] * desired_sr])
                
                # this subject forgot to rate for one trial --> exclude this:
                if (sub_id == 'sub-cr18') & (i == 3) & (trial == 4):
                    print('One trial was excluded, no rating was given')
                else:
                    epochs_all.append(epoch)

            # load event file
            event_file = os.path.join(raw_path, sub_id, 'beh', sub_id + '_task-' + task + '_run-' + str(i + 1) +
                                      '_events.csv')
            trial_infos = pd.read_csv(event_file, index_col=None)
            
            if (sub_id == 'sub-cr18') & (i == 3):
                trial_infos.drop(3, axis=0, inplace=True)
                
            # stack all meta data tables
            trial_infos_all = pd.concat([trial_infos_all, trial_infos])

        epochs_all = np.vstack(epochs_all)
        assert epochs_all.shape[0] == len(trial_infos_all), "Epochs and trial info are not the same size"

        save_name = sub_out_path + os.sep + 'rating_' + task + '_epo'
        np.save(save_name, epochs_all, allow_pickle=True)
        save_name = sub_out_path + os.sep + 'rating_' + task + '_all_events.csv'
        trial_infos_all.to_csv(save_name)

    if plot:
        save_name = sub_out_path + os.sep + 'rating_' + task + '_epo.npy'
        epochs = np.load(save_name, allow_pickle=True)
        
        save_name = sub_out_path + os.sep + 'rating_' + task + '_all_events.csv'
        trial_infos = pd.read_csv(save_name, index_col=False, header=0)
        #trial_infos = trial_infos_all
        
        x = np.linspace(interval[0], interval[1], epochs.shape[1])
        intensity_data = epochs[trial_infos.rating_type == 'intensity']
        unpleasantness_data = epochs[trial_infos.rating_type == 'unpleasantness']
        data = pd.DataFrame({'Time': np.tile(x, intensity_data.shape[0] + unpleasantness_data.shape[0]),
                             'Rating': np.concatenate([intensity_data.flatten(), unpleasantness_data.flatten()]),
                             'Type': ['Intensity']*intensity_data.size + ['Unpleasantness']*unpleasantness_data.size})
        
        plt.figure(figsize=(7,6))
        sns.lineplot(data=data, x='Time', y='Rating', hue='Type', errorbar='sd')
        plt.show()

if group_plot:
    final_num_subs = 0
    for s, isub in enumerate(subjects):
        sub_id = 'sub-cr' + '{:02d}'.format(isub)
        print(sub_id)
        if isub in excluded:
            print('subject {} was excluded'.format(isub))
        else:
            final_num_subs += 1
    
            # load epochs
            sub_out_path = os.path.join(output_path, sub_id, 'rating')
            save_name = sub_out_path + os.sep + 'rating_' + task + '_epo.npy'
            epochs_all = np.load(save_name, allow_pickle=True)
            
            save_name = sub_out_path + os.sep + 'rating_' + task + '_all_events.csv'
            trial_infos = pd.read_csv(save_name, index_col=False, header=0)
            
            intensity_data = epochs_all[trial_infos.rating_type == 'intensity']
            unpleasantness_data = epochs_all[trial_infos.rating_type == 'unpleasantness']
        
            # average conditions of interest
            intensity_mean_all = np.nanmean(intensity_data, axis=0)
            unpleasantness_mean_all = np.nanmean(unpleasantness_data, axis=0)
            mean_all = np.nanmean(epochs_all, axis=0)
            if final_num_subs == 1:
                all_intensity_means = intensity_mean_all
                all_unpleasantness_means = unpleasantness_mean_all
                all_means = mean_all
            else:
                all_intensity_means = np.vstack((all_intensity_means, intensity_mean_all))
                all_unpleasantness_means = np.vstack((all_unpleasantness_means, unpleasantness_mean_all))
                all_means = np.vstack((all_means, mean_all))
    
    
    x = np.linspace(interval[0], interval[1], all_means.shape[1])
    
    stat_interval = [0, 30]
    start = int((stat_interval[0] - interval[0]) * desired_sr)
    stop = int((stat_interval[1] - interval[0]) * desired_sr)
    print('Testing intensity ratings')
    t_obs, clusters, cluster_pv, H0 = mne.stats.permutation_cluster_1samp_test(all_intensity_means[:, start:stop] -
                                                                               0,
                                                                               threshold=None, n_permutations=10000,
                                                                               seed=1990, out_type='mask', tail=1)
    for i_c, c in enumerate(clusters):
        c = c[0]
        if cluster_pv[i_c] < 0.05:
            print(f'Cluster {i_c+1} is significant (p={cluster_pv[i_c]})')
            print(f'It extends from {x[c.start + start]} to {x[c.stop + start]}')
        else:
            print(f'Cluster {i_c+1} is not significant (p={cluster_pv[i_c]})')
            print(f'It extends from {x[c.start + start]} to {x[c.stop + start]}')
    
    print('Testing unpleasantness ratings')
    t_obs, clusters, cluster_pv, H0 = mne.stats.permutation_cluster_1samp_test(all_unpleasantness_means[:, start:stop] -
                                                                               0,
                                                                               threshold=None, n_permutations=10000,
                                                                               seed=1992, out_type='mask', tail=1)
    for i_c, c in enumerate(clusters):
        c = c[0]
        if cluster_pv[i_c] < 0.05:
            print(f'Cluster {i_c+1} is significant (p={cluster_pv[i_c]})')
            print(f'It extends from {x[c.start + start]} to {x[c.stop + start]}')
        else:
            print(f'Cluster {i_c+1} is not significant (p={cluster_pv[i_c]})')
            print(f'It extends from {x[c.start + start]} to {x[c.stop + start]}')


    sns.set_context("poster")
    sns.set_style("ticks")

    fig2, ax2 = plt.subplots(figsize=(10, 8))
    x = np.linspace(interval[0], interval[1], all_means.shape[1])
    
    mean_intensity = np.mean(all_intensity_means, axis=0)
    sem_intensity = np.std(all_intensity_means, axis=0) / np.sqrt(all_intensity_means.shape[0])

    mean_unpleasantness = np.mean(all_unpleasantness_means, axis=0)
    sem_unpleasantness = np.std(all_unpleasantness_means, axis=0) / np.sqrt(all_unpleasantness_means.shape[0])
    
    rect = patches.Rectangle([0, 0], 30, 1, transform=ax2.get_xaxis_transform(),
                             edgecolor=None, facecolor='gray', alpha=0.2)
    ax2.add_patch(rect)
    #plt.axvline(x=0, color='black', linestyle='--')
    #plt.axvline(x=30, color='black', linestyle='--')

    # ax2.plot(x, mean_rep, color=my_color)
    # ax2.fill_between(x, mean_rep - sem_rep, mean_rep + sem_rep, alpha=0.2, color=my_color)
    ax2.plot(x, mean_intensity, color=color1, label='intensity')
    ax2.fill_between(x, mean_intensity - sem_intensity, 
                     mean_intensity + sem_intensity, 
                     alpha=0.2, color=color1)
    ax2.plot(x, mean_unpleasantness, color=color1, linestyle='dotted', label='unpleasantness')
    ax2.fill_between(x, mean_unpleasantness - sem_unpleasantness, 
                     mean_unpleasantness + sem_unpleasantness, 
                     alpha=0.2, color=color1)
    ax2.set(ylabel='Rating')
    ax2.set(xlabel='Time (seconds)')
    ax2.legend()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    #fig2.suptitle('Average continuous rating')
    plt.savefig(os.path.join(result_path, 'rating_mean.svg'), 
                transparent=True, format='svg', bbox_inches='tight')
    plt.show()
