#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analysis of pupil dilation responses
----------------------------------------
low-pass filter pupil data with 5th order butterworth cut-off 0.01 Hz,
opens GUI for manual blink correction,
epoch data according to events.tsv
and plot resulting responses for each subject and over all subjects

authors:
--------
Ulrike Horn

contact:
--------
uhorn@cbs.mpg.de

date:
-----
15th Jan 2025
"""

import os
import glob
import mne
import mne.stats
import numpy as np
import pandas as pd
import json
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pupil_GUI import interpolate_blinks_GUI
import seaborn as sns

raw_path = '/data/pt_02661_raw/Heatpain_behav/data/raw_data'
output_path = '/data/pt_02661_raw/Heatpain_behav/data/derivatives'
result_path = '/data/pt_02661_raw/Heatpain_behav/data/derivatives/results'

interval = [-5, 35]  # interval for cutting and plotting the data
baseline_interval = [-0.1, 0]
task = 'continrating'

color3 = [65/255, 171/255, 93/255]

subjects = np.arange(1, 21)
excluded = []

preprocess = False
manual = False
epoch = False
delete_bad_trials = False
subject_plot = False
count_bad = False
group_plot = True

np.random.seed(1990)

new_rc_params = {"font.family": 'Arial', "font.size": 12, "font.serif": [],
                 "svg.fonttype": 'none'}
mpl.rcParams.update(new_rc_params)

# how many seconds you want to display in the manual correction GUI
time_window = 20
interpolation_buffer = [0.05, 0.2]

if not os.path.exists(result_path):
    os.makedirs(result_path)

for isub in subjects:
    sub_id = 'sub-cr' + '{:02d}'.format(isub)
    print(sub_id)
    
    sub_out_path = os.path.join(output_path, sub_id, 'pupil')
    if not os.path.exists(sub_out_path):
        os.makedirs(sub_out_path)
        

    physio_files = glob.glob(os.path.join(raw_path, sub_id, 'beh', '*' + task +
                                          '*eyetracking_physio.asc'))
    physio_files = np.array(sorted(physio_files))

    if preprocess:
        for i, p in enumerate(physio_files):
            raw_et = mne.io.read_raw_eyelink(p, create_annotations=["blinks", "messages"])
            SR = raw_et.info['sfreq']
            channels = raw_et.info['ch_names']
            events_from_annot, event_dict = mne.events_from_annotations(raw_et,
                                                                        event_id={'Heat onset': 1})
            # collect the raw data (we need this one for manual correction display)
            data_raw = raw_et.get_data()
            data_pupil_raw = data_raw[channels.index('pupil_right'), :].transpose()

            # interpolate data within certain window before and after (usually 0.05s before and 0.2s after)
            mne.preprocessing.eyetracking.interpolate_blinks(raw_et,
                                                             buffer=(interpolation_buffer[0],
                                                                     interpolation_buffer[1]),
                                                             interpolate_gaze=True)

            # we want to store the information where the data has been interpolated
            # get blink events and durations
            blink_events_from_annot, blink_event_dict = mne.events_from_annotations(raw_et,
                                                                                    event_id={'blink': 5})
            blink_durs = []
            for a in raw_et.annotations:
                if a['description'] == 'blink':
                    blink_durs.append(a['duration'])

            # create a trigger column
            # compare with saved events file
            if (isub == 1) & (i > 0):
                event_file = os.path.join(raw_path, sub_id, 'beh', 
                                          sub_id + '_task-' + 
                                          task + '_run-' + str(i + 2) +
                                          '_events.csv')
            else:
                event_file = os.path.join(raw_path, sub_id, 'beh', 
                                          sub_id + '_task-' + 
                                          task + '_run-' + str(i + 1) +
                                          '_events.csv')
            trial_infos = pd.read_csv(event_file)
            
            assert len(trial_infos) == len(events_from_annot), "Event trigger and trial info are not the same size"
            
            trigger = np.zeros((len(data_pupil_raw), 1))
            trigger[events_from_annot[:, 0]] = 1

            # create a blink column
            blinks = np.zeros((len(data_pupil_raw), 1))
            for b, blink_dur in enumerate(blink_durs):
                onset = blink_events_from_annot[b, 0]
                blinks[onset:onset + int(blink_dur * SR)] = 1

            # the blinks were interpolated in a region around the actual blink
            # extend the sections accordingly in this column
            def extend_ones_sections(arr, samples_before, samples_after):
                # Ensure the array is a numpy array
                arr = np.asarray(arr)
                # find start and end of sections
                starts = []
                ends = []
                in_section = False
                for i, val in enumerate(arr):
                    if val == 1 and not in_section:
                        starts.append(i)
                        in_section = True
                    elif val == 0 and in_section:
                        ends.append(i)
                        in_section = False
                if in_section:  # Handle case where the last section goes till the end
                    ends.append(len(arr))

                # extend all starts and ends
                extended = np.zeros_like(arr)
                for start, end in zip(starts, ends):
                    extended[max(0, start - samples_before):min(len(arr), end + samples_after)] = 1
                return extended


            extended_blinks = extend_ones_sections(blinks, int(interpolation_buffer[0] * SR),
                                                   int(interpolation_buffer[1] * SR))

            # filter with butterworth low pass
            cutoff_freq = 2
            iir_params = dict(order=1, ftype='butter', output='sos')
            filtered = raw_et.filter(h_freq=cutoff_freq, l_freq=None,
                                     method="iir", iir_params=iir_params)
            data = filtered.get_data()
            data_pupil = data[channels.index('pupil_right'), :].transpose()

            # save data and blink array
            # filtered.save(os.path.join(sub_out_path, 'filtered_interpolated_data.fif'))
            df = pd.DataFrame({
                'raw_data': data_pupil_raw,
                'interpolated_data': data_pupil,
                'interpolated_bool': np.squeeze(extended_blinks),
                'event_trigger': np.squeeze(trigger)})
            save_name = os.path.join(sub_out_path, 'pupil_' + task + '_filtered_pupil_' +
                                     'run_' + str(i + 1) + '.csv')
            df.to_csv(save_name)
            # save the sampling frequency in a json
            my_json = {'SamplingFrequency': SR}
            json_data = json.dumps(my_json, indent=4)
            with open(os.path.join(sub_out_path, 'pupil_' + task + '_info_run_' +
                                                 str(i + 1) + '.json'), 'w') as outfile:
                outfile.write(json_data)

    if manual:
        for i, p in enumerate(physio_files):
            print(p)
            # read accompanying json file
            with open(os.path.join(sub_out_path, 'pupil_' + task + '_info_run_' +
                                                 str(i + 1) + '.json'), 'r') as f:
                json_dict = json.load(f)
            sr = json_dict['SamplingFrequency']
            interval_size = int(sr * time_window)
            save_name = os.path.join(sub_out_path, 'pupil_' + task + '_filtered_pupil_' +
                                     'run_' + str(i + 1) + '.csv')
            df = pd.read_csv(save_name)
            save_name = os.path.join(sub_out_path, 'pupil_' + task + '_cleaned_pupil_' +
                                     'run_' + str(i + 1) + '.csv')
            raw_data = df['raw_data'].to_numpy()
            filt_data = df['interpolated_data'].to_numpy()
            interpolated_bool = df['interpolated_bool'].to_numpy()
            # open GUI and interpolate data at missed blinks
            my_GUI = interpolate_blinks_GUI(filt_data, raw_data, interpolated_bool,
                                            interval_size, sr, save_name)
            # open created file again and add previously stored info again
            df_new = pd.read_csv(save_name, index_col=0)
            df_new['event_trigger'] = df['event_trigger']
            df_new.to_csv(save_name, index=False)

    if epoch:
        trial_infos_all = pd.DataFrame()
        all_epochs = []

        for i, p in enumerate(physio_files):
            # search for filtered + manually corrected data
            try:
                save_name = os.path.join(sub_out_path, 'pupil_' + task +
                                         '_cleaned_pupil_run_' + str(i + 1) + '.csv')
                df = pd.read_csv(save_name)
                filt_data = df['interpolated_data'].to_numpy()
            except FileNotFoundError:
                print("Error: Manually corrected data does not appear to exist")
                raise

            interpolated_bool = df['interpolated_bool'].to_numpy()
            # read accompanying json file
            with open(os.path.join(sub_out_path, 'pupil_' + task + '_info_run_' +
                                                 str(i + 1) + '.json'), 'r') as f:
                json_dict = json.load(f)
            sr = json_dict['SamplingFrequency']

            # read events.tsv
            if (isub == 1) & (i > 0):
                event_file = os.path.join(raw_path, sub_id, 'beh', 
                                          sub_id + '_task-' + 
                                          task + '_run-' + str(i + 2) +
                                          '_events.csv')
            else:
                event_file = os.path.join(raw_path, sub_id, 'beh', 
                                          sub_id + '_task-' + 
                                          task + '_run-' + str(i + 1) +
                                          '_events.csv')
            trial_infos = pd.read_csv(event_file)
            if not sum(df['event_trigger']) == len(trial_infos):
                raise Exception("Error: File lengths don't match")

            # stack all meta data tables
            trial_infos_all = pd.concat([trial_infos_all, trial_infos])

            onsets = df[df['event_trigger'] == 1].index

            # epoching
            if len(onsets) > 0:
                for e, event in enumerate(onsets):
                    pd_epo = filt_data[int((event + interval[0] * sr)):
                                       int((event+ interval[1] * sr))]
                    pd_epo = pd_epo - filt_data[int(event)]
                    blink_epo = interpolated_bool[int((event + interval[0] * sr)):
                                                  int((event + interval[1] * sr))]
                    if e == 0:
                        epochs_run = pd_epo
                        blink_epochs_run = blink_epo
                    else:
                        epochs_run = np.vstack((epochs_run, pd_epo))
                        blink_epochs_run = np.vstack((blink_epochs_run, blink_epo))
                # stack all runs together
                if i == 0:
                    epochs_all = epochs_run
                    blink_epochs_all = blink_epochs_run
                else:
                    epochs_all = np.vstack((epochs_all, epochs_run))
                    blink_epochs_all = np.vstack((blink_epochs_all, blink_epochs_run))
            else:
                print('one run had no valid trials')

        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_all_epochs')
        np.save(save_name, epochs_all)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_blinks_all_epochs')
        np.save(save_name, blink_epochs_all)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_all_events.csv')
        trial_infos_all.to_csv(save_name, index=False)

        # also save z-scored epochs for each subject
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_all_epochs_zscored')
        epochs_all_z = stats.zscore(epochs_all, axis=None, nan_policy='omit')
        np.save(save_name, epochs_all_z)

    if delete_bad_trials:
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_all_epochs.npy')
        epochs_all = np.load(save_name)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_all_epochs_zscored.npy')
        epochs_all_z = np.load(save_name)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_all_events.csv')
        events_all = pd.read_csv(save_name)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_blinks_all_epochs.npy')
        blinks_all = np.load(save_name)
        num_trials = epochs_all.shape[0]
        num_samples = epochs_all.shape[1]
        bad = np.zeros(num_trials)
        for trial in range(num_trials):
            if sum(blinks_all[trial]) > 0.5 * num_samples:
                bad[trial] = 1
        print('{} trials had to be excluded because of too much interpolation'.format(sum(bad)))
        blocks = np.unique(events_all['run'])
        for block in blocks:
            if sum(bad[events_all['run'] == block]) > 0.5 * len(bad[events_all['run'] == block]):
                bad[events_all['run'] == block] = 1
        print('{} trials had to be excluded when we delete whole blocks.'.format(sum(bad)))
        bad_df = pd.DataFrame(bad)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_bad_events.csv')
        bad_df.to_csv(save_name, index=False, header=False)
        events_all = events_all.drop(np.where(bad == 1)[0]).reset_index(drop=True)
        epochs_all = np.delete(epochs_all, np.where(bad == 1)[0], axis=0)
        epochs_all_z = np.delete(epochs_all_z, np.where(bad == 1)[0], axis=0)
        blinks_all = np.delete(blinks_all, np.where(bad == 1)[0], axis=0)
        print('We have {} good trials remaining.'.format(epochs_all.shape[0]))

        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_epochs')
        np.save(save_name, epochs_all)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_epochs_zscored')
        np.save(save_name, epochs_all_z)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_blinks_good_epochs')
        np.save(save_name, blinks_all)
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_events.csv')
        events_all.to_csv(save_name, index=False)

    if subject_plot:
        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_epochs.npy')
        epochs_all = np.load(save_name)

        save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_events.csv')
        events = pd.read_csv(save_name)
        
        # plotting average response over all blocks
        sns.set_context('talk')
        fig1, ax1 = plt.subplots(figsize=(7, 6))
        x = np.linspace(interval[0], interval[1], epochs_all.shape[1])
        mean_hp = np.nanmean(epochs_all, axis=0)
        sem_hp = np.nanstd(epochs_all, axis=0) / np.sqrt(epochs_all.shape[0])
        ax1.plot(x, mean_hp, color='red')
        ax1.fill_between(x, mean_hp - sem_hp, mean_hp + sem_hp, color='red', alpha=0.2)
        plt.ylabel(r'Pupil dilation (a.u.)', color='black', fontsize=20)
        plt.xlabel('Time (s)', color='black', fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        # To specify the number of ticks on both or any single axes
        plt.locator_params(axis='y', nbins=4)
        plt.locator_params(axis='x', nbins=6)
        ax1.axvline(x=0, linestyle="--", c="black")
        ax1.tick_params(color='black', labelcolor='black')
        plt.tight_layout()
        for spine in ax1.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(1.5)
        plt.savefig(os.path.join(result_path, '{}_pupil_dilation_mean.png'.format(sub_id)))
        plt.show()
        
        
if count_bad:
    n_trial_all = 0
    n_trial_bad = 0
    for isub in subjects:
        sub_id = 'sub-cr' + '{:02d}'.format(isub)
        print('----------' + sub_id + '--------')
        sub_out_path = os.path.join(output_path, sub_id, 'pupil')
        bad_df = pd.read_csv(sub_out_path + os.sep + 'pupil_' + task + '_bad_events.csv',
                             header=None, names=['bad'])
        n_trial_all += len(bad_df)
        n_trial_bad += sum(bad_df['bad'])
        print('{} out of {} trials are bad'.format(sum(bad_df['bad']), len(bad_df)))
    print('Overall from all subjects we have {} epochs'.format(n_trial_all))
    print('Out of these {} epochs are bad and had to be discarded'.format(n_trial_bad))
    print('This equals {} %'.format(round(n_trial_bad/n_trial_all*100, 3)))

if group_plot:
    sns.set_style('white')
    sns.set_context('talk')
    final_num_subs = 0
    result_df = pd.DataFrame([])
    # collect subjects that are no outliers
    for isub in subjects:
        if isub in excluded:
            print('subject {} was excluded'.format(isub))
        else:
            sub_id = 'sub-cr' + '{:02d}'.format(isub)
            print(sub_id)
            final_num_subs += 1
            sub_out_path = os.path.join(output_path, sub_id, 'pupil')
            save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_epochs_zscored.npy')
            epochs = np.load(save_name)
            save_name = os.path.join(sub_out_path, 'pupil_' + task + '_good_events.csv')
            df = pd.read_csv(save_name)
            mean_all = np.nanmean(epochs, axis=0)
            if final_num_subs == 1:
                all_means = mean_all
            else:
                all_means = np.vstack((all_means, mean_all))
    
    x = np.linspace(interval[0], interval[1], all_means.shape[1])
    
    sr = all_means.shape[1]/(interval[1]-interval[0])
    stat_interval = [0, 30]
    start = int((stat_interval[0] - interval[0]) * sr)
    stop = int((stat_interval[1] - interval[0]) * sr)
    t_obs, clusters, cluster_pv, H0 = mne.stats.permutation_cluster_1samp_test(all_means[:, start:stop] -
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
            
    sns.set_context("poster")
    sns.set_style("ticks")
    
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    mean_rep = np.mean(all_means, axis=0)
    sem_rep = np.std(all_means, axis=0) / np.sqrt(all_means.shape[0])

    rect = patches.Rectangle([0, 0], 30, 1, transform=ax2.get_xaxis_transform(),
                             edgecolor=None, facecolor='gray', alpha=0.2)
    ax2.add_patch(rect)
    #plt.axvline(x=0, color='black', linestyle='--')
    #plt.axvline(x=30, color='black', linestyle='--')
    
    ax2.plot(x, mean_rep, color=color3)
    ax2.fill_between(x, mean_rep - sem_rep, mean_rep + sem_rep, 
                     alpha=0.4, color=color3)
    ax2.set(ylabel='Pupil dilation (z-scored)')
    ax2.set(xlabel='Time (seconds)')

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # fig2.suptitle('Average pupil dilation')
    plt.savefig(os.path.join(result_path, 'pdr_mean_zscore.svg'), 
                transparent=True, format='svg', bbox_inches='tight')
    plt.show()

