# Analysis of skin conductance

import numpy as np
import pandas as pd
import os
import glob
import mne
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

color2 = [35/255, 139/255, 69/255]

if not os.path.exists(output_path):
    os.makedirs(output_path)

if not os.path.exists(result_path):
    os.makedirs(result_path)

interval = [-5, 35]  # interval for cutting and plotting the data
baseline_interval = [-0.1, 0]
new_sr = 100  # Hz that data will be downsampled to
task = 'continrating'

subjects = np.arange(1, 21)

excluded = []

preprocess = False
outlier = False
plot = False
group_plot = True

new_rc_params = {"font.family": 'Arial', "font.size": 12, "font.serif": [],
                 "svg.fonttype": 'none'}
mpl.rcParams.update(new_rc_params)

for isub in subjects:
    sub_id = 'sub-cr' + '{:02d}'.format(isub)
    print(sub_id)
    physio_files = glob.glob(os.path.join(raw_path, sub_id, 'beh', sub_id + '_task-' + task + '*physio.vhdr'))
    physio_files = np.array(sorted(physio_files))

    sub_out_path = os.path.join(output_path, sub_id, 'scr')
    if not os.path.exists(sub_out_path):
        os.makedirs(sub_out_path)

    if preprocess:
        trial_infos_all = pd.DataFrame()
        all_epochs = []
        for i, p in enumerate(physio_files):
            raw = mne.io.read_raw(p, preload=True)
            channels = raw.info['ch_names']
            print(channels)
            sr = raw.info['sfreq']
            raw.pick_channels(ch_names=["GSR_MR_100_left_hand"])
            raw.set_channel_types(mapping={"GSR_MR_100_left_hand": 'eeg'})  # convert channel type from MISC to eeg

            # load event file
            event_file = os.path.join(raw_path, sub_id, 'beh', sub_id + '_task-' + task + '_run-' + str(i + 1) +
                                      '_events.csv')
            trial_infos = pd.read_csv(event_file, index_col=None)

            # downsampling to 100 Hz and extraction of onsets
            gsr_downsampled = raw.resample(new_sr, npad="auto")
            events, event_ids = mne.events_from_annotations(gsr_downsampled, regexp='^Stimulus')
            onsets = events[:, 0]

            # 1st order butterworth filter (iir) band-pass with 0.0159 and 5 Hz
            l_freq = 0.0159  # lower cutoff frequency
            #h_freq = 5  # upper cutoff frequency
            h_freq = 1  # upper cutoff frequency
            iir_params = dict(order=1, ftype='butter', output='sos')
            filtered = gsr_downsampled.filter(l_freq=l_freq, h_freq=h_freq, picks=["GSR_MR_100_left_hand"],
                                              method="iir", iir_params=iir_params)

            # cut into epochs from -1 to 8 using the given baseline for correction
            epochs = mne.Epochs(filtered, events, event_id=event_ids,
                                baseline=(baseline_interval[0], baseline_interval[1]),
                                tmin=interval[0], tmax=interval[1], preload=True)
            all_epochs.append(epochs)

            # for this subject BrainAmp was started too late and the baseline period is not in the data
            # --> MNE excludes it
            if ((sub_id == 'sub-cr01') & (i == 2)) | ((sub_id == 'sub-cr10') & (i == 1)):
                trial_infos = trial_infos.iloc[1:, :]

            # stack all meta data tables
            trial_infos_all = pd.concat([trial_infos_all, trial_infos])

        epochs_all = mne.concatenate_epochs(all_epochs)

        assert len(epochs_all) == len(trial_infos_all), "Epochs and trial info are not the same size"

        save_name = sub_out_path + os.sep + 'scr_' + task + '_epo.fif'
        epochs_all.save(save_name, overwrite=True)
        save_name = sub_out_path + os.sep + 'scr_' + task + '_all_events.csv'
        trial_infos_all.to_csv(save_name, index=False)

    if outlier:
        save_name = sub_out_path + os.sep + 'scr_' + task + '_epo.fif'
        epochs = mne.read_epochs(save_name, preload=True)
        epochs_all = epochs.get_data()  # n_epochs x n_channels x n_times
        epochs_all = np.squeeze(epochs_all)  # epochs x time points
        epochs_all = epochs_all * 1e6  # convert unit: V -> µS
        save_name = sub_out_path + os.sep + 'scr_' + task + '_all_events.csv'
        events_all = pd.read_csv(save_name)
        discarded = 0
        for i, block in enumerate(events_all.run.unique()):
            epochs_run = epochs_all[events_all['run'] == block]
            events_run = events_all.loc[events_all['run'] == block, :].copy()
            # within interval of 1 to 8 s after onset there should be something larger than baseline in more than 25%
            # of the trials of a block
            responses = np.zeros([len(events_run)])
            for j in range(len(events_run)):
                any_response_interval = [1, 30]  # relative to onset 0 which is dependent on epoch cutting interval
                baseline_interval = [-1, 0]
                # find maximum in response interval
                max_ind = np.argmax(epochs_run[j, -interval[0] * new_sr + any_response_interval[0] * new_sr:
                                                  -interval[0] * new_sr + any_response_interval[1] * new_sr])
                # maximum index in absolute interval terms
                abs_max_ind = max_ind - interval[0] * new_sr + any_response_interval[0] * new_sr
                amp = epochs_run[j, abs_max_ind]
                # baseline value
                base = np.mean(epochs_run[j, -interval[0] * new_sr + baseline_interval[0] * new_sr:
                                             -interval[0] * new_sr + baseline_interval[1] * new_sr])
                if amp - base > 0.01:
                    responses[j] = 1
            if sum(responses) < 0.25 * len(responses):
                events_run['discard_block'] = 1
                print('block discarded')
                discarded += 1
            else:
                events_run['discard_block'] = 0
                print('good')
            if i == 0:
                events_all_new = events_run
            else:
                events_all_new = pd.concat([events_all_new, events_run])
        if discarded/len(events_all.run.unique()) > 0.5:
            print('subject needs to be excluded')
            events_all_new['discard_subject'] = 1
        else:
            events_all_new['discard_subject'] = 0
        
        save_name = sub_out_path + os.sep + 'scr_' + task + '_all_events_outlier.csv'
        events_all_new.to_csv(save_name, index=False)

    if plot:
        save_name = sub_out_path + os.sep + 'scr_' + task + '_epo.fif'
        epochs = mne.read_epochs(save_name, preload=True)
        epochs_all = epochs.get_data()  # n_epochs x n_channels x n_times
        epochs_all = np.squeeze(epochs_all)  # epochs x time points
        epochs_all = epochs_all * 1e6  # convert unit: V -> µS
        # save_name = os.path.join(output_path, sub_id, 'scr', 'scr_all_epochs.npy')
        # epochs_all = np.load(save_name)

        # plotting average response over all blocks
        fig1, ax1 = plt.subplots(figsize=(7, 6))
        x = np.linspace(interval[0], interval[1], epochs_all.shape[1])
        mean_hp = np.nanmean(epochs_all, axis=0)
        sem_hp = np.nanstd(epochs_all, axis=0) / np.sqrt(epochs_all.shape[0])
        ax1.plot(x, mean_hp, color='red')
        ax1.fill_between(x, mean_hp - sem_hp, mean_hp + sem_hp, color='red', alpha=0.2)
        plt.ylabel(r'Skin conductance amplitude ($\mu$S)', color='black', fontsize=20)
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
        plt.savefig(os.path.join(result_path, '{}_scr_mean_onset_only.png'.format(sub_id)))
        plt.show()

if group_plot:
    final_num_subs = 0
    for s, isub in enumerate(subjects):
        sub_id = 'sub-cr' + '{:02d}'.format(isub)
        print(sub_id)
        save_name = sub_out_path + os.sep + 'scr_' + task + '_all_events_outlier.csv'
        events_all = pd.read_csv(save_name)

        if events_all['discard_subject'].unique()[0]:
            print('Subject {} excluded during SCR analysis'.format(isub))
            continue
        
        final_num_subs += 1

        # load epochs
        sub_out_path = os.path.join(output_path, sub_id, 'scr')
        save_name = sub_out_path + os.sep + 'scr_' + task + '_epo.fif'
        epochs = mne.read_epochs(save_name, preload=True)
        epochs_all = epochs.get_data()  # n_epochs x n_channels x n_times
        epochs_all = np.squeeze(epochs_all)  # epochs x time points
        epochs_all = epochs_all * 1e6  # convert unit: V -> µS

        # z-transform the time series to account for subjects variance in amplitude
        #epochs_all = stats.zscore(epochs_all, axis=None, nan_policy='omit')

        # average conditions of interest
        mean_all = np.nanmean(epochs_all, axis=0)
        if final_num_subs == 1:
            all_means = mean_all
        else:
            all_means = np.vstack((all_means, mean_all))
    
    x = np.linspace(interval[0], interval[1], all_means.shape[1])
    
    stat_interval = [1, 30]
    start = int((stat_interval[0] - interval[0]) * new_sr)
    stop = int((stat_interval[1] - interval[0]) * new_sr)
    t_obs, clusters, cluster_pv, H0 = mne.stats.permutation_cluster_1samp_test(all_means[:, start:stop] -
                                                                               0,
                                                                               threshold=None, n_permutations=10000,
                                                                               seed=1991, out_type='mask', tail=1)
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

    ax2.plot(x, mean_rep, color=color2)
    ax2.fill_between(x, mean_rep - sem_rep, mean_rep + sem_rep, 
                     alpha=0.4, color=color2)
    ax2.set(ylabel=r'Skin conductance ($\mu$S)')
    ax2.set(xlabel='Time (seconds)')
    plt.ylim([-0.1, 0.6])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # fig2.suptitle('Average skin conductance response')
    plt.savefig(os.path.join(result_path, 'scr_mean.svg'), 
                transparent=True, format='svg', bbox_inches='tight')
    plt.show()
