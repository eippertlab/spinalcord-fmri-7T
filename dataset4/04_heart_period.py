import pandas as pd
import numpy as np
import os
import glob
import mne
from scipy.interpolate import RegularGridInterpolator
from scipy import signal
#from systole.utils import to_epochs
import matplotlib as mpl
# mpl.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from heart_GUI import r_peak_GUI

main_path = '/data/pt_02661_raw/Heatpain_behav/data'
raw_path = os.path.join(main_path, 'raw_data')
output_path = os.path.join(main_path, 'derivatives')
result_path = os.path.join(output_path, 'results')

new_sr = 10
baseline = (-5, 0)  # interval for baseline correction
interval = (-5, 35)  # interval for epoch
# how many seconds you want to display in the manual correction GUI
time_window = 20
task = 'continrating'

subjects = np.arange(1, 21)
excluded = [7, 12, 17]

color4 = [116/255, 196/255, 118/255]

preprocess = False
manual = False
hb2hp = False
epoch = False
plot = False
group_plot = True

new_rc_params = {"font.family": 'Arial', "font.size": 12, "font.serif": [],
                 "svg.fonttype": 'none'}
mpl.rcParams.update(new_rc_params)

all_subs_epochs = []
for isub in subjects:
    sub_id = 'sub-cr' + '{:02d}'.format(isub)
    print(sub_id)
    physio_files = glob.glob(os.path.join(raw_path, sub_id, 'beh', sub_id + '_task-' + task + '*physio.vhdr'))
    physio_files = np.array(sorted(physio_files))

    sub_out_path = os.path.join(output_path, sub_id, 'hpr')
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
            raw.pick_channels(ch_names=["ECG"])
            raw.set_channel_types(mapping={"ECG": 'ecg'})  # convert channel type from MISC to eeg
            raw_data = raw.get_data()
            raw_data = np.squeeze(raw_data)
            
            # extraction of onsets
            events, event_ids = mne.events_from_annotations(raw, regexp='^Stimulus')
            onsets = events[:, 0]
            
            # load event file
            event_file = os.path.join(raw_path, sub_id, 'beh', sub_id + '_task-' + task + '_run-' + str(i + 1) +
                                      '_events.csv')
            trial_infos = pd.read_csv(event_file, index_col=None)
            
            # find ECG events
            ecg_events, ch_ecg, avg_pulse = mne.preprocessing.find_ecg_events(raw)
            r_peaks = ecg_events[:, 0]
            
            # find maximum peaks as sometimes these are not at the real peak yet
            range = [-0.2 * sr, 0.2 * sr]
            for ii, peak in enumerate(r_peaks):
                bgn = int(peak + range[0])
                enn = int(peak + range[1])
                # you would think that you need to add bgn again to the index ind as
                # it could be relative to the interval start
                # but with the slicing the indices are called like in the original dataframe
                if enn <= len(raw_data) and bgn > 0:
                    temp = raw_data[bgn:enn]
                    ind = np.where(temp == np.max(temp))[0]
                    r_peaks[ii] = bgn + ind[0]
            # due to this step you sometimes assign the same maximum twice
            # delete duplicates
            r_peaks = np.unique(r_peaks)

            # sometimes it defines another peak (e.g. on T-waves)
            # find out if the median difference is somewhat reasonable
            ibi = np.diff(r_peaks)
            median_ibi = np.median(ibi) / sr
            # if there are too many peaks take the median peak amplitude and
            # remove all the peaks that are below that
            if median_ibi < 0.55:
                amps = raw_data[r_peaks]
                median_amp = np.median(amps)
                r_peaks = r_peaks[amps > median_amp]
            
            onsets_down = onsets / sr * new_sr
            
            if ((sub_id == 'sub-cr10') & (i == 1)):
                trial_infos = trial_infos.iloc[1:, :]
            trial_infos['onset_downsampled'] = onsets_down
            
            save_name = os.path.join(sub_out_path, 'hp_{}_raw_rpeaks_run_{}'.format(task, i+1))
            np.save(save_name, r_peaks)
            print('Peaks saved for run {}'.format(i+1))
            
            save_name = sub_out_path + os.sep + 'hp_{}_events_run_{}.csv'.format(task, i+1)
            trial_infos.to_csv(save_name, index=False)  
    
    if manual:
        for i, p in enumerate(physio_files):
            print('Displaying run {}'.format(i+1))
            # load physio file and select ECG channel
            raw = mne.io.read_raw(p, preload=True)
            channels = raw.info['ch_names']
            print(channels)
            sr = raw.info['sfreq']
            raw.pick_channels(ch_names=["ECG"])
            raw.set_channel_types(mapping={"ECG": 'ecg'})
            raw_data = raw.get_data()
            raw_data = np.squeeze(raw_data)

            # read previously created R peaks
            save_name = os.path.join(sub_out_path, 'hp_{}_raw_rpeaks_run_{}.npy'.format(task, i + 1))
            r_peaks = np.load(save_name)

            # put into GUI to manually correct
            interval_size = int(sr * time_window)
            save_name = os.path.join(sub_out_path, 'hp_{}_cleaned_rpeaks_run_{}.npy'.format(task, i + 1))
            my_GUI = r_peak_GUI(raw_data, r_peaks, interval_size, sr, save_name)
            
    if hb2hp:
        for i, p in enumerate(physio_files):
            # load physio file and select ECG channel
            raw = mne.io.read_raw(p, preload=True)
            sr = raw.info['sfreq']
            raw.pick_channels(ch_names=["ECG"])
            raw.set_channel_types(mapping={"ECG": 'ecg'})
            raw_data = raw.get_data()
            raw_data = np.squeeze(raw_data)

            # read saved cleaned R-peaks
            # read previously created R peaks
            save_name = os.path.join(sub_out_path, 'hp_{}_cleaned_rpeaks_run_{}.npy'.format(task, i + 1))
            r_peaks = np.load(save_name)

            # convert into heart period
            num_samples = raw_data.shape[0]
            hb = np.array([r/sr for r in r_peaks])
            ibi = np.diff(hb)
            idx = (ibi > 0.5) & (ibi < 1.6)  # limits for realistic values
            hp = ibi * 1000  # in ms
            new_time = np.arange(1.0 / new_sr, (num_samples - 1) / sr, 1 / new_sr)
            X = hb[np.where(idx)[0] + 1]
            Y = hp[idx]
            interfunc = RegularGridInterpolator([X], Y, method='linear', bounds_error=False, fill_value=None)
            heart_period = interfunc(new_time)

            # another bandpass filter (Castagnetti et al 2016 Psychophysiology; Paulus et al., 2016)
            low_cut = 0.01
            high_cut = 2
            low = low_cut / (sr / 2)
            high = high_cut / (sr / 2)
            b, a = signal.butter(2, [low_cut, high_cut], btype='bandpass', fs=new_sr)
            heart_period = signal.filtfilt(b, a, heart_period)

            # save the resulting heart period
            # read previously created R peaks
            save_name = os.path.join(sub_out_path, 'heart_period_{}_run_{}'.format(task, i + 1))
            np.save(save_name, heart_period)    
        
        
    if epoch:
        trial_infos_all = pd.DataFrame()
        for i, p in enumerate(physio_files):
            # load heart period
            save_name = os.path.join(sub_out_path, 'heart_period_{}_run_{}.npy'.format(task, i + 1))
            heart_period = np.load(save_name)
            
            # load event file
            save_name = sub_out_path + os.sep + 'hp_{}_events_run_{}.csv'.format(task, i+1)
            trial_infos = pd.read_csv(save_name, index_col=False)
            # for this subject BrainAmp was started too late and the baseline period is not in the data
            # --> exclude it
            if ((sub_id == 'sub-cr01') & (i == 2)):
                trial_infos = trial_infos.iloc[1:, :]
            
            # stack all meta data tables
            trial_infos_all = pd.concat([trial_infos_all, trial_infos])
            
            onsets = trial_infos['onset_downsampled']
            
            if len(onsets) > 0:
                # if a recording is a bit too short
                if any(onsets > len(heart_period)):
                    print('Data is too short to build epochs of all trigger events')
                    onsets = onsets[0:-1]
                    events = trial_infos[:-1]
                
                for j, onset in enumerate(onsets):
                    curr_epoch = heart_period[int(onset + interval[0]*new_sr):int(onset + interval[1]*new_sr)]
                    curr_epoch = curr_epoch - curr_epoch[-interval[0]*new_sr]
                    # stack all runs together
                    if (j == 0):
                        epochs_run = curr_epoch
                    else:
                        epochs_run = np.vstack((epochs_run, curr_epoch))
                # epochs_run, rejected = to_epochs(signal=heart_period, triggers_idx=np.array(onsets * new_sr, dtype=int),
                #                                  sfreq=new_sr, apply_baseline=(baseline[0], baseline[1]),
                #                                  tmin=interval[0], tmax=interval[1])

                # stack all runs together
                if (i == 0):
                    epochs_all = epochs_run
                else:
                    epochs_all = np.vstack((epochs_all, epochs_run))
                    
        save_name = os.path.join(sub_out_path, 'hp_' + task + '_all_epochs')
        np.save(save_name, epochs_all)
        save_name = os.path.join(sub_out_path, 'hp_' + task + '_all_events.csv')
        trial_infos_all.to_csv(save_name, index=False)

    if plot:
        print('plotting')
        save_name = sub_out_path + os.sep + 'hp_' + task + '_all_epochs.npy'
        epochs = np.load(save_name)
    
        # plotting average response over all blocks
        fig1, ax1 = plt.subplots(figsize=(7, 6))
        x = np.linspace(interval[0], interval[1], epochs_all.shape[1])
        mean_hp = np.nanmean(epochs_all, axis=0)
        sem_hp = np.nanstd(epochs_all, axis=0) / np.sqrt(epochs_all.shape[0])
        ax1.plot(x, mean_hp, color='red')
        ax1.fill_between(x, mean_hp - sem_hp, mean_hp + sem_hp, color='red', alpha=0.2)
        plt.ylabel(r'Heart period response (ms)', color='black', fontsize=20)
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
        plt.savefig(os.path.join(result_path, '{}_hpr_mean.png'.format(sub_id)))
        plt.show()
    
if group_plot:

    for s, isub in enumerate(subjects):
        sub_id = 'sub-cr' + '{:02d}'.format(isub)
        if isub in excluded:
            print('subject {} was excluded'.format(isub))
        else:
            print(sub_id)
    
            # load epochs
            save_name = os.path.join(output_path, sub_id, 'hpr', 'hp_' + task + '_all_epochs.npy')
            epochs_all = np.load(save_name)
    
            # z-transform the time series to account for subjects variance in amplitude
            # epochs_all = stats.zscore(epochs_all, axis=None, nan_policy='omit')
    
            # average
            mean_all = np.nanmean(epochs_all, axis=0)
            if s == 0:
                all_means = mean_all
            else:
                all_means = np.vstack((all_means, mean_all))

    sns.set_context("poster")
    sns.set_style("ticks")

    fig1, ax1 = plt.subplots(figsize=(10, 8))
    x = np.linspace(interval[0], interval[1], all_means.shape[1])
    mean_hp = np.nanmean(all_means, axis=0)
    sem_hp = np.nanstd(all_means, axis=0) / np.sqrt(all_means.shape[0])

    rect = patches.Rectangle([0, 0], 30, 1, transform=ax1.get_xaxis_transform(),
                             edgecolor=None, facecolor='gray', alpha=0.2)
    ax1.add_patch(rect)
    #plt.axvline(x=0, color='black', linestyle='--')
    #plt.axvline(x=30, color='black', linestyle='--')
    
    ax1.plot(x, mean_hp, color=color4)
    ax1.fill_between(x, mean_hp - sem_hp, mean_hp + sem_hp, 
                     alpha=0.4, color=color4)
    ax1.set(ylabel='Heart period (ms)')
    ax1.set(xlabel='Time (seconds)')

    # ax1.axvline(x=0, linestyle="--", c=my_color)
    # ax1.axvline(x=30, linestyle="--", c=my_color)
    
    plt.ylim([-30, 30])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    #fig1.suptitle('Average heart period response')
    plt.savefig(os.path.join(result_path, 'hpr_mean.svg'), 
                transparent=True, format='svg', bbox_inches='tight')
    plt.show()

