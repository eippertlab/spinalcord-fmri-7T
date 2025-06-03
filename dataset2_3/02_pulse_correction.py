# Displays the pulse data to manually correct R-peaks
# needs helper function class_GUI_hb.py
# code by Ulrike Horn

import glob
import os
import sys
import json
import numpy as np
import pandas as pd
import scipy.signal
from helper_functions.class_GUI_hb import r_peak_GUI

# subjects = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

task = 'longpain'

interval_size = 20  # plotted interval in s

data_path = '/data/pt_02661_raw/Heatpain/raw_data/'
output_path = '/data/pt_02661_raw/Heatpain/derivatives/'

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
    # Dataset 1 and 2 are organized in different folders
    if isub < 17:
        func_path = data_path + 'dataset2' + os.sep + sub + os.sep + 'func'
    else:
        func_path = data_path + 'dataset3' + os.sep + sub + os.sep + 'func'
    # load pulse data for each run
    for i in range(num_runs):
        run = i + 1
        print('Starting with run: ' + str(run))
        if isub < 17:
            pulse_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) + '*_recording-pulsoxy_physio.tsv.gz')
        else:
            pulse_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) + '*_recording-pulsoxyrespiration_physio.tsv.gz')

        if len(pulse_file) == 1:
            pulse_file = pulse_file[0]
        else:
            raise Exception('Too many or too few pulseoxy files for this run: ' + str(pulse_file))

        save_path = output_path + sub + os.sep + 'func' + os.sep + 'run-' + str(run)
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        tmp = pulse_file.split('.')
        name = tmp[0]

        with open(name + '.json', 'r') as f:
            json_dict = json.load(f)
        sr = json_dict['SamplingFrequency']
        column_names = json_dict['Columns']
        start_time = json_dict['StartTime']

        df = pd.read_csv(pulse_file, header=None, sep='\t', compression='gzip', names=column_names)
        raw_data = df['cardiac']

        min_height = min(raw_data) + (max(raw_data) - min(raw_data))/3  # at least 1/3 of distance between highest and lowest data point
        peaks = scipy.signal.find_peaks(raw_data, distance=sr/2, height=min_height)[0]  # detects only peaks which
        # are min. 0.5s apart

        bad_starts = np.array([])
        bad_ends = np.array([])

        save_name = save_path + os.sep + 'cleaned_ecg_events_run_' + str(run) + '.csv'
        my_GUI = r_peak_GUI(raw_data, peaks, bad_starts, bad_ends, interval_size*sr, sr, save_name)

        # put the trigger back in for some checking
        new_df = pd.read_csv(save_name)
        new_df['trigger'] = df['trigger']
        new_df.to_csv(save_name, index=False)

