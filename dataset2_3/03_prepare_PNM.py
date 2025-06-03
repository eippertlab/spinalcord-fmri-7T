# loads the saved pulse peak data
# if there are bad intervals interpolate or fill them with the mean (depending on length)
# load respiration data and plot it to determine data quality
# if there are several options the user can choose the best respiration data
# as output a file physio.txt is created that contains the columns:
# trigger
# cardiac
# respiration
# to be used in PNM.

# code by Ulrike Horn and Lisa-Marie Pohle

import glob
import json
import os
import statistics
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal, stats

# subjects = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

task = 'longpain'

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

    for irun in range(num_runs):
        run = irun + 1
        print('Run: ' + str(run))
        save_path = output_path + sub + os.sep + 'func' + os.sep + 'run-' + str(run)

        # load info about pulse data
        if isub < 17:
            json_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) + '*_recording-pulsoxy_physio.json')[0]
        else:
            json_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) + '*_recording-pulsoxyrespiration_physio.json')[0]
        with open(json_file, 'r') as f:
            json_dict = json.load(f)
        sr_pulse = json_dict['SamplingFrequency']
        start_time_pulse = json_dict['StartTime']

        # load corrected pulse data
        saved_file = save_path + os.sep + 'cleaned_ecg_events_run_' + str(run) + '.csv'
        pulse_data = pd.read_csv(saved_file)
        # correct bad intervals if there are any
        starts = np.where(pulse_data.bad_intervals_1start_2end == 1)[0]
        if starts.size == 0:
            bad_interval_pulse = np.zeros(len(pulse_data))
            print('Subject: ' + sub + ', run: ' + str(run) + ' --- No bad intervals found!')
        else:
            bad_interval_pulse = np.zeros(len(pulse_data))
            starts_plot = np.subtract(starts, 2*sr_pulse)  # start plotting 2s earlier
            # except when it is at the end of the data
            starts_plot[starts_plot < 0] = 0

            ends = np.where(pulse_data.bad_intervals_1start_2end == 2)[0]
            ends_plot = np.add(ends, 2*sr_pulse)  # end plotting 2s later
            # except when it is at the end of the data
            ends_plot[ends_plot > len(pulse_data)] = len(pulse_data) - 1
            for interval in range(len(starts)):
                time = np.arange(starts_plot[interval], ends_plot[interval]+1)
                time = np.divide(time, sr_pulse)
                raw_data = pulse_data.loc[starts_plot[interval]:ends_plot[interval], 'cardiac']
                plt.plot(time, raw_data)
                plt.show(block=True)

                check = False
                while not check:
                    todo = input('What do you want to do with the data? [interpolate, ignore]')
                    if todo == 'interpolate':
                        # if there are any peaks in bad interval, delete them first
                        pulse_data.loc[starts[interval]:ends[interval], 'peaks'] = 0

                        # interpolate peaks based on median IBI
                        beats = np.where(pulse_data.peaks == 1)[0]
                        ibis = np.diff(beats)
                        median_ibi = statistics.median(ibis)
                        print('Median IBI is: ' + str(median_ibi/sr_pulse))
                        # if the last peak before is defined, then interpolate normally
                        # (just go until the end instead of until the next peak)
                        # if only the peak after is defined (because beginning is bad) then
                        # change interpolation to go in the other direction
                        try:
                            last_peak_before = pulse_data.loc[:starts[interval]].index[pulse_data.loc[:starts[interval]].peaks == 1][-1]
                            last_peak_before_missing = False
                        except IndexError:
                            print("There is no peak before the interval")
                            last_peak_before_missing = True
                            last_peak_before = starts[interval]
                        try:
                            first_peak_after_missing = False
                            first_peak_after = pulse_data.loc[ends[interval]:].index[pulse_data.loc[ends[interval]:].peaks == 1][0]
                        except IndexError:
                            print("There is no peak after the interval")
                            first_peak_after_missing = True
                            first_peak_after = ends[interval]

                        # how much time is missing and how much peaks can fit in there
                        time_to_interpolate = (first_peak_after-last_peak_before)/sr_pulse
                        print('Time to interpolate: ' + str(time_to_interpolate))
                        peak_count = round((first_peak_after - last_peak_before) / median_ibi - 1)
                        # -1 to get to number of peaks instead of number of ibis
                        print('No of peaks to interpolate: ' + str(peak_count))
                        # start from first peak that is there and go backwards
                        # with median IBI
                        if last_peak_before_missing:
                            used_ibi = median_ibi
                            for i in range(peak_count):
                                pulse_data.at[int(first_peak_after - (i+1)*used_ibi), 'peaks'] = 1
                        # start from last peak that is there and go forwards
                        # with median IBI
                        elif first_peak_after_missing:
                            used_ibi = median_ibi
                            for i in range(peak_count):
                                pulse_data.at[int(last_peak_before + (i+1)*used_ibi), 'peaks'] = 1
                        # if both peaks are there choose an interbeat-interval so that peaks are
                        # equally distributed between them
                        else:
                            used_ibi = round((first_peak_after-last_peak_before)/(peak_count+1))
                            for i in range(peak_count):
                                pulse_data.at[int(last_peak_before + (i+1)*used_ibi), 'peaks'] = 1

                        # if the interpolated time window is too large or
                        # the used IBI does not at all match the median IBI, then you should not trust the interpolation
                        # save a text file with this bad interval
                        # later you can remove these time points from the analysis
                        if (time_to_interpolate > 45.0) | (abs(used_ibi - median_ibi) > 0.3 * sr_pulse):
                            bad_interval_pulse[starts[interval]:ends[interval]] = 1
                            print("The interpolation should not be trusted. Will be saved as bad interval in text file.")

                        # plot interpolated result
                        plt.plot(time, raw_data)
                        plt.vlines(np.divide(pulse_data.loc[starts_plot[interval]:ends_plot[interval]].index[pulse_data.loc[starts_plot[interval]:ends_plot[interval]].peaks == 1], sr_pulse),
                                   ymin=min(raw_data), ymax=max(raw_data), color='red')
                        plt.show(block=True)
                        check = True

                    elif todo == 'ignore':
                        check = True  # take data as it is
                    else:
                        print('Please enter one of the given terms!')

        print("Loading respiration data")
        if isub < 17:
            resp_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) +
                                  '*_recording-respirationskinconductance_physio.tsv.gz')
            json_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) +
                                  '*_recording-respirationskinconductance_physio.json')
        else:
            # here we recorded with Biopac both pulse and respiration
            # but have a Brainamp recording with only resp as backup
            resp_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) +
                                  '*_recording-pulsoxyrespiration_physio.tsv.gz')
            json_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) +
                                  '*_recording-pulsoxyrespiration_physio.json')
            alternative_resp_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) +
                                              '*_recording-respiration_physio.tsv.gz')
            alternative_json_file = glob.glob(func_path + os.sep + '*' + task + '*' + str(run) +
                                              '*_recording-respiration_physio.json')

        if len(resp_file) == 1:
            resp_file = resp_file[0]
            json_file = json_file[0]
        else:
            raise Exception('Too many or too few respiration files for this run: ' + str(resp_file))

        # plot respiration data and see whether anything is off
        with open(json_file, 'r') as f:
            json_dict = json.load(f)
        sr_resp = json_dict['SamplingFrequency']
        start_time_resp = json_dict['StartTime']
        column_names = json_dict['Columns']
        resp_data = pd.read_csv(resp_file, header=None, sep='\t', compression='gzip', names=column_names)
        # here data needs to be inverted
        if isub < 17:
            resp_data.respiratory = resp_data.respiratory * (-1)
        # here data is missing because recording was stopped too early, fill with mean value
        if ((sub == 'sub-sspr01') & (run == 2)) | ((sub == 'sub-sspr03') & (run == 4)):
            bad_interval_resp = np.where(resp_data['respiratory'].notna(), 0, 1)
            print("You had missing respiratory data for " + str(sum(bad_interval_resp) / sr_resp) + " s")
            resp_data['respiratory'].fillna(value=np.mean(resp_data['respiratory']), inplace=True)

        plt.plot(resp_data.respiratory)
        # detrend the data
        resp_data.respiratory = signal.detrend(resp_data.respiratory)
        plt.plot(resp_data.respiratory)
        plt.show(block=True)

        # for the first 16 subjects you anyway have to take that data
        # for the others you can check the alternative if this data is off
        if isub < 17:
            take_first_resp = True
        else:
            check = False
            while not check:
                take_this_resp = input('Do you want to use the shown respiration data (y/n)?')
                if take_this_resp == 'y':
                    take_first_resp = True
                    check = True
                # if not show alternative data
                elif take_this_resp == 'n':
                    take_first_resp = False
                    if len(alternative_resp_file) == 1:
                        alternative_resp_file = alternative_resp_file[0]
                        alternative_json_file = alternative_json_file[0]
                    else:
                        raise Exception('Too many or too few respiration files for this run: ' + str(alternative_resp_file))
                    with open(alternative_json_file, 'r') as f:
                        json_dict = json.load(f)
                    sr_resp_alt = json_dict['SamplingFrequency']
                    column_names_resp_alt = json_dict['Columns']
                    start_time_resp_alt = json_dict['StartTime']

                    resp_data_alt = pd.read_csv(alternative_resp_file, header=None, sep='\t', compression='gzip',
                                                names=column_names_resp_alt)
                    # again this is a Brainamp recording that you need to invert first
                    resp_data_alt.respiratory = resp_data_alt.respiratory * (-1)
                    plt.plot(resp_data_alt.respiratory)
                    plt.show(block=True)

                    check_inner = False
                    while not check_inner:
                        take_that_resp = input('Do you want to use the alternative respiration data instead (y/n)?')
                        if take_that_resp == 'y':
                            take_first_resp = False
                            check_inner = True
                        elif take_that_resp == 'n':
                            print("Then I assume the alternative data is even worse and will "
                                  "take the previously shown data instead")
                            take_first_resp = True
                            check_inner = True
                        else:
                            print('Please enter y or n!')
                    check = True
                else:
                    print('Please enter y or n!')

        if take_first_resp:
            chosen_resp_data = resp_data
            start_time_resp = start_time_resp
            sr_resp = sr_resp
        else:
            chosen_resp_data = resp_data_alt
            start_time_resp = start_time_resp_alt
            sr_resp = sr_resp_alt

        if (start_time_pulse == 0) & (start_time_resp == 0):
            # check that both data have the same triggers and TRs
            # in these subjects recordings stopped too early, fill with mean value
            # and save as bad interval
            if (sub == 'sub-sspr01') & (run == 2):
                print('Here we know that something is not matching')
            else:
                assert sum(chosen_resp_data['trigger']) == sum(pulse_data['trigger']), "Trigger numbers not matching!"
                assert (len(chosen_resp_data) / sr_resp / sum(chosen_resp_data['trigger']) ==
                        len(pulse_data) / sr_pulse / sum(pulse_data['trigger'])), "TRs or data lengths are not matching"
        else:
            # in sub 18 run 4 I started the biopac too late, interpolate
            if start_time_pulse < 0:
                print('Interpolating missing data')
                assert (len(chosen_resp_data) / sr_resp / sum(chosen_resp_data['trigger']) ==
                        len(pulse_data) / sr_pulse / sum(pulse_data['trigger'])), "TRs or data lengths are not matching"
                TR = len(chosen_resp_data) / sr_resp / sum(chosen_resp_data['trigger'])  # in s
                beats = np.where(pulse_data.peaks == 1)[0]
                ibis = np.diff(beats)
                median_ibi = statistics.median(ibis)
                print('Median IBI is: ' + str(median_ibi / sr_pulse))
                # fill data with Nan to the required length and then assume that there is a bad interval
                # from beginning until the first proper peak, interpolate as above
                assert (sum(chosen_resp_data['trigger']) - sum(pulse_data['trigger'])) * TR == -start_time_pulse, \
                    "Start time is weird"
                pulse_data_extended = np.hstack([np.empty(int(-start_time_pulse*sr_pulse)), pulse_data['cardiac']])
                peaks_extended = np.hstack([np.zeros(int(-start_time_pulse*sr_pulse)), pulse_data['peaks']])
                trigger_extended = np.hstack([np.zeros(int(-start_time_pulse * sr_pulse)), pulse_data['trigger']])
                trigger_missing = int(-start_time_pulse/TR)
                for t in range(trigger_missing + 1):
                    trigger_extended[int(TR*sr_pulse*t)] = 1
                assert sum(chosen_resp_data['trigger']) == sum(trigger_extended), "Trigger numbers not matching!"
                starts = 0
                first_peak_after = int(np.where(peaks_extended == 1)[0][0])
                last_peak_before = 0

                # how much time is missing and how much peaks can fit in there
                time_to_interpolate = (first_peak_after - last_peak_before) / sr_pulse
                print('Time to interpolate: ' + str(time_to_interpolate))
                peak_count = round((first_peak_after - last_peak_before) / median_ibi - 1)
                # -1 to get to number of peaks instead of number of ibis
                print('No of peaks to interpolate: ' + str(peak_count))
                # start from first peak that is there and go backwards
                # with median IBI
                for i in range(peak_count):
                    peaks_extended[int(first_peak_after - (i + 1) * median_ibi)] = 1
                pulse_data = pd.DataFrame([])
                pulse_data['peaks'] = peaks_extended
                pulse_data['trigger'] = trigger_extended
                pulse_data['cardiac'] = pulse_data_extended

        # normalize respiration data so that the mean is 0
        chosen_resp_data['respiratory'] = stats.zscore(chosen_resp_data['respiratory'])

        # if the sampling rates do not match, downsample the one with the higher SR
        if sr_resp != sr_pulse:
            if sr_resp > sr_pulse:
                resp_down = signal.resample(chosen_resp_data['respiratory'], len(pulse_data))
                pnm = pd.DataFrame({'trigger': pulse_data['trigger'],
                                    'cardiac': pulse_data['peaks'],
                                    'respiratory': resp_down})
                final_sr = sr_pulse
            else:
                pulse_down = signal.resample(pulse_data['cardiac'], len(chosen_resp_data))
                peaks = np.where(pulse_data['peaks'] == 1)[0]
                peaks_down = peaks / sr_pulse * sr_resp
                # x = np.arange(0, len(pulse_down))
                # xticks = np.divide(x, sr_pulse)
                # plt.plot(xticks, pulse_down)
                peaks_down_event_s = np.divide(peaks_down, sr_pulse)
                # plt.scatter(peaks_down_event_s, pulse_down[np.round(peaks_down).astype(int)],
                # s=80, facecolors='none', edgecolors='r')
                # plt.show(block=True)
                peaks_column = np.array([0] * len(chosen_resp_data))
                peaks_column[np.round(peaks_down).astype(int)] = 1
                pnm = pd.DataFrame({'trigger': chosen_resp_data['trigger'],
                                    'cardiac': peaks_column,
                                    'respiratory': chosen_resp_data['respiratory']})
                final_sr = sr_resp
        # otherwise you can just put them together
        else:
            pnm = pd.DataFrame({'trigger': pulse_data['trigger'],
                                'cardiac': pulse_data['peaks'],
                                'respiratory': chosen_resp_data['respiratory']})
            final_sr = sr_pulse

        pnm.to_csv(save_path + os.sep + 'physio.txt', sep=' ', header=False, index=False)
        with open(save_path + os.sep + 'physio_sr.txt', 'w') as f:
            f.write(str(final_sr))

        # bad intervals for respiration were all not relevant, so save only pulse
        # if there is one
        if sum(bad_interval_pulse) > 0:
            bad_interval_pulse = pd.DataFrame(bad_interval_pulse)
            bad_interval_pulse.to_csv(save_path + os.sep + 'physio_bad.txt', header=False, index=False)
            print("Saved bad interval in text file")
