#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Description:
------------
Helper function to display ECG data for R peak correction
Opens a GUI in which EKG raw data and previously detected ECG events (= R peaks) 
will be shown in an interval of a user-specific size (recommended: 20s, should be
specified in samples)
Users can add and delete peaks and mark time periods with bad signal as bad interval


Usage:
------
from class_GUI_hb import r_peak_GUI
my_GUI = r_peak_GUI(raw_data, r_peaks, interval_size, sampling_rate, save_name)


Authors: 
--------
Alexandra John translated Matlab code by Ulrike Horn into Python
(supervised by Ulrike Horn)
Code edited by Lisa-Marie Pohle to include bad intervals

"""

import math
import tkinter as tk
from tkinter import *
from matplotlib.figure import Figure
import numpy as np
from matplotlib.offsetbox import AnchoredText
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import pandas as pd


class r_peak_GUI:

    def __init__(self, raw_data, ecg_events, bad_interval_starts, bad_interval_ends, interval_size, sampling_frequency, save_path):
        self.root = tk.Tk()
        self.idx = 0
        self.interval_start = 0
        self.interval_size = interval_size
        self.interval_end = self.interval_start + self.interval_size
        self.idx_max = math.ceil((len(raw_data[:])) / self.interval_size) - 1
        self.click_collector = []
        self.ecg_events = ecg_events
        self.bad_interval_starts = bad_interval_starts
        self.bad_interval_ends = bad_interval_ends
        self.sampling_frequency = sampling_frequency
        self.raw_data = raw_data
        self.save_path = save_path
        self.root.title("R-PEAK DETECTION")
        self.root.geometry("1000x600")

        self.root.columnconfigure(0, pad=3, weight=1)
        self.root.columnconfigure(1, pad=3, weight=1)
        self.root.columnconfigure(2, pad=3, weight=1)
        self.root.columnconfigure(3, pad=3, weight=1)
        self.root.columnconfigure(4, pad=3, weight=1)
        self.root.rowconfigure(0, pad=3, weight=20)
        self.root.rowconfigure(1, pad=3, weight=1)
        self.root.rowconfigure(2, pad=3, weight=1)

        # matplotlib plot graph
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        x = np.arange(self.interval_start, self.interval_end)
        xticks = np.divide(x, self.sampling_frequency)  # to convert sampling rate in sec
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # 0.499 1.102 s
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        [self.ylim_min, self.ylim_max] = self.ax.get_ylim()

        # connect matplotlib figure to tkinter window
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)  # A tk.DrawingArea.
        self.canvas.draw()
        # self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.get_tk_widget().grid(column=0, columnspan=5, row=0, sticky='nsew')

        # connect mouse_click to fig
        self.canvas.mpl_connect("button_press_event", self.mouse_click)

        # button close
        # master.destroy
        self.button1 = tk.Button(self.root, text="close", command=lambda: self.closing())
        self.button1.grid(row=3, column=2)
        # self.button1.pack(side="bottom", pady='10')

        # button next
        self.button2 = tk.Button(self.root, text="next", command=lambda: self.next_timeframe())
        self.button2.grid(row=2, column=4)
        # self.button2.pack(side="right")

        # button back
        self.button3 = tk.Button(self.root, text="back", command=lambda: self.previous_timeframe())
        self.button3.grid(row=2, column=0)
        # self.button3.pack(side="left")

        # button delete
        self.button4 = tk.Button(self.root, text="delete peak", command=lambda: self.delete_r_peaks())
        self.button4.grid(row=2, column=2)
        # self.button4.pack(side="bottom")

        # button add
        self.button5 = tk.Button(self.root, text="add R-peak", command=lambda: self.add_r_peak())
        self.button5.grid(row=1, column=2)
        # self.button5.pack(side="bottom")

        # button add starting point of bad interval
        self.button6 = tk.Button(self.root, text="add start of bad interval", command=lambda: self.add_start_bad())
        self.button6.grid(row=1, column=1)
        # self.button6.pack(side='left', padx='10')

        # button add ending point of bad interval
        self.button7 = tk.Button(self.root, text="add end of bad interval", command=lambda: self.add_end_bad())
        self.button7.grid(row=1, column=3)
        # self.button7.pack(side='bottom')

        # button remove starting point of bad interval
        self.button8 = tk.Button(self.root, text="remove start of bad interval",
                                 command=lambda: self.delete_start_bad())
        self.button8.grid(row=2, column=1)
        # self.button8.pack(side='left', padx = '10')

        # button remove ending  point of bad interval
        self.button9 = tk.Button(self.root, text="remove end of bad interval", command=lambda: self.delete_end_bad())
        self.button9.grid(row=2, column=3)
        # self.button8.pack(side='bottom')

        # define what happens when clicking on closing icon
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.root.mainloop()

    def closing(self):

        if len(self.bad_interval_starts) == len(self.bad_interval_ends):
            output_data = pd.DataFrame({'cardiac': self.raw_data, 'peaks': 0, 'bad_intervals_1start_2end': 0})
            output_data.loc[self.ecg_events, 'peaks'] = 1
            output_data.loc[self.bad_interval_starts, 'bad_intervals_1start_2end'] = 1
            output_data.loc[self.bad_interval_ends, 'bad_intervals_1start_2end'] = 2
            output_data.to_csv(self.save_path, index=False)
            self.root.destroy()
        else:
            tk.messagebox.showwarning(title='Attention!', message='The number of starting points of bad intervals is '
                                                                  'unequal to the number of endpoints. Please check '
                                                                  'the data again!')

    def mouse_click(self, event):

        global x_coord
        x_coord = event.xdata
        if x_coord:
            print("selected x coordinate:", x_coord.round(3))

            # collecting data points -> it only keeps 1 data point
            self.click_collector.append(x_coord)
            if len(self.click_collector) > 1:
                del self.click_collector[0]

            # update plot
            self.fig.clear()
            x = np.arange(self.interval_start, self.interval_end)  # 0 1 2 3
            xticks = np.divide(x, self.sampling_frequency)  # 0 0.001 0.002 s
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
            x_events = self.ecg_events[
                (self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]  # 499 1102
            x_events_s = np.divide(x_events, self.sampling_frequency)  # 0.499 1.102 s
            self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
            x_bads_start = self.bad_interval_starts[
                (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
            x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
            x_bads_end = self.bad_interval_ends[
                (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
            x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
            self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
            self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
            for x in x_bads_start_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
            for x in x_bads_end_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
            self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
            self.ax.set_ylim(self.ylim_min, self.ylim_max)
            self.canvas.draw()

            for i, idx in enumerate(self.click_collector[:]):
                # plot a red dot
                self.ax.scatter(self.click_collector[i],
                                self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red',
                                marker='+')
                self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
                self.canvas.draw()

    def next_timeframe(self):

        # increase self.idx index and adapt interval_start and interval_end
        if self.idx < self.idx_max - 1:
            self.idx = self.idx + 1
            self.interval_start = self.idx * self.interval_size
            self.interval_end = self.interval_start + self.interval_size
            # plot next time interval
            self.fig.clf()
            x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
            xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
            x_events = self.ecg_events[
                (self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
            x_events_s = np.divide(x_events, self.sampling_frequency)  # to convert sampling rate in sec
            self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
            x_bads_start = self.bad_interval_starts[
                (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
            x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
            x_bads_end = self.bad_interval_ends[
                (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
            x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
            self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
            self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
            for x in x_bads_start_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
            for x in x_bads_end_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
            self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
            [self.ylim_min, self.ylim_max] = self.ax.get_ylim()
            self.canvas.draw()
        elif self.idx < self.idx_max:  # last interval
            self.idx = self.idx + 1
            self.interval_start = self.idx * self.interval_size
            self.interval_end = len(self.raw_data[:])
            self.fig.clf()
            x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
            xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
            x_events = self.ecg_events[
                (self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
            x_events_s = np.divide(x_events, self.sampling_frequency)  # to convert sampling rate in sec
            self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
            x_bads_start = self.bad_interval_starts[
                (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
            x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
            x_bads_end = self.bad_interval_ends[
                (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
            x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
            self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
            self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
            for x in x_bads_start_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
            for x in x_bads_end_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
            self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
            [self.ylim_min, self.ylim_max] = self.ax.get_ylim()
            self.canvas.draw()

    def previous_timeframe(self):

        # decrease self.idx index and adapt interval_start and interval_end
        if self.idx > 0:  # to make sure that going back is limited
            self.idx = self.idx - 1
            self.interval_start = self.idx * self.interval_size
            self.interval_end = self.interval_start + self.interval_size
            # plot previous time interval
            self.fig.clf()
            x = np.arange(self.interval_start, self.interval_end)
            xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
            x_events = self.ecg_events[
                (self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
            x_events_s = np.divide(x_events, self.sampling_frequency)  # to convert sampling rate in sec
            self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
            x_bads_start = self.bad_interval_starts[
                (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
            x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
            x_bads_end = self.bad_interval_ends[
                (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
            x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
            self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
            self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
            for x in x_bads_start_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
            for x in x_bads_end_s:
                self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
            self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
            [self.ylim_min, self.ylim_max] = self.ax.get_ylim()
            self.canvas.draw()
            self.canvas.get_tk_widget().grid(column=0, columnspan=5, row=0, sticky='nsew')

    def add_r_peak(self):

        # takes value from mouseclick (click_collector -> in sec)
        value = int(self.click_collector[0] * self.sampling_frequency)  # convert sec in sampling rate
        # search for nearest max around y_value
        range_max = int(0.2 * self.sampling_frequency)  # interval to search for maximum peak
        # the search interval can be shorter if we are at the beginning or end of a data set
        if value - range_max < 0:
            new_x_max_y = np.argmax(self.raw_data[0:value + range_max])
            value_in_range = value
        elif value + range_max > len(self.raw_data):
            new_x_max_y = np.argmax(self.raw_data[value - range_max:len(self.raw_data)])
            value_in_range = range_max
        else:
            new_x_max_y = np.argmax(self.raw_data[value - range_max:value + range_max])
            value_in_range = range_max

        # calculated maximum in relation to whole data
        new_x_max_y = (new_x_max_y - value_in_range) + value

        # add value to ecg_events array
        self.ecg_events = np.append(self.ecg_events, new_x_max_y)

        # update plot
        self.fig.clear()
        x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
        xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # to sec
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.canvas.draw()

        for i, idx in enumerate(self.click_collector[:]):
            # plot a red dot
            self.ax.scatter(self.click_collector[i],
                            self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red', marker='+')
            self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
            self.canvas.draw()

        return self.ecg_events

    def delete_r_peaks(self):

        # takes value from mouseclick (click_collector -> in sec)
        value = int(self.click_collector[0] * self.sampling_frequency)  # convert sec in sampling rate

        print("value:", value)
        # finds nearest value in ecg_events
        nearest_value = self.ecg_events[np.abs(self.ecg_events - value).argmin()]
        print("nearest:", nearest_value)

        # search for index of this value and then remove it
        idx_remove = np.where(
            self.ecg_events == nearest_value)  # saves 1. array containing rod idxand 2. with column idx
        idx_remove = idx_remove[0]  # just need row idx array
        print("idx_remove", idx_remove)
        self.ecg_events = np.delete(self.ecg_events, idx_remove)

        # update plot
        self.fig.clear()
        x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
        xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # to convert sampling rate in sec
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.canvas.draw()

        for i, idx in enumerate(self.click_collector[:]):
            # plot a red dot
            self.ax.scatter(self.click_collector[i],
                            self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red', marker='+')
            self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
            self.canvas.draw()

        return self.ecg_events

    def add_start_bad(self):

        # takes value from mouseclick (click_collector -> in sec)
        value = int(self.click_collector[0] * self.sampling_frequency)  # convert sec in sampling rate

        # add value to bad_interval_starts array
        self.bad_interval_starts = np.append(self.bad_interval_starts, value)

        # update plot
        self.fig.clear()
        x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
        xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # to sec
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.canvas.draw()

        for i, idx in enumerate(self.click_collector[:]):
            # plot a red dot
            self.ax.scatter(self.click_collector[i],
                            self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red', marker='+')
            self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
            self.canvas.draw()

        return self.bad_interval_starts

    def add_end_bad(self):

        # takes value from mouseclick (click_collector -> in sec)
        value = int(self.click_collector[0] * self.sampling_frequency)  # convert sec in sampling rate

        # add value to bad_interval_ends array
        self.bad_interval_ends = np.append(self.bad_interval_ends, value)

        # update plot
        self.fig.clear()
        x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
        xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # to sec
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.canvas.draw()

        for i, idx in enumerate(self.click_collector[:]):
            # plot a red dot
            self.ax.scatter(self.click_collector[i],
                            self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red', marker='+')
            self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
            self.canvas.draw()

        return self.bad_interval_ends

    def delete_start_bad(self):

        # takes value from mouseclick (click_collector -> in sec)
        value = int(self.click_collector[0] * self.sampling_frequency)  # convert sec in sampling rate

        print("value:", value)
        # finds nearest value in ecg_events
        nearest_value = self.bad_interval_starts[np.abs(self.bad_interval_starts - value).argmin()]
        print("nearest:", nearest_value)

        # search for index of this value and then remove it
        idx_remove = np.where(
            self.bad_interval_starts == nearest_value)  # saves 1. array containing rod idxand 2. with column idx
        idx_remove = idx_remove[0]  # just need row idx array
        print("idx_remove", idx_remove)
        self.bad_interval_starts = np.delete(self.bad_interval_starts, idx_remove)

        # update plot
        self.fig.clear()
        x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
        xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # to convert sampling rate in sec
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.canvas.draw()

        for i, idx in enumerate(self.click_collector[:]):
            # plot a red dot
            self.ax.scatter(self.click_collector[i],
                            self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red', marker='+')
            self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
            self.canvas.draw()

        return self.bad_interval_starts

    def delete_end_bad(self):

        # takes value from mouseclick (click_collector -> in sec)
        value = int(self.click_collector[0] * self.sampling_frequency)  # convert sec in sampling rate

        print("value:", value)
        # finds nearest value in ecg_events
        nearest_value = self.bad_interval_ends[np.abs(self.bad_interval_ends - value).argmin()]
        print("nearest:", nearest_value)

        # search for index of this value and then remove it
        idx_remove = np.where(
            self.bad_interval_ends == nearest_value)  # saves 1. array containing rod idxand 2. with column idx
        idx_remove = idx_remove[0]  # just need row idx array
        print("idx_remove", idx_remove)
        self.bad_interval_ends = np.delete(self.bad_interval_ends, idx_remove)

        # update plot
        self.fig.clear()
        x = np.arange(self.interval_start, self.interval_end)  # to make ongoing x axis
        xticks = np.divide(x, self.sampling_frequency)  # to convert in sec
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(xticks, self.raw_data[self.interval_start:self.interval_end])
        x_events = self.ecg_events[(self.interval_start <= self.ecg_events) & (self.ecg_events <= self.interval_end)]
        x_events_s = np.divide(x_events, self.sampling_frequency)  # to convert sampling rate in sec
        self.ax.scatter(x_events_s, self.raw_data[x_events], s=80, facecolors='none', edgecolors='r')
        x_bads_start = self.bad_interval_starts[
            (self.interval_start <= self.bad_interval_starts) & (self.bad_interval_starts <= self.interval_end)]
        x_bads_start_s = np.divide(x_bads_start, self.sampling_frequency)
        x_bads_end = self.bad_interval_ends[
            (self.interval_start <= self.bad_interval_ends) & (self.bad_interval_ends <= self.interval_end)]
        x_bads_end_s = np.divide(x_bads_end, self.sampling_frequency)
        self.ax.vlines(x_bads_start_s, min(self.raw_data), max(self.raw_data), colors='blue')
        self.ax.vlines(x_bads_end_s, min(self.raw_data), max(self.raw_data), colors='blue')
        for x in x_bads_start_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 0.5, 'Start Bad Interval', color='blue')
        for x in x_bads_end_s:
            self.ax.text(x + 0.1, max(self.raw_data) - 1, 'End Bad Interval', color='blue')
        self.ax.add_artist(AnchoredText("[{}/{}]".format(self.idx + 1, self.idx_max + 1), 'lower right'))
        self.ax.set_ylim(self.ylim_min, self.ylim_max)
        self.canvas.draw()

        for i, idx in enumerate(self.click_collector[:]):
            # plot a red dot
            self.ax.scatter(self.click_collector[i],
                            self.raw_data[int(self.click_collector[i] * self.sampling_frequency)], c='red', marker='+')
            self.ax.axvline(x=self.click_collector[i], c='red', alpha=0.4)
            self.canvas.draw()

        return self.bad_interval_ends

    def on_closing(self):
        if len(self.bad_interval_starts) == len(self.bad_interval_ends):
            close_action = tk.messagebox.askyesnocancel("Quit", "Do you want to save your data on closing?")
            if close_action:
                output_data = pd.DataFrame({'cardiac': self.raw_data, 'peaks': 0, 'bad_intervals_1start_2end': 0})
                output_data.loc[self.ecg_events, 'peaks'] = 1
                output_data.loc[self.bad_interval_starts, 'bad_intervals_1start_2end'] = 1
                output_data.loc[self.bad_interval_ends, 'bad_intervals_1start_2end'] = 2
                output_data.to_csv(self.save_path, index=False)
                self.root.destroy()
            elif close_action == None:
                pass
            else:
                self.root.destroy()
        else:
            tk.messagebox.showwarning(title='Attention!',
                                      message='The number of starting points of bad intervals is unequal to the number of endpoints. Please check the data again!')
