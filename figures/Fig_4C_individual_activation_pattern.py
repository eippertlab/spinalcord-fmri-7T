import os
import numpy as np
import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

data_path = '/data/pt_02661_raw/Heatpain/derivatives'

subjects = np.arange(1, 42)
# delete subjects that have been measured twice
subjects = np.delete(subjects, 11)
subjects = np.delete(subjects, 0)

region = 'dh_left'

collect_data = True
make_plot = True

# C5 ends at z=882 and C8 starts at z=756 --> take only these
min_z = 755
max_z = 883
C6 = [819, 849]
C7 = [788, 818]

if collect_data:
    df_all = pd.DataFrame()
    for isub in subjects:
        sub_id = 'sub-sspr' + '{:02d}'.format(isub)
        print(sub_id)
        fname = os.path.join(data_path, sub_id, 'copes_templatespace',
                             'block_design_2mm_susan', 'extracted', region + '_cope1_all.txt')
        df = pd.read_csv(fname, header=None, sep='\s+').transpose()
        df = df.rename(columns={0: "x", 1: "y", 2: "z", 3: "value"})
        avgs = df.groupby('z').mean().reset_index()
        avgs.drop(['x', 'y'], axis=1, inplace=True)

        selection = avgs[(avgs['z'] > min_z) & (avgs['z'] < max_z)].copy()
        selection['sub'] = isub

        # z-score data within a subject
        # selection['value'] = stats.zscore(selection['value'])

        # if necessary smooth it a bit more
        # from scipy.ndimage import gaussian_filter1d
        # selection['smoothed_beta'] = np.transpose(gaussian_filter1d([selection['value']], sigma=1.5))

        df_all = pd.concat([df_all, selection])

    # set all negative betas to 0
    # df_all.loc[df_all['value'] < 0, 'value'] = 0

    # create array
    pivot_df = df_all.pivot(index='z', columns='sub', values='value')
    my_array = pivot_df.to_numpy()

    pd.to_pickle(my_array, os.path.join(data_path, 'results', 'heat_map_data_raw_beta_{}.pkl'.format(region)))

if make_plot:
    my_array = pd.read_pickle(os.path.join(data_path, 'results', 'heat_map_data_raw_beta_{}.pkl'.format(region)))
    segment_df = pd.DataFrame({'label': ['C5', 'C6', 'C7', 'C8'],
                               'start': [850, 819, 788, 756],
                               'end': [882, 849, 818, 787]})

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap
    sns.set_style('ticks')
    sns.set_context("talk")
    f, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(my_array, cmap=cmap, vmin=-150, vmax=150, center=0, 
                mask=(my_array == 0),
                linewidths=.5, xticklabels=np.arange(1, len(subjects) + 1),
                yticklabels=False, cbar_kws={'label': 'beta'})
    middle_points = []
    for i, segment in segment_df.iterrows():
        plt.axhline(y=segment['start'] - min_z, color='black', linestyle='--')
        middle_points.append(segment['start'] - min_z + (segment['end'] - segment['start'])/2)
    plt.axhline(y=segment_df.loc[segment_df['label'] == 'C5', 'end'][0] - min_z, color='black', linestyle='--')
    ax.set_ylim([max_z - min_z, -1])
    ax.set_yticks(middle_points, labels=segment_df['label'])
    ax.invert_yaxis()
    ax.set_xlabel('Participant')
    ax.set_ylabel('Segment')
    ax.set_xticklabels([])
    f.savefig(os.path.join(data_path, 'results', 'heat_map_' + region + '_single_subjects_raw_beta.svg'),
              transparent=True, format='svg', bbox_inches='tight')
    plt.show()
