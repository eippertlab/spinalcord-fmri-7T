"""
Run the accompanying bash file first!
"""

import numpy as np
import os
import pandas as pd
import scipy.io as sio
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt
import seaborn as sns

main_path = '/data/pt_02661_raw/Heatpain/derivatives'
data_path = os.path.join(main_path, 'zshim_calc')

if not os.path.exists(data_path):
    os.makedirs(data_path)

color_green_bright = [161/255, 217/255, 155/255]
color_green = [49/255, 163/255, 84/255]

# which steps to run
collect_data = False
statistics = True
plot = False

# dataset 3
subjects = np.arange(17, 42)
numSlices = 16

if collect_data:
    # collect all subjects' slice wise values for Y = Z-shim, N=No z-shim
    subjValsY = np.zeros((numSlices, len(subjects)))
    subjValsN = np.zeros((numSlices, len(subjects)))
    sub_counter = 0
    for isub in subjects:
        sub_id = 'sub-sspr' + '{:02d}'.format(isub)
        print(f"Working on subject {sub_id}")
        zshim_path = os.path.join(main_path, sub_id, 'func', 'zshim')
        # Load data
        dataY = np.loadtxt(os.path.join(zshim_path, 'zshim_recon_signal_all.txt'))[2:4, :]
        dataN = np.loadtxt(os.path.join(zshim_path, 'no_zshim_signal_all.txt'))[2:4, :]
        sliceValsY = np.zeros(16)
        sliceValsN = np.zeros(16)
        for s in range(0, numSlices):
            firstY = np.where(dataY[0, :] == s)[0][0]
            firstN = np.where(dataN[0, :] == s)[0][0]
            lastY = np.where(dataY[0, :] == s)[0][-1]
            lastN = np.where(dataN[0, :] == s)[0][-1]
            sliceValsY[s] = np.mean(dataY[1, firstY:lastY+1])
            sliceValsN[s] = np.mean(dataN[1, firstN:lastN+1])
        subjValsY[:, sub_counter] = sliceValsY
        subjValsN[:, sub_counter] = sliceValsN
        sub_counter += 1
    np.save(os.path.join(data_path, 'valuesY'), subjValsY)
    np.save(os.path.join(data_path, 'valuesN'), subjValsN)

if statistics:
    # Load the data
    valuesY = np.load(os.path.join(data_path, 'valuesY.npy'))
    valuesN = np.load(os.path.join(data_path, 'valuesN.npy'))

    # Calculate percent of average intensity increase (one value / subject)
    groupValsM = np.column_stack((np.mean(valuesN, axis=0), np.mean(valuesY, axis=0)))
    percInc = groupValsM[:, 1] / groupValsM[:, 0] - 1
    resPercInc = [np.mean(percInc), np.min(percInc), np.max(percInc)]
    print('The increases are on average {}, from minimum {} to maximum {} increase'.format(round(resPercInc[0], 3),
                                                                                           round(resPercInc[1], 3),
                                                                                           round(resPercInc[2], 3)))
    hM, pM = ttest_rel(groupValsM[:, 1], groupValsM[:, 0])
    print('The hypothesis that the average intensities are the same ' + \
          'is rejected with a p-value of {} and a t-value of {}'.format(pM, hM))
    differences = groupValsM[:, 1] - groupValsM[:, 0]
    mean_diff = np.mean(differences)
    std_dev = np.std(differences, ddof=1)
    cohen_d_value = mean_diff/std_dev
    print('Cohens D: {}'.format(cohen_d_value))
    meanNoM, meanYesM = groupValsM[:, 0], groupValsM[:, 1]

    # Calculate percent of signal variation decrease across slices
    groupValsV = np.column_stack((np.var(valuesN, axis=0, ddof=1), np.var(valuesY, axis=0, ddof=1)))
    percDec = 1 - groupValsV[:, 1] / groupValsV[:, 0]
    resPercDec = [np.mean(percDec), np.min(percDec), np.max(percDec)]
    print('The decreases are on average {}, from minimum {} to maximum {} decrease'.format(round(resPercDec[0], 3),
                                                                                           round(resPercDec[1], 3),
                                                                                           round(resPercDec[2], 3)))
    hV, pV = ttest_rel(groupValsV[:, 0], groupValsV[:, 1])
    print('The hypothesis that the intensity variations are the same ' + \
          'is rejected with a p-value of {} and a t-value of {}'.format(pV, hV))
    differences = groupValsV[:, 0] - groupValsV[:, 1]
    mean_diff = np.mean(differences)
    std_dev = np.std(differences, ddof=1)
    cohen_d_value = mean_diff/std_dev
    print('Cohens D: {}'.format(cohen_d_value))
    meanNoV, meanYesV = groupValsV[:, 0], groupValsV[:, 1]

    # Save the results matrix
    results = pd.DataFrame({'No_M': meanNoM, 'Yes_M': meanYesM, 'No_V': meanNoV, 'Yes_V': meanYesV})
    results.to_csv(os.path.join(data_path, 'resultsMatrix.csv'), index=None)

if plot:
    sns.set_context("poster")
    sns.set_style("ticks")
    jitter = 0.1
    # order: violinplot, boxplot, raw data, raw data, boxplot, violinplot
    positions = [0.5, 0.75, 1.25, 2.25, 2.75, 3]
    width_violin = 0.8
    width_box = 0.2

    df = pd.read_csv(os.path.join(data_path, 'resultsMatrix.csv'), index_col=None)
    df.reset_index(inplace=True)
    df_long = pd.wide_to_long(df, stubnames=['No', 'Yes'], i='index', j='Stat', sep='_', suffix='\\w+').reset_index()

    for stats in ['M', 'V']:
        curr_df = df_long[df_long['Stat'] == stats]
        curr_df.reset_index(inplace=True, drop=True)

        df_x_jitter = pd.DataFrame(np.random.normal(loc=positions[2], scale=jitter, size=(curr_df.values.shape[0], 2)),
                                   columns=['nozshim', 'zshim'])
        df_x_jitter['zshim'] += positions[3] - positions[2]

        fig, ax = plt.subplots(1, 1, figsize=(5, 10))
        # if stats == 'M':
        #     fig.suptitle("Average signal intensity", fontsize=30)
        # else:
        #     fig.suptitle("Signal intensity variation", fontsize=30)

        v1 = ax.violinplot(curr_df['No'], positions=[positions[0]], showmeans=False, showextrema=False,
                           showmedians=False, side="low", widths=width_violin)  # `side` param requires matplotlib 3.9+
        for pc in v1['bodies']:
            pc.set_facecolor(color_green_bright)
            pc.set_edgecolor('black')
            pc.set_linewidth(1.2)
            pc.set_alpha(0.7)

        a = sns.boxplot(ax=ax, x=np.repeat(positions[1], curr_df.values.shape[0]),
                        y=curr_df['No'], color=color_green_bright, width=width_box, native_scale=True, showfliers=False)
        for patch in a.patches:
            r, g, b, al = patch.get_facecolor()
            patch.set_facecolor((r, g, b, 0.7))

        b = sns.boxplot(ax=ax, x=np.repeat(positions[4], curr_df.values.shape[0]),
                        y=curr_df['Yes'], color=color_green, width=width_box, native_scale=True, showfliers=False)
        for patch in b.patches:
            r, g, b, al = patch.get_facecolor()
            patch.set_facecolor((r, g, b, 0.7))

        v2 = ax.violinplot(curr_df['Yes'], positions=[positions[5]], showmeans=False, showextrema=False,
                           showmedians=False, side="high", widths=width_violin)  # `side` param requires matplotlib 3.9+
        for pc in v2['bodies']:
            pc.set_facecolor(color_green)
            pc.set_edgecolor('black')
            pc.set_linewidth(1.2)
            pc.set_alpha(0.7)

        df_no = pd.DataFrame({'jitter': df_x_jitter['nozshim'],
                              'value': curr_df['No']})
        df_yes = pd.DataFrame({'jitter': df_x_jitter['zshim'],
                               'value': curr_df['Yes']})
        a = sns.scatterplot(ax=ax, data=df_no, x='jitter', y='value', color=color_green_bright,
                            zorder=100, edgecolor="black", alpha=0.6, legend=False)
        b = sns.scatterplot(ax=ax, data=df_yes, x='jitter', y='value', color=color_green,
                            zorder=100, edgecolor="black", alpha=0.6)
        for idx in curr_df.index:
            ax.plot(df_x_jitter.loc[idx, ['nozshim', 'zshim']], curr_df.loc[idx, ['No', 'Yes']], color='gray',
                    linewidth=1.0, linestyle='-', zorder=-1)
        ax.set_xticks([0.5, 3.5])
        ax.set_xticklabels(['No z-shim', 'Z-shim'])
        ax.tick_params(axis='x', which='major')
        # ax.set_ylim(0, 22)
        sns.despine()
        ax.set_xlabel('Acquisition')
        if stats == 'M':
            plt.ylim([200, 900])
            ax.set_ylabel('Mean signal intensity (a.u.)')
        else:
            plt.ylim([0, 18000])
            ax.set_ylabel('Across-slice signal intensity variation (a.u.)')

        plt.subplots_adjust(bottom=0.15, left=0.2)

        fig.savefig(os.path.join(main_path, 'results', 
                                 'Zshim_comparison_{}.svg'.format(stats)), 
                    transparent=True, bbox_inches='tight', format='svg')
        plt.show()
