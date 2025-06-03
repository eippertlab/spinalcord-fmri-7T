import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt
import nibabel as nib

output_path = '/data/pt_02661_raw/Heatpain/derivatives/'

subjects = np.arange(1, 42)

color_lilac = [117/255, 112/255, 179/255]
color_green = [49/255, 163/255, 84/255]

df_all = pd.DataFrame([])
for isub in subjects:
    sub_id = 'sub-sspr' + '{:02d}'.format(isub)
    func_path = output_path + os.sep + sub_id + os.sep + 'func'
    if (isub == 6) | (isub == 25) | (isub == 27) | (isub == 29):
        num_runs = 2
    elif isub == 31:
        num_runs = 3
    else:
        num_runs = 4
    for irun in range(num_runs):
        raw_tsnr_file = os.path.join(func_path, 'run-' + str(irun + 1), 
                                     'bold_tsnr.nii.gz')
        raw_tsnr = nib.load(raw_tsnr_file)
        raw_tsnr_img = raw_tsnr.get_fdata()
        
        denoised_tsnr_file = os.path.join(func_path, 'run-' + str(irun + 1), 
                                          'bold_denoised_tsnr.nii.gz')
        denoised_tsnr = nib.load(denoised_tsnr_file)
        denoised_tsnr_img = denoised_tsnr.get_fdata()
        
        cord_mask = os.path.join(func_path, 'second_level', 
                                 'MOCO_MEAN_deepseg_manual.nii.gz')
        mask = nib.load(cord_mask)
        mask_img = mask.get_fdata()
        
        raw_mean = np.mean(raw_tsnr_img[mask_img == 1])
        denoised_mean = np.mean(denoised_tsnr_img[mask_img == 1])
        
        df = pd.DataFrame([raw_mean, denoised_mean], columns=['tSNR'])
        df['processing'] = ['raw', 'denoised']
        df['subject'] = sub_id
        df['run'] = irun+1
        if isub > 16:
            df['new'] = 1
        else:
            df['new'] = 0
        df_all = pd.concat([df_all, df])

# average over runs
avgs = df_all.groupby(['subject', 'processing']).mean()
avgs.reset_index(inplace=True)
avgs.drop('run', axis=1, inplace=True)

result = avgs.groupby(['processing', 'new']).mean(['tSNR'])
print(result)
denoised_old = result['tSNR']['denoised'][0.0]
denoised_new = result['tSNR']['denoised'][1.0]
raw_old = result['tSNR']['raw'][0.0]
raw_new = result['tSNR']['raw'][1.0]

perc_increase_old = ((denoised_old - raw_old) / raw_old) * 100
print("This equals an increase of {} in dataset 2".format(round(perc_increase_old, 3)))
result = ttest_rel(avgs[(avgs['new']==0) & (avgs['processing'] == 'denoised')]['tSNR'],
                   avgs[(avgs['new']==0) & (avgs['processing'] == 'raw')]['tSNR'])
print(result)
differences = avgs[(avgs['new']==0) & (avgs['processing'] == 'denoised')]['tSNR'].values - avgs[(avgs['new']==0) & (avgs['processing'] == 'raw')]['tSNR'].values
mean_diff = np.mean(differences)
std_dev = np.std(differences, ddof=1)
cohen_d_value = mean_diff/std_dev
print('Cohens D: {}'.format(cohen_d_value))

perc_increase_new = ((denoised_new - raw_new) / raw_new) * 100
print("This equals an increase of {} in dataset 3".format(round(perc_increase_new, 3)))
result = ttest_rel(avgs[(avgs['new']==1) & (avgs['processing'] == 'denoised')]['tSNR'], 
                   avgs[(avgs['new']==1) & (avgs['processing'] == 'raw')]['tSNR'])
print(result)
differences = avgs[(avgs['new']==1) & (avgs['processing'] == 'denoised')]['tSNR'].values - avgs[(avgs['new']==1) & (avgs['processing'] == 'raw')]['tSNR'].values
mean_diff = np.mean(differences)
std_dev = np.std(differences, ddof=1)
cohen_d_value = mean_diff/std_dev
print('Cohens D: {}'.format(cohen_d_value))

sns.set_context("poster")
sns.set_style("ticks")
jitter = 0.1
alpha = 0.7

# order: violinplot, boxplot, boxplot, raw data, raw data, boxplot, boxplot, violinplot
positions = [0.5, 0.75, 1, 1.5, 2.5, 3, 3.25, 2.5]
width_violin = 0.8
width_box = 0.2

pivot_df = avgs.pivot_table(index=['subject', 'new'], columns=['processing'], values='tSNR')
pivot_df.reset_index(inplace=True, drop=False)
df_x_jitter = pd.DataFrame(np.random.normal(loc=positions[3], scale=jitter, size=(pivot_df.values.shape[0], 2)),
                           columns=['raw', 'denoised'])
df_x_jitter['denoised'] += positions[4] - positions[3]

fig, ax = plt.subplots(1, 1, figsize=(8, 10))

v1 = sns.violinplot(pivot_df, y='raw', x=positions[0], order=positions,
                    split=True, inner=None, hue=True, hue_order=[True, False],
                    palette=['gray', 'gray'], linecolor='black', linewidth=1.2, legend=None)

box1 = sns.boxplot(ax=ax, x=np.repeat(positions[1], len(pivot_df[pivot_df['new'] == 0])),
                   y=pivot_df[pivot_df['new'] == 0]['raw'],
                   color=color_lilac, width=width_box,
                   boxprops=dict(edgecolor='black', linewidth=1.2),
                   whiskerprops=dict(color='black', linewidth=1.2),
                   native_scale=True, showfliers=False)
for patch in box1.patches:
    r, g, b, al = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.7))

box2 = sns.boxplot(ax=ax, x=np.repeat(positions[2], len(pivot_df[pivot_df['new'] == 1])),
                   y=pivot_df[pivot_df['new'] == 1]['raw'],
                   color=color_green, width=width_box,
                   boxprops=dict(edgecolor='black', linewidth=1.2),
                   whiskerprops=dict(color='black', linewidth=1.2),
                   native_scale=True, showfliers=False)
for patch in box2.patches:
    r, g, b, al = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.7))

v2 = sns.violinplot(pivot_df, y='denoised', x=positions[7], order=positions,
                    split=True, inner=None, hue=False, hue_order=[True, False],
                    palette=['gray', 'gray'],
                    linecolor='black', linewidth=1.2, legend=None)

box3 = sns.boxplot(ax=ax, x=np.repeat(positions[5], len(pivot_df[pivot_df['new'] == 0])),
                   y=pivot_df[pivot_df['new'] == 0]['denoised'],
                   color=color_lilac, width=width_box,
                   boxprops=dict(edgecolor='black', linewidth=1.2),
                   whiskerprops=dict(color='black', linewidth=1.2),
                   native_scale=True, showfliers=False)
for patch in box3.patches:
    r, g, b, al = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.7))

box4 = sns.boxplot(ax=ax, x=np.repeat(positions[6], len(pivot_df[pivot_df['new'] == 1])),
                   y=pivot_df[pivot_df['new'] == 1]['denoised'],
                   color=color_green, width=width_box,
                   boxprops=dict(edgecolor='black', linewidth=1.2),
                   whiskerprops=dict(color='black', linewidth=1.2),
                   native_scale=True, showfliers=False)
for patch in box4.patches:
    r, g, b, al = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.7))

df_raw = pd.DataFrame({'jitter': df_x_jitter['raw'], 'tsnr': pivot_df['raw'], 'group': pivot_df['new']})
df_denoised = pd.DataFrame({'jitter': df_x_jitter['denoised'], 'tsnr': pivot_df['denoised'], 'group': pivot_df['new']})
a = sns.scatterplot(ax=ax, data=df_raw, x='jitter', y='tsnr', hue='group',
                    palette={0: color_lilac, 1: color_green},
                    zorder=100, edgecolor="black", alpha=0.6, legend=False)
b = sns.scatterplot(ax=ax, data=df_denoised, x='jitter', y='tsnr', hue='group',
                    palette={0: color_lilac, 1: color_green},
                    zorder=100, edgecolor="black", alpha=0.6)
for idx in pivot_df.index:
    ax.plot(df_x_jitter.loc[idx, ['raw', 'denoised']], pivot_df.loc[idx, ['raw', 'denoised']], color='gray',
            linewidth=1.0, linestyle='-', zorder=-1)
ax.set_xticks([0.5, 3.5])
ax.set_xticklabels(['raw', 'denoised'])
ax.set_ylim(0, 22)
sns.despine()
ax.set_xlabel('Processing')
ax.set_ylabel('tSNR')
handles, labels = b.get_legend_handles_labels()

labels = ['Dataset 2', 'Dataset 3']
plt.legend(handles=handles, labels=labels, frameon=False)
plt.subplots_adjust(bottom=0.15, left=0.2)

fig.savefig(output_path + os.sep + 'results' + os.sep + 'TSNR_both_datasets.svg', 
            transparent=True, bbox_inches='tight', format='svg')
plt.show()

# also collect the spatial smoothing values from 3dFWHMx
df_all = pd.DataFrame([])
for isub in subjects:
    sub_id = 'sub-sspr' + '{:02d}'.format(isub)
    func_path = output_path + os.sep + sub_id + os.sep + 'func'
    if (isub == 6) | (isub == 25) | (isub == 27) | (isub == 29):
        num_runs = 2
    elif isub == 31:
        num_runs = 3
    else:
        num_runs = 4
    for irun in range(num_runs):
        smooth_cord = float(pd.read_csv(func_path + os.sep + 'run-' + str(irun+1) + os.sep +
                                        'bold_acf_detrend_cord_stdout.txt', header=None)[0][1].split()[3])
        smooth_cord_denoised = float(pd.read_csv(func_path + os.sep + 'run-' + str(irun + 1) + os.sep +
                                                 'bold_denoised_acf_detrend_cord_stdout.txt',
                                                 header=None)[0][1].split()[3])
        df = pd.DataFrame([smooth_cord, smooth_cord_denoised], columns=['smoothing'])
        df['processing'] = ['raw', 'denoised']
        df['subject'] = sub_id
        df['run'] = irun+1
        if isub > 16:
            df['new'] = 1
        else:
            df['new'] = 0
        df_all = pd.concat([df_all, df])

# average over runs
avgs = df_all.groupby(['subject', 'processing']).mean()
avgs.reset_index(inplace=True)
avgs.drop('run', axis=1, inplace=True)

result = avgs.groupby(['processing', 'new']).mean(['smoothing'])
print(result)
denoised_old = result['smoothing']['denoised'][0.0]
denoised_new = result['smoothing']['denoised'][1.0]
raw_old = result['smoothing']['raw'][0.0]
raw_new = result['smoothing']['raw'][1.0]

perc_increase_old = ((denoised_old - raw_old) / raw_old) * 100
print("This equals an increase of {} in dataset 2".format(round(perc_increase_old, 3)))
result = ttest_rel(avgs[(avgs['new']==0) & (avgs['processing'] == 'denoised')]['smoothing'],
                   avgs[(avgs['new']==0) & (avgs['processing'] == 'raw')]['smoothing'])
print(result)
perc_increase_new = ((denoised_new - raw_new) / raw_new) * 100
print("This equals an increase of {} in dataset 3".format(round(perc_increase_new, 3)))
result = ttest_rel(avgs[(avgs['new']==1) & (avgs['processing'] == 'denoised')]['smoothing'],
                   avgs[(avgs['new']==1) & (avgs['processing'] == 'raw')]['smoothing'])
print(result)

