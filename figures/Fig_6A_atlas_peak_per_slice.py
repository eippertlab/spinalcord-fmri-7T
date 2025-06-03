import pandas as pd
import numpy as np
import os
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


data_path = '/data/pt_02661_raw/Heatpain/derivatives'
result_path = '/data/pt_02661_raw/Heatpain/derivatives/results/'

design = 'block'
smooth = '2mm_susan'  # 2mm_susan no_smooth
dataset = 'dataset1'  # both

t_map_file = os.path.join(data_path,
                          'group_analysis', dataset, design + '_design_' + smooth +
                          '_vasa/PAM50_cord/cope1_tstat1.nii.gz')

p_map_file = os.path.join(data_path,
                          'group_analysis', dataset, design + '_design_' + smooth +
                          '_vasa/PAM50_cord/cope1_tfce_p_tstat1.nii.gz')

mask_dir = os.path.join(data_path, 'group_analysis', dataset, 'masks')

rois = ['L1L2_left', 'L1L2_right',
        'L3L4_left', 'L3L4_right',
        'L5L6_left', 'L5L6_right',
        'L7ToL10_left', 'L7ToL10_right']

color_ventral = [5/255, 113/255, 176/255]
color_superficial = [252/255, 141/255, 89/255]
color_middle = [227/255, 74/255, 51/255]
color_deep = [179/255, 0/255, 0/255]

color_order = [color_superficial, color_middle, color_deep, color_ventral]

# threshold the data like you did for the activation images
t_map = nib.load(t_map_file)
t_map_img = t_map.get_fdata()

p_map = nib.load(p_map_file)
p_map_img = p_map.get_fdata()

# set all t-values where p>0.01 = 1-0.99 uncorrected to -99, so they don't get selected
t_map_img[p_map_img < 0.99] = -99

# open each atlas ROI and determine the maxima per slice
# then compare for each slice to determine which ROI has the highest maximum = maximum overall
slice_maxima = {roi: [] for roi in rois}
for roi in rois:
    atlas = os.path.join(mask_dir, 'atlas_C6_' + roi + '_PAM50_cut.nii.gz')
    mask = nib.load(atlas)
    mask_img = mask.get_fdata()

    for slice_idx in range(t_map_img.shape[2]):
        t_map_slice = t_map_img[:, :, slice_idx]
        mask_slice = mask_img[:, :, slice_idx]

        masked_values = t_map_slice[mask_slice > 0]
        if masked_values.size > 0:
            peak_value = np.nanmax(masked_values)
        else:
            peak_value = np.nan  # for empty slices
        slice_maxima[roi].append(peak_value)

# determine who 'won'
mask_wins = {roi: 0 for roi in rois}
for slice_idx in range(t_map_img.shape[2]):
    slice_peaks = {roi: maxima[slice_idx] for roi, maxima in slice_maxima.items()}
    valid_peaks = {roi: peak for roi, peak in slice_peaks.items() if not np.isnan(peak)}
    if valid_peaks:
        winning_mask = max(slice_peaks, key=slice_peaks.get)
        mask_wins[winning_mask] += 1

for roi, win_count in mask_wins.items():
    print(f"Mask {roi} won {win_count} times.")

df_wins = pd.DataFrame(list(mask_wins.items()), columns=['Mask', 'Win_count'])
print(df_wins)

df_wins[['Mask', 'Side']] = df_wins['Mask'].str.extract(r'(.+)_(left|right)')

sns.set_theme(style='ticks')
sns.set_context('talk')
fig, axs = plt.subplots(1, 2, figsize=(10, 8), sharey=True)
sns.barplot(data=df_wins[df_wins['Side'] == 'right'], x='Mask', y='Win_count', hue='Mask',
            ax=axs[0], palette=color_order)
axs[0].set_title('Right side')
axs[0].set_ylabel('# gray matter maxima per slice')
sns.barplot(data=df_wins[df_wins['Side'] == 'left'], x='Mask', y='Win_count', hue='Mask',
            ax=axs[1], palette=color_order)
axs[1].set_title('Left side')
axs[1].set_ylabel('# gray matter slice maxima')
fig.suptitle('Number of t-value peaks in each atlas - {} design'.format(design))
plt.tight_layout()
fig.savefig(os.path.join(result_path, 'Atlas_Peak_value_counts_{}_{}.png'.format(design, smooth)), transparent=True)
plt.show()

fig2, axs2 = plt.subplots(1, 1, figsize=(3, 5), sharey=True)
g = sns.barplot(data=df_wins[df_wins['Side'] == 'left'], y='Mask', x='Win_count', hue='Mask',
                ax=axs2, order=['L7ToL10', 'L5L6', 'L3L4', 'L1L2'], palette=color_order)
axs2.set_title('Left side')
axs2.set_xlabel('# gray matter slice maxima')
axs2.set_ylabel('')
fig2.suptitle('{} design'.format(design))
plt.tight_layout()
axs2.spines['top'].set_visible(False)
axs2.spines['right'].set_visible(False)
fig2.savefig(os.path.join(result_path, 
                          'Atlas_Peak_value_counts_{}_swapped_{}.svg'.format(design, smooth)), 
             transparent=True, format='svg', bbox_inches='tight')
plt.show()
