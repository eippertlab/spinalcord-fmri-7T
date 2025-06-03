import os
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

main_path = '/data/pt_02661_raw/Heatpain/derivatives'

sns.set_context("poster")
sns.set_style("ticks")

my_color = 'black'
plot_levels = True

design = 'block'
dataset = 'both'
thresh = 0.01


data_dir = os.path.join(main_path, 'group_analysis', dataset,
                        design + '_design_2mm_susan_vasa', 'PAM50_cord')

p_map_file = os.path.join(data_dir,'cope1_tfce_p_tstat1.nii.gz')
p_map = nib.load(p_map_file)
p_map_img = p_map.get_fdata()

img_data = np.zeros(p_map_img.shape)
img_data[p_map_img > (1-thresh)] = 1

# how was the original image cut in terms of template space coordinates (see group analysis commands)
if dataset == 'both':
    cut_low = 790
    num_slices = 65
elif dataset == 'dataset2':
    cut_low = 785
    num_slices = 70
elif dataset == 'dataset3':
    cut_low = 790
    num_slices = 69
else:
    raise ValueError('Please provide a valid dataset name')

if img_data.shape[2] != num_slices:
    raise ValueError('The provided number of slices does not match the image!')

active_voxel = []
for z in range(img_data.shape[2]):
    temp = img_data[:, :, z]
    active_voxel_per_slice = sum(sum(temp > 0))
    active_voxel.append(active_voxel_per_slice)

act_vox_df = pd.DataFrame({'num_voxel': active_voxel, 'slice': np.arange(1, img_data.shape[2] + 1)})
act_vox_df['slice_template'] = act_vox_df['slice'] + cut_low - 1

# apply gaussian filter for smoother line
act_vox_df['smoothed_count'] = np.transpose(gaussian_filter1d([act_vox_df['num_voxel']], sigma=1.5))

fig, ax = plt.subplots(1, 1, figsize=(10, 4))
plt.plot(act_vox_df['slice_template'], act_vox_df['smoothed_count'], color=my_color)

my_ylim = ax.get_ylim()
ax.set_ylim([-0.4, my_ylim[1]])

# show with a rectangle where each segment starts and ends
if plot_levels:
    y_min, y_max = ax.get_ylim()
    segment_df = pd.DataFrame({'label': ['C5', 'C6', 'C7'],
                               'start': [850, 819, 788], 'end': [882, 849, 818]})

    for i, row in segment_df.iterrows():
        print(row['label'])
        if row['start'] < cut_low + num_slices:
            # is the end still inside?
            if row['end'] < cut_low + num_slices:
                real_length = row['end'] - row['start'] + 1
                plot_it = True
            else:
                real_length = cut_low + num_slices - 1 - row['start']
                plot_it = True
            # is the start still inside? if not change length again by what is missing
            if row['start'] > cut_low:
                real_start = row['start']
                plot_it = True
            else:
                real_start = cut_low
                real_length = real_length - (cut_low - row['start'])
                plot_it = False
            print(real_start)
            print(real_length)
            if plot_it:
                plt.axvline(x=real_start, ls='--', color=my_color)

ax.set_xlim(cut_low + num_slices, cut_low - 1)  # decreasing
#plt.gca().invert_yaxis()

for label in ax.get_xticklabels(which='major'):
    label.set(rotation=90, horizontalalignment='center')
for label in ax.get_yticklabels(which='major'):
    label.set(rotation=90, horizontalalignment='center')
plt.xlabel('Slice')
plt.ylabel('Voxel count')
plt.subplots_adjust(bottom=0.15, left=0.1)

ax.spines['bottom'].set_color(my_color)
ax.spines['top'].set_color(my_color)
ax.spines['right'].set_color(my_color)
ax.spines['left'].set_color(my_color)
ax.xaxis.label.set_color(my_color)
ax.yaxis.label.set_color(my_color)
ax.tick_params(colors=my_color, which='both')  # 'both' refers to minor and major axes
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig(os.path.join(main_path, 'results',
                         'density_' + design + '_2mm_susan_vasa_' + 
                         dataset + '_uncorr_tmap_cord_only_levels_' +
                         str(plot_levels) + '.svg'), bbox_inches='tight', 
            transparent=True, format='svg')
plt.show()

