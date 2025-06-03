#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ulrike Horn
29.11.2024
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import nibabel as nib

project_dir = '/data/pt_02661_raw/Heatpain/derivatives/group_analysis'
out_dir = '/data/pt_02661_raw/Heatpain/derivatives/results/spatial_specificity'

dataset = 'both'
design = 'block_design_2mm_susan_vasa'
with_legend = True
thresh = 0.05

data_dir = os.path.join(project_dir, dataset, design, 'PAM50_cord')
mask_dir = os.path.join(project_dir, dataset, 'masks')

vr_color = [130/255, 193/255, 240/255]  # bright blue
vl_color = [42/255, 77/255, 105/255]  # dark blue
dr_color = [244/255, 125/255, 125/255]  # bright red
dl_color = [178/255, 34/255, 34/255]  # dark red
r_color = [187/255, 159/255, 183/255] # mix of bright colors
l_color = [110/255, 56/255, 70/255]  # mix of dark colors
v_color = [86/255, 135/255, 173/255] # mix of blues
d_color = [211/255, 80/255, 80/255] # mix of reds

sns.set_context("poster")
sns.set_style("ticks")

my_color = 'black'

rois_horns = ['dh_left', 'dh_right', 'vh_left', 'vh_right']
levels = ['C5', 'C6', 'C7']

p_map_file = os.path.join(data_dir, 'cope1_tfce_p_tstat1.nii.gz')

# threshold the data like you did for the activation images
p_map = nib.load(p_map_file)
p_map_img = p_map.get_fdata()

counts = []
rois = []
curr_levels = []
mask_sizes = []
for roi in rois_horns:
    for level in levels:
        mask = nib.load(os.path.join(mask_dir, roi + '_' + level + '_cut.nii.gz'))
        mask_img = mask.get_fdata()
        mask_sizes.append(np.sum(mask_img>0))
        
        masked_p_map = p_map_img[mask_img > 0]
        values = masked_p_map.flatten()
        p_values = 1 - values
        
        counts.append(sum(p_values < thresh))
        rois.append(roi)
        curr_levels.append(level)

df = pd.DataFrame({'val': counts,'roi': rois, 'level': curr_levels, 
                   'mask_size': mask_sizes})
df['perc'] = df['val'] / df['mask_size'] * 100

custom_palette = {
    'dh_left': dl_color,
    'dh_right': dr_color,
    'vh_left': vl_color,
    'vh_right': vr_color
}
mpl.rcParams['text.color'] = my_color
mpl.rcParams['xtick.color'] = my_color
mpl.rcParams['ytick.color'] = my_color
mpl.rcParams['axes.labelcolor'] = my_color
mpl.rcParams['axes.edgecolor'] = my_color

if dataset == 'both':
    fig, ax = plt.subplots(1, 1, figsize=(9, 8))
else:
    fig, ax = plt.subplots(1, 1, figsize=(9, 5))

sns.barplot(data=df, x='level', y='perc', hue='roi', palette=custom_palette, hue_order=rois_horns)
plt.title('Active voxel in each horn')
plt.xlabel('Level')
plt.ylabel('% active voxel')

if 'onset' in design:
    plt.ylim(0, 50)
elif 'block' in design:
    plt.ylim((0, 100))
    
if with_legend:
    sns.move_legend(ax, "upper left")
    legend = ax.legend()
    legend.set_title('')
    proper_labels = ['left dorsal', 'right dorsal', 'left ventral', 'right ventral']
    for t, l in zip(legend.texts, proper_labels):
        t.set_text(l)
    legend.get_frame().set_linewidth(0)
    legend.get_frame().set_alpha(0)
else:
    ax.get_legend().remove()
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(os.path.join(out_dir, 'barchart_horns_' + design + 
                         '_2mm_susan_vasa_' + dataset + '_uncorr_tmap_percent.svg'),
            transparent=True, format='svg', bbox_inches='tight')

plt.show()

df[['d_v', 'left_right']] = df['roi'].str.extract(r'([dv]h)_(left|right)')

# left vs right
fig2, ax2 = plt.subplots(1, 1, figsize=(9, 8))
df_LR = df.groupby(['left_right', 'level'])[['val', 'mask_size']].sum()
df_LR['perc'] = df_LR['val'] / df_LR['mask_size'] * 100
custom_palette = {
    'left': l_color,
    'right': r_color
    }
sns.barplot(data=df_LR, x='level', y='perc', hue='left_right', 
            palette=custom_palette, hue_order=['left', 'right'])
plt.title('Active voxel left vs right')
plt.xlabel('Level')
plt.ylabel('% active voxel')
if 'onset' in design:
    plt.ylim(0, 50)
elif 'block' in design:
    plt.ylim((0, 100))
plt.tight_layout()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig2.savefig(os.path.join(out_dir, 'barchart_horns_' + design + 
                         '_2mm_susan_vasa_' + dataset + '_uncorr_tmap_percent_LR.svg'),
            transparent=True, format='svg', bbox_inches='tight')
plt.show()

# dorsal vs ventral
fig3, ax3 = plt.subplots(1, 1, figsize=(9, 8))
df_DV = df.groupby(['d_v', 'level'])[['val', 'mask_size']].sum()
df_DV['perc'] = df_DV['val'] / df_DV['mask_size'] * 100
custom_palette = {
    'dh': d_color,
    'vh': v_color
    }
sns.barplot(data=df_DV, x='level', y='perc', hue='d_v', 
            palette=custom_palette, hue_order=['dh', 'vh'])
plt.title('Active voxel dorsal vs ventral')
plt.xlabel('Level')
plt.ylabel('% active voxel')
if 'onset' in design:
    plt.ylim(0, 50)
elif 'block' in design:
    plt.ylim((0, 100))
plt.tight_layout()
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
fig3.savefig(os.path.join(out_dir, 'barchart_horns_' + design + 
                         '_2mm_susan_vasa_' + dataset + '_uncorr_tmap_percent_DV.svg'),
            transparent=True, format='svg', bbox_inches='tight')
plt.show()
