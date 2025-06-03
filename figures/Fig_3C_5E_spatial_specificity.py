#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code by Alice Dabbagh,
changed by Ulrike Horn
29.10.2024
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import nibabel as nib

project_dir = '/data/pt_02661_raw/Heatpain/derivatives/group_analysis'
out_dir = '/data/pt_02661_raw/Heatpain/derivatives/results/spatial_specificity'

dataset = 'both'
design = 'block_design_2mm_susan_vasa'
stat = 'tfce'
thresh = 0.05

data_dir = os.path.join(project_dir, dataset, design, 'PAM50_cord')
mask_dir = os.path.join(project_dir, dataset, 'masks')

t_map_file = os.path.join(data_dir,'cope1_tstat1.nii.gz')
t_map = nib.load(t_map_file)
t_map_img = t_map.get_fdata()

p_map_file = os.path.join(data_dir,'cope1_tfce_p_tstat1.nii.gz')
p_map = nib.load(p_map_file)
p_map_img = p_map.get_fdata()

cord_file = os.path.join(data_dir, 'PAM50_cord_cut.nii.gz')
cord = nib.load(cord_file)
cord_img = cord.get_fdata()

gm_file = os.path.join(mask_dir, 'PAM50_gm_bin_cut.nii.gz')
gm = nib.load(gm_file)
gm_img = gm.get_fdata()

# overall values:
gm_t = t_map_img[gm_img == 1]
print('The avg t-value in gray matter is {}'.format(np.mean(gm_t)))
num_vox_gm = np.sum(gm_img == 1)
wm_t = t_map_img[(gm_img == 0) & (cord_img == 1)]
print('The avg t-value in white matter is {}'.format(np.mean(wm_t)))
num_vox_wm = np.sum((gm_img == 0) & (cord_img == 1))

# only significantly activated:
gm_sign_t_mean = np.mean(t_map_img[(p_map_img > (1-thresh)) & (gm_img == 1)])
print('The avg t-value of the significant voxels in gray matter is {}'.format(gm_sign_t_mean))
num_sign_vox_gm = np.sum((p_map_img > (1-thresh)) & (gm_img == 1))
wm_sign_t_mean = np.mean(t_map_img[(p_map_img > (1-thresh)) & (gm_img == 0) & (cord_img == 1)])
print('The avg t-value of the significant voxels in white matter is {}'.format(wm_sign_t_mean))
num_sign_vox_wm = np.sum((p_map_img > (1-thresh)) & (gm_img == 0) & (cord_img == 1))

print('{} % of voxels in gray matter is significantly activated'.format(round(num_sign_vox_gm/num_vox_gm*100, 2)))
print('{} % of voxels in white matter is significantly activated'.format(round(num_sign_vox_wm/num_vox_wm*100, 2)))

x, y, z = np.indices(p_map_img.shape)
x_flat = x.flatten()
y_flat = y.flatten()
z_flat = z.flatten()
values_flat = p_map_img.flatten()
t_values_flat = t_map_img.flatten()

df = pd.DataFrame({'x': x_flat, 'y': y_flat, 'z': z_flat, 
                   'val': values_flat, 'tval': t_values_flat})
df['roi'] = 'cord'

x_gm = x[gm_img == 1]
y_gm = y[gm_img == 1]
z_gm = z[gm_img == 1]
values_gm = p_map_img[gm_img == 1]
t_values_gm = t_map_img[gm_img == 1]

df_gm = pd.DataFrame({'x': x_gm, 'y': y_gm, 'z': z_gm, 
                      'val': values_gm, 'tval': t_values_gm})
df_gm['roi'] = 'gm'

x_wm = x[(gm_img == 0) & (cord_img == 1)]
y_wm = y[(gm_img == 0) & (cord_img == 1)]
z_wm = z[(gm_img == 0) & (cord_img == 1)]
values_wm = p_map_img[(gm_img == 0) & (cord_img == 1)]
t_values_wm = t_map_img[(gm_img == 0) & (cord_img == 1)]

df_wm = pd.DataFrame({'x': x_wm, 'y': y_wm, 'z': z_wm, 
                      'val': values_wm, 'tval': t_values_wm})
df_wm['roi'] = 'wm'

data = pd.concat([df, df_gm, df_wm])

data['pval'] = 1 - data["val"]
data_thresh = data[data['pval'] < thresh]
data_thresh = data_thresh.drop(["z"], axis=1)
x = data_thresh["x"]
y = data_thresh["y"]

# examplary slice for outline
example_slice = 52

slice_cord = cord_img[:, :, example_slice]
slice_cord_array = pd.DataFrame(columns=["x", "y"])
for x in range(slice_cord.shape[0]):
    for y in range(slice_cord.shape[1]):
        if slice_cord[x, y] == 1:
            slice_cord_array = pd.concat([slice_cord_array, 
                                          pd.DataFrame({"x": float(x), "y": float(y)}, index=[0])],
                                         ignore_index=True)
slice_cord_array["x"] = pd.to_numeric(slice_cord_array['x'], downcast='float')
slice_cord_array["y"] = pd.to_numeric(slice_cord_array['y'], downcast='float')

# prep axes
whole = data_thresh[data_thresh['roi'] == 'cord']
gm = data_thresh[data_thresh['roi'] == 'gm']
wm = data_thresh[data_thresh['roi'] == 'wm']

x1 = gm["x"].values
# x2 = whole["x"].values
x2 = wm["x"].values

y1 = gm["y"].values
# y2 = whole["y"].values
y2 = wm["y"].values


# Spatial specificity of BOLD responses.
mpl.rcParams['pdf.fonttype'] = 42
sns.set(style="white")  

color_gm = "red"

with sns.plotting_context('paper', font_scale=2):
    fig, axScatter = plt.subplots(figsize=(12, 7))
        
    sns.regplot(data=wm, marker='o', x='x', y='y',
                fit_reg=False, x_jitter=0.5, seed=np.random.seed(2),
                y_jitter=0.5, color="grey", ax=axScatter,
                scatter_kws={'alpha': 0.3, 's': 8}, label="White matter")
    sns.regplot(data=gm, marker='o', x='x', y='y',
                fit_reg=False, x_jitter=0.5, seed=np.random.seed(22),
                y_jitter=0.5, color=color_gm, ax=axScatter,
                scatter_kws={'alpha': 0.3, 's': 8}, label="Gray matter")
  
    axScatter.set_xlim([55, 85])
    axScatter.set_ylim([62, 83])
    sns.despine(left=True, bottom=True)
    
    divider = make_axes_locatable(axScatter)
    axHistx = divider.append_axes("bottom", 1.2, pad=0.2, sharex=axScatter)
    axHisty = divider.append_axes("left", 1.2, pad=0.2, sharey=axScatter)
    
    # make some labels invisible
    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)

    # add density plots, x
    sns.kdeplot(x=np.array(x2), bw_method=0.5, ax=axHistx, color="grey", linewidth=0.6)
    sns.kdeplot(x=np.array(x1), bw_method=0.5, ax=axHistx, color=color_gm, linewidth=0.6)
     
    sns.kdeplot(y=np.array(y2), bw_method=0.5, ax=axHisty, color="grey", linewidth=0.6)
    sns.kdeplot(y=np.array(y1), bw_method=0.5, ax=axHisty, color=color_gm, linewidth=0.6)
    
    # add contour
    sns.kdeplot(data=slice_cord_array, x="x", y="y", levels=1, bw_method=0.2, 
                color="k", linewidths=1, ax=axScatter)
    axHisty.invert_xaxis()
    axHistx.invert_yaxis()

    # despine axes
    axHistx.set_axis_off()
    axHisty.set_axis_off()
    axScatter.legend(bbox_to_anchor=(1.3, 1), loc='upper left', 
                     borderaxespad=0, markerscale=5)
    
    # set the labels
    axScatter.set_xticks([59, 81])
    axScatter.set_xticklabels([' ', '  '])
    axScatter.set(xlabel=None)
    axScatter.set_yticks([64, 79])
    axScatter.set_yticklabels([' ', ' '])
    axScatter.set(ylabel=None)
    
    mpl.axes.Axes.set_aspect(axScatter, aspect="equal")
    # save
    plt.savefig(os.path.join(out_dir, 
                             'spatial_spec_{}_{}_{}_{}_gm.svg'.format(dataset, 
                                                                      design, 
                                                                      stat, 
                                                                      thresh)), 
                bbox_inches='tight', format="svg", dpi=300, transparent=True)
    plt.show()
