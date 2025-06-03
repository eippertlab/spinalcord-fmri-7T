import pandas as pd
import numpy as np
import os
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib import colormaps
import seaborn as sns

data_path = '/data/pt_02661_raw/Heatpain/derivatives'
result_path = '/data/pt_02661_raw/Heatpain/derivatives/results/even_odd'

dataset = 'both'

if dataset == 'both':
    subjects = np.arange(1, 42)
    # delete subjects that have been measured twice and the one that had no MEGRE scan
    subjects = np.delete(subjects, 32)
    subjects = np.delete(subjects, 11)
    subjects = np.delete(subjects, 0)
elif dataset == 'dataset1':
    subjects = np.arange(1, 17)
elif dataset == 'dataset2':
    subjects = np.arange(17, 42)
    subjects = np.delete(subjects, 16)

to_plot = 'avg'  # avg peak
stat = 'zstat'  # cope zstat

plot_single_sub = False
plot_correlations = True

even_color = [216/255, 179/255, 101/255]
odd_color = [90/255, 180/255, 172/255]
mix_color = [0.6, 0.7, 0.5]

sns.set_context("poster")
sns.set_style("ticks")

selected_slices_even = []
selected_slices_odd = []
correlations = []
voxelwise_correlations = []

for isub in subjects:
    sub_id = 'sub-sspr' + '{:02d}'.format(isub)
    print(sub_id)

    # even - statistical map
    stat_map_file = os.path.join(data_path, sub_id, 'func', 'second_level', 
                                 'block_design_2mm_susan_even_odd.gfeat',
                                 'cope1.feat', 'stats', stat + '1.nii.gz')
    stat_map = nib.load(stat_map_file)
    even_img = stat_map.get_fdata()
    
    # odd - statistical map
    stat_map_file = os.path.join(data_path, sub_id, 'func', 'second_level', 
                                 'block_design_2mm_susan_even_odd.gfeat',
                                 'cope2.feat', 'stats', stat + '1.nii.gz')
    stat_map = nib.load(stat_map_file)
    odd_img = stat_map.get_fdata()
    
    # individual dorsal horn mask
    mask_file = os.path.join(data_path, sub_id, 'anat_t2s', 
                             'MEGRE_RMS_MOCO_MEAN_gmseg_EPIspace_bspline_left_dorsal.nii.gz')
    mask = nib.load(mask_file)
    mask_img = mask.get_fdata()
    
    # extract average and peak activation of DH left within each slice
    even_peaks = []
    odd_peaks = []
    even_avgs = []
    odd_avgs = []
    num_slices = mask_img.shape[2]
    for slice_idx in range(num_slices):
        mask_slice = mask_img[:, :, slice_idx]
        even_slice = even_img[:, :, slice_idx]
        odd_slice = odd_img[:, :, slice_idx]
        
        even_values = even_slice[mask_slice > 0]
        odd_values = odd_slice[mask_slice > 0]
        
        if len(even_values) > 0:
            even_peak = np.max(even_values)
            odd_peak = np.max(odd_values)
            even_avg = np.mean(even_values)
            odd_avg = np.mean(odd_values)
        else:
            even_peak = np.nan
            odd_peak = np.nan
            even_avg = np.nan
            odd_avg = np.nan
        even_peaks.append(even_peak)
        odd_peaks.append(odd_peak)
        even_avgs.append(even_avg)
        odd_avgs.append(odd_avg)
    
    df_new = pd.DataFrame({'slice': np.arange(1, num_slices + 1),
                           'value_even_peak': even_peaks,
                           'value_odd_peak': odd_peaks,
                           'value_even_avg': even_avgs,
                           'value_odd_avg': odd_avgs})
    
    df_long = pd.wide_to_long(df_new, stubnames='value', i='slice', j='extracted', sep='_',
                              suffix='\w+').reset_index()
    df_long[['type', 'extracted']] = df_long['extracted'].str.split('_', expand=True)

    # correlate the two
    # if one of the values is NaN somewhere then exclude these
    df_even = df_long[(df_long['extracted'] == to_plot) & (df_long['type'] == 'even')]
    df_odd = df_long[(df_long['extracted'] == to_plot) & (df_long['type'] == 'odd')]
    df_even = df_even.dropna(axis=0)
    df_odd = df_odd.dropna(axis=0)
    correlation = np.corrcoef(df_even['value'], df_odd['value'])[0][1]
    correlations.append(correlation)
    
    if plot_single_sub:
        # check which slice is in which segment 
        # by taking a warped template segment mask
        # and asking where this mask has its maximum and minimum z value
        mask_dir = os.path.join(data_path, sub_id, 'masks_subjectspace')
        mask_indices = {}
        for segment in ['C5', 'C6', 'C7', 'C8']:
            mask = nib.load(os.path.join(mask_dir, segment + '.nii.gz'))
            mask_img = mask.get_fdata()
            indices = np.argwhere(mask_img == 1)
            if indices.size == 0:
                min_z = np.nan
                max_z = np.nan
            else:
                # here slices are coded from 0 to 14 instead of 1 to 15
                z_indices = indices[:, 2]
                min_z = np.min(z_indices) + 1
                max_z = np.max(z_indices) + 1
            mask_indices[segment] = {'min_z': min_z, 'max_z': max_z}
    
        # plot the activation values
        fig, ax = plt.subplots(1, 1, figsize=(4, 9))
        sns.lineplot(data=df_long[df_long['extracted'] == to_plot], 
                     y='slice', x='value', hue='type', orient='y', 
                     palette=[even_color, odd_color], hue_order=['even', 'odd'])
        ax.set_title("Subject {} correlation {}".format(sub_id, round(correlation, 3)))
    
        my_lims = ax.get_xlim()
        ax.set_xlim([my_lims[0]-0.5, my_lims[1]+0.5])
        
        ax.set_ylim([0, 17])
        
        # go from the bottom to the top
        segments = ['C8', 'C7', 'C6', 'C5']
        for s, segment in enumerate(segments):
            max_z = mask_indices[segment]['max_z']
            # if this segment is in there
            if not np.isnan(max_z):
                # and the next one is also still in there
                if s + 1 < len(segments):
                    # then determine where their border is
                    min_z_next = mask_indices[segments[s+1]]['min_z']
                    border = (max_z + min_z_next)/2
                    plt.axhline(y=border, color='black', linestyle='--')
        
        plt.legend()
        sns.despine()
        ax.get_legend().remove()
        #plt.tight_layout()
        plt.savefig(os.path.join(result_path, 
                                 'subject_space_even_odd_activation_{}_dhl_{}_{}.svg'.format(sub_id, stat, to_plot)),
                    transparent=True, format='svg', bbox_inches='tight')
        plt.show()
    
print(correlations)
print(np.mean(correlations))

if plot_correlations:
    corr_df = pd.DataFrame({'correlation': correlations})
    fig, ax = plt.subplots(1, 1, figsize=(12,4))
    ax.bar(corr_df.index + 1, corr_df['correlation'], color='darkgray')
    #ax.set_title('Correlation for each participant')
    #ax.set_xlabel('Subject index')
    ax.set_ylabel('Correlation')
    ax.set_ylim([-0.5, 1])
    #plt.tight_layout()
    sns.despine()
    plt.savefig(os.path.join(result_path, 
                             'subject_space_even_odd_correlation_dhl_{}_{}.svg'.format(stat, to_plot)),
                transparent=True, format='svg', bbox_inches='tight')
    plt.show()
    
    positions = [-0.4, 1]
    fig2, ax2 = plt.subplots(1, 1, figsize=(3, 4))
    v1 = ax2.violinplot(corr_df['correlation'], positions=[positions[0]], 
                        showmeans=False, showextrema=False,
                        showmedians=False, side="low", widths=0.8)
    for pc in v1['bodies']:
        pc.set_facecolor('gray')
        pc.set_edgecolor('black')
        pc.set_linewidth(1.2)
        pc.set_alpha(0.7)
    # `side` param requires matplotlib 3.9+
    a = sns.boxplot(ax=ax2, y=corr_df['correlation'], 
                    x=np.repeat(positions[1], corr_df.values.shape[0]),
                    showfliers=False, width=0.2, color='gray')
    for patch in a.patches:
        r, g, b, al = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.7))
    
    #ax2.set_title('Correlation over all participants')
    ax2.set_xlabel('')
    ax2.set_ylabel('Correlation')
    ax2.set_ylim([-0.5, 1])
    ax2.set_xlim([-1,1])
    #plt.tight_layout()
    sns.despine()
    plt.savefig(os.path.join(result_path, 
                             'subject_space_even_odd_correlation_avg_dhl_{}_{}.svg'.format(stat, to_plot)),
                transparent=True, format='svg', bbox_inches='tight')
    plt.show()
