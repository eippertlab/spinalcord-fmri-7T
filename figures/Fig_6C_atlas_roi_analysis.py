import os
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp, ttest_rel, normaltest

main_path = '/data/pt_02661_raw/Heatpain/derivatives'
result_path = os.path.join(main_path, 'results')

color_ventral = [5/255, 113/255, 176/255]
color_superficial = [252/255, 141/255, 89/255]
color_middle = [227/255, 74/255, 51/255]
color_deep = [179/255, 0/255, 0/255]

sns.set_context("talk")
sns.set_style("ticks")

my_color = 'black'

design = 'block'  # block onset
smooth = '2mm_susan'  # 2mm_susan no_smooth
dataset = 'dataset2'

atlas_names = ['L7ToL10', 'L5L6', 'L3L4', 'L1L2']

df = pd.DataFrame([])
mask_dir = os.path.join(main_path, 'group_analysis', dataset, 'masks')

cope_file = os.path.join(main_path,
                         'group_analysis', dataset, design + '_design_' + 
                         smooth + '_vasa/cope1_merged_cut.nii.gz')
cope = nib.load(cope_file)
cope_img = cope.get_fdata()

for atlas in atlas_names:
    
    atlas_file = os.path.join(mask_dir, 'atlas_C6_' + atlas + '_left_PAM50_cut.nii.gz')
    mask = nib.load(atlas_file)
    mask_img = mask.get_fdata()
    
    avg_values = []
    for sub in range(cope_img.shape[3]):
        sub_cope = cope_img[:, :, :, sub]
        avg_sub = np.mean(sub_cope[mask_img > 0])
        avg_values.append(avg_sub)
    
    tmp = pd.DataFrame({atlas: avg_values})
    df = pd.concat([df, tmp], axis=1)

df.reset_index(inplace=True)
melted = df.melt(var_name='Atlas', value_name='avg_cope', id_vars='index')
df.to_csv(os.path.join(result_path, 'atlas_roi_single_sub_data_{}.csv'.format(design)), index=None)

for atlas in atlas_names:
    print('Test for {} > 0:'.format(atlas))
    data = melted[melted['Atlas']==atlas]['avg_cope']
    print(normaltest(data))
    print(ttest_1samp(data, 0, alternative='greater'))
    mean_diff = np.mean(data) - 0
    std_dev = np.std(data, ddof=1)
    cohen_d_value = mean_diff/std_dev
    print('Cohens D: {}'.format(cohen_d_value))


print('Test for L1L2 different from L5L6:')
print(ttest_rel(melted[melted['Atlas']=='L1L2']['avg_cope'], 
                melted[melted['Atlas']=='L5L6']['avg_cope']))
differences = melted[melted['Atlas']=='L1L2']['avg_cope'].values - melted[melted['Atlas']=='L5L6']['avg_cope'].values
mean_diff = np.mean(differences)
std_dev = np.std(differences, ddof=1)
cohen_d_value = mean_diff/std_dev
print('Cohens D: {}'.format(cohen_d_value))

color_order = [color_ventral, color_deep, color_middle, color_superficial]


fig, ax = plt.subplots(1, 1, figsize=(6, 5))
sns.boxplot(y='Atlas', x='avg_cope', data=melted, whis=[0, 100], width=0.6, hue='Atlas', palette=color_order)
sns.stripplot(y='Atlas', x='avg_cope', data=melted, size=5, color='.3', linewidth=0)
# plt.title('Individual mean cope in atlas ROI {}'.format(design))
plt.xlabel('Mean beta value')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(os.path.join(result_path, 'atlas_left_cope_avg_' + dataset + '_' +
                         design + '_' + smooth + '.png'),
            transparent=True, format='png', bbox_inches='tight')
plt.show()

