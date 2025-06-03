import pandas as pd
import numpy as np
import nibabel as nib
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import t

data_path = '/data/pt_02661_raw/Heatpain/derivatives/group_analysis'
result_path = '/data/pt_02661_raw/Heatpain/derivatives/results/vasa'

if not os.path.exists(result_path):
    os.makedirs(result_path)

sns.set_context('poster')
sns.set_style('ticks')

design = 'block'
thresh = 0.01
mask = 'dh_left_C6'
vasa_path = os.path.join(data_path, 'both', design + '_design_2mm_susan_vasa')
normal_path = os.path.join(data_path, 'both', design + '_design_2mm_susan')

# vasa t and p stat
vasa_t_file = nib.load(os.path.join(vasa_path, mask, 'cope1_tstat1.nii.gz'))
vasa_t = vasa_t_file.get_fdata()

vasa_p_file = nib.load(os.path.join(vasa_path, mask, 'cope1_tfce_p_tstat1.nii.gz'))
vasa_p = vasa_p_file.get_fdata()

# no vasa t and p stat
t_file = nib.load(os.path.join(normal_path, mask, 'cope1_tstat1.nii.gz'))
no_vasa_t = t_file.get_fdata()

p_file = nib.load(os.path.join(normal_path, mask, 'cope1_tfce_p_tstat1.nii.gz'))
p = p_file.get_fdata()

# select only the voxels that are both significant
p_value_mask = np.logical_and(vasa_p > 1-thresh, p > 1-thresh)

vasa_tvals = vasa_t[p_value_mask].flatten()
tvals = no_vasa_t[p_value_mask].flatten()

df_thresh = pd.DataFrame({'without': tvals, 'vasa': vasa_tvals})
long_df = pd.melt(df_thresh, value_vars=['without', 'vasa'], var_name='type', value_name='values')

print('The t stats changed from without vasa to vasa by {}'
      ' on average'.format(round(np.mean(df_thresh['vasa'] - df_thresh['without']), 2)))

df_thresh['perc_incr'] = ((df_thresh['vasa'] - df_thresh['without']) / df_thresh['without']) * 100
print('This equals {}'
      ' %.'.format(round(np.mean(df_thresh['perc_incr']), 2)))
print('The range of improvement is from {} to {} %.'.format(round(np.min(df_thresh['perc_incr']), 2),
                                                            round(np.max(df_thresh['perc_incr']), 2)))

fig, axs = plt.subplots(1, 1, figsize=(6, 7))
# fig.suptitle("VasA", fontsize=30)
g = sns.scatterplot(x="without", y="vasa", data=df_thresh, 
                    alpha=0.5, c='gray', 
                    edgecolors='black', s=85)
axs.set_aspect('equal', 'box')
axs.margins(0.1)
xlim = axs.get_xlim()
ylim = axs.get_ylim()
new_min = np.round(min(xlim[0], ylim[0]), 1)
new_max = np.ceil(max(xlim[1], ylim[1]))
plt.xlim(new_min, new_max)
plt.ylim(new_min, new_max)

axs.spines['top'].set_visible(False)
axs.spines['right'].set_visible(False)
axs.set_xlabel('T-values wihout VasA', fontsize=25)
axs.set_ylabel('T-values with VasA', fontsize=25)
diag_line, = axs.plot(axs.get_xlim(), axs.get_ylim(), ls="--", c=".3")
# plt.title('T-values {} design {}'.format(design, mask))
plt.tight_layout()
plt.savefig(os.path.join(result_path, 
                         'vasa_tvalues_scatter_{}_{}.svg'.format(design, mask)), 
            transparent=True, bbox_inches='tight', format='svg')
plt.show()
