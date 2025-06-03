import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

data_path = '/data/pt_02661_raw/Heatpain/derivatives'

subjects = np.arange(1, 42)
# delete subjects that have been measured twice
subjects = np.delete(subjects, 11)
subjects = np.delete(subjects, 0)

region = 'dh_left'

# C5 ends at z=882 and C8 starts at z=756 --> take only these
min_z = 755
max_z = 883
C6 = [819, 849]
C7 = [788, 818]

peaks = []
for isub in subjects:
    sub_id = 'sub-sspr' + '{:02d}'.format(isub)
    fname = os.path.join(data_path, sub_id, 'copes_templatespace',
                         'block_design_2mm_susan', 'extracted', region + '_cope1_all.txt')
    df = pd.read_csv(fname, header=None, sep='\s+').transpose()
    df = df.rename(columns={0: "x", 1: "y", 2: "z", 3: "value"})
    max_voxel = df.loc[df['value'].idxmax()]
    peaks.append(max_voxel.z)

sns.set_style('ticks')
sns.set_context("talk")
f, ax = plt.subplots(figsize=(3, 6))

x_value = 1

# Define the levels with start and end points
segment_df = pd.DataFrame({'label': ['C5', 'C6', 'C7', 'C8'],
                           'start': [850, 819, 788, 756],
                           'end': [882, 849, 818, 787]})

# plt.scatter([x_value] * len(peaks), peaks, color='black', alpha=0.3)
data = pd.DataFrame({'x': np.repeat(x_value, len(peaks)), 'y': peaks})
sns.regplot(data=data, marker='o', x='x', y='y',
            fit_reg=False, x_jitter=0.05, seed=np.random.seed(99),
            color='black', ax=ax,
            scatter_kws={'alpha': 0.3, 's': 64})
for i, segment in segment_df.iterrows():
    plt.axhline(y=segment['start']+0.5, color='gray', linewidth=1, linestyle='--')
    #plt.axhline(y=segment['end'], color='gray', linewidth=0.5, linestyle='--')
    plt.text(x_value + 0.5, (segment['start'] + segment['end']) / 2, segment['label'], 
             va='center', ha='right')

# also show the end of C5
plt.axhline(y=882+0.5, color='gray', linewidth=1, linestyle='--')

ax.set_xlim([0.5, 1.5])
ax.set_ylim([min_z, max_z])
# Labels and title
plt.xlabel('')
plt.ylabel('Inividual peak location z-direction')
# plt.title('Individual max beta peaks in left dorsal horn differ in their z position')
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
f.savefig(os.path.join(data_path, 'results', 'individual_max_peak_' + region + '_raw_beta.svg'),
          transparent=True, format='svg', bbox_inches='tight')
plt.show()

# how many are in each segment?
for i, segment in segment_df.iterrows():
    num_peaks = np.sum(np.logical_and(data['y'] >= segment['start'], 
                                      data['y'] < segment['end']))
    print('In segment {} there are {} peaks.'.format(segment['label'], num_peaks))

