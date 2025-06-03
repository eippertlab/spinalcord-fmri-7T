# Ulrike Horn
# uhorn@cbs.mpg.de

import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

output_path = '/data/pt_02661_raw/SatPad_Pilot_Merve/data/'
result_path = '/data/pt_02661_raw/Heatpain/derivatives/results/SatPad/'

subjects = [5]

color_green = [27/255, 158/255, 119/255]
color_orange = [217/255, 95/255, 2/255]

for isub in subjects:
    sub = 'sub-' + str("{:02d}".format(isub))
    print('Starting with subject: ' + sub)

    # Load data
    on_values = np.loadtxt(os.path.join(output_path, sub, 'SAT-ON',
                                        'output', 'active_FM_AS_T1_CroppedLabels.txt'))[3, :]
    off_values = np.loadtxt(os.path.join(output_path, sub, 'SAT-OFF',
                                         'output', 'active_FM_AS_T1_CroppedLabels.txt'))[3, :]

    # Set up the figure
    sns.set_context('poster')
    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=(7, 4))
    # Plot histograms
    ax.hist(on_values, bins=np.arange(min(on_values), max(on_values) + 0.01, 0.01),
            color=color_green, edgecolor='black', linewidth=1, alpha=0.5, label='ON')
    ax.hist(off_values, bins=np.arange(min(off_values), max(off_values) + 0.01, 0.01),
            color=color_orange, edgecolor='black', linewidth=1, alpha=0.5, label='OFF')
    # Set limits and labels
    ax.set_xlim([-0.1, 0.2])
    ax.set_ylim([0, 300])
    ax.set_xlabel('Ppm')
    ax.set_ylabel('Count')
    sns.despine()
    plt.tight_layout()
    fig.savefig(os.path.join(result_path, 'histogram_PPM_active_sub-05.svg'),
                bbox_inches='tight', transparent=True, format='svg')
    plt.show()
