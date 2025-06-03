import os
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt


sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
sns.set_context('talk')

main_path = '/data/pt_02661_raw/Heatpain/derivatives'

color_ventral = [5/255, 113/255, 176/255]
color_superficial = [252/255, 141/255, 89/255]
color_middle = [227/255, 74/255, 51/255]
color_deep = [179/255, 0/255, 0/255]

my_color = 'black'

smooth = '2mm_susan'  # 2mm_susan no_smooth

thresh = 0.05
atlas_names = ['L7ToL10_left', 'L5L6_left', 'L3L4_left', 'L1L2_left']

for design in ['block', 'onset']:
    # for dataset in ['both']:
    for dataset in ['dataset1', 'dataset2']:

        data = []
        mask_path = os.path.join(main_path, 'group_analysis', dataset, 'masks')

        data_path = os.path.join(main_path, 'group_analysis', dataset,
                                 design + '_design_' + smooth + '_vasa')
        # read t map and then mask it with different masks
        t_map = nib.load(os.path.join(data_path, 'PAM50_cord', 'cope1_tstat1.nii.gz'))
        t_map_img = t_map.get_fdata()
        for atlas in atlas_names:
            atlas_file = os.path.join(mask_path, 'atlas_C6_' + atlas + '_PAM50_cut.nii.gz')
            mask = nib.load(atlas_file)
            mask_img = mask.get_fdata()

            masked_values = t_map_img[mask_img == 1]

            tmp = pd.DataFrame({'atlas': atlas, 'tval': masked_values})

            data.append(tmp)

        # Concatenate all dataframes at once
        data = pd.concat(data, ignore_index=True)
        # data_thresh = data[data['tval'] > 0]

        color_order = [color_ventral, color_deep, color_middle, color_superficial]

        # Initialize the FacetGrid object
        g = sns.FacetGrid(data, row="atlas", hue="atlas", aspect=4, height=1.2, palette=color_order)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, "tval",
              bw_adjust=0.5, clip_on=False,
              fill=True, alpha=1, linewidth=1.5)
        g.map(sns.kdeplot, "tval", clip_on=False, color="black", lw=1, bw_adjust=.5)

        # draw a vertical line at 0
        g.map(plt.axvline, x=0, linewidth=1, linestyle='--', color=my_color)

        # passing color=None to refline() uses the hue mapping
        g.refline(y=0, linewidth=0.2, linestyle="-", color=None, clip_on=False)

        # Define and use a simple function to label the plot in axes coordinates
        def label(x, color, label):
            ax = plt.gca()
            # ax.text(0, .2, label, fontweight="bold", color=color,
            #         ha="left", va="center", transform=ax.transAxes)
            ax.xaxis.label.set_color(my_color)
            ax.yaxis.label.set_color(my_color)
            ax.tick_params(colors=my_color, which='both')  # 'both' refers to minor and major axes

        g.map(label, "tval")

        # Set the subplots to overlap
        g.figure.subplots_adjust(hspace=-.25)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[], ylabel="", xlabel='T-values')

        g.despine(bottom=True, left=True)

        g.set(xlim=(-2.5, 5.5))

        g.savefig(os.path.join(main_path,
                               'results', 'atlas_thresh_tvalues_{}_{}_{}_{}.svg'.format(design, dataset,
                                                                                        my_color, smooth)),
                  transparent=True, format='svg', bbox_inches='tight')
        plt.show()
