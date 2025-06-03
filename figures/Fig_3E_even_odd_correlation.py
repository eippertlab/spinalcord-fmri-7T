import pandas as pd
import numpy as np
import os
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from scipy import stats

project_dir = '/data/pt_02661_raw/Heatpain/derivatives/group_analysis'
out_dir = '/data/pt_02661_raw/Heatpain/derivatives/results/even_odd'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dataset = 'both'
design = 'block_design_2mm_susan_vasa_even_odd'
my_color = 'black'

vr_color = [130/255, 193/255, 240/255]  # bright blue
vl_color = [42/255, 77/255, 105/255]  # dark blue
dr_color = [244/255, 125/255, 125/255]  # bright red
dl_color = [178/255, 34/255, 34/255]  # dark red


def get_dataframe(mask):
    group_path = (os.path.join(project_dir, dataset, design, 'PAM50_cord'))
    mask_path = os.path.join(project_dir, dataset, 'masks')
    even_path = os.path.join(group_path, 'cope1_tstat1.nii.gz')
    odd_path = os.path.join(group_path, 'cope2_tstat1.nii.gz')
    even = nib.load(even_path)
    even_data = even.get_fdata()
    odd = nib.load(odd_path)
    odd_data = odd.get_fdata()

    my_mask = os.path.join(mask_path, mask + '_cut.nii.gz')
    mask = nib.load(my_mask)
    mask_data = mask.get_fdata()

    # use the mask to mask the cope image
    even_masked = even_data.copy()
    even_masked[mask_data == 0] = np.nan
    odd_masked = odd_data.copy()
    odd_masked[mask_data == 0] = np.nan

    # flatten data
    even_flat = even_masked.flatten()
    odd_flat = odd_masked.flatten()

    # does it take out the same voxel?
    even = even_flat[~np.isnan(even_flat)]
    odd = odd_flat[~np.isnan(odd_flat)]

    result_df = pd.DataFrame({'even': even, 'odd': odd})
    return result_df


sns.set_context("poster")
sns.set_theme(style="white")

mpl.rcParams['text.color'] = my_color
mpl.rcParams['xtick.color'] = my_color
mpl.rcParams['ytick.color'] = my_color
mpl.rcParams['axes.labelcolor'] = my_color

my_limits = [-3.5, 6.5]
fig, axs = plt.subplots(4, 4, figsize=(10, 8),
                        gridspec_kw={'hspace': 0.1,
                                     'wspace': 0.1,
                                     'width_ratios': [1, 5, 5, 1],
                                     'height_ratios': [1, 5, 5, 1]})

# 1. Ventral horn right
df = get_dataframe('vh_right')
r, p = stats.pearsonr(df['even'], df['odd'])
# axs[0, 1].title.set_text('Ventral horn right')
# Upper part charts
sns.histplot(data=df, x="even", kde=True, bins=10, ax=axs[0, 1], color=vr_color)
axs[0, 1].set_xlim(my_limits)

axs[0, 0].axis("off")
axs[0, 1].axis("off")
axs[1, 1].axis("off")
axs[1, 2].axis("off")
axs[0, 2].axis("off")
axs[0, 3].axis("off")
axs[1, 3].axis("off")
# Left part charts
sns.histplot(data=df, y='odd', kde=True, bins=10, ax=axs[1, 0], color=vr_color)
axs[1, 0].invert_xaxis()
for side in ['top', 'right', 'bottom', 'left']:
    axs[1, 0].spines[side].set_visible(False)
# axs[1, 0].tick_params(axis='x', which='both', labelbottom=False, bottom=False)
axs[1, 0].get_xaxis().set_visible(False)
axs[1, 0].set_ylim(my_limits)

# Linear regression for the middle plot
axs[1, 1].axhline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
axs[1, 1].axvline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
# scatter middle part
# sns.scatterplot(x='even', y='odd', data=df, ax=axs[1, 1], edgecolor=None, alpha=0.2, zorder=2)
# sns.regplot(x="even", y="odd", data=df, color='black', marker='+', ax=axs[1, 1], scatter=True)
sns.regplot(data=df, marker='o', x='even', y='odd',
            fit_reg=True, ax=axs[1, 1], ci=0,
            scatter_kws={'alpha': 0.3, 's': 8, 'color': vr_color}, line_kws={'alpha': 1, 'color': my_color})
axs[1, 1].set_xlim(my_limits)
axs[1, 1].set_ylim(my_limits)

rect = patches.Rectangle(
    (axs[1, 1].get_position().x0, axs[1, 1].get_position().y0),
    axs[1, 1].get_position().width,
    axs[1, 1].get_position().height,
    linewidth=1.5, edgecolor=my_color, facecolor='none', transform=fig.transFigure
)
fig.patches.append(rect)

# correlation text
axs[1, 1].text(5, -2, 'r = ' + str(np.round(r, 3)), horizontalalignment='right', color=my_color)
axs[1, 1].text(-1, 5, 'Ventral horn right')


# 2. Ventral horn left
df = get_dataframe('vh_left')
r, p = stats.pearsonr(df['even'], df['odd'])
# axs[0, 2].title.set_text('Ventral horn left')
# Upper part charts
sns.histplot(data=df, x="even", kde=True, bins=10, ax=axs[0, 2], color=vl_color)
axs[0, 2].set_xlim(my_limits)
# Right part charts
sns.histplot(data=df, y='odd', kde=True, bins=10, ax=axs[1, 3], color=vl_color)
axs[1, 3].set_ylim(my_limits)

# Linear regression for the middle plot
axs[1, 2].axhline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
axs[1, 2].axvline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
sns.regplot(data=df, marker='o', x='even', y='odd',
            fit_reg=True, ax=axs[1, 2], ci=0,
            scatter_kws={'alpha': 0.3, 's': 8, 'color': vl_color}, line_kws={'alpha': 1, 'color': my_color})
axs[1, 2].set_xlim(my_limits)
axs[1, 2].set_ylim(my_limits)

rect = patches.Rectangle(
    (axs[1, 2].get_position().x0, axs[1, 2].get_position().y0),
    axs[1, 2].get_position().width,
    axs[1, 2].get_position().height,
    linewidth=1.5, edgecolor=my_color, facecolor='none', transform=fig.transFigure
)
fig.patches.append(rect)
axs[1, 2].text(5, -2, 'r = ' + str(np.round(r, 3)), horizontalalignment='right', color=my_color)
axs[1, 2].text(-1, 5, 'Ventral horn left')

# 3. Dorsal horn right
df = get_dataframe('dh_right')
r, p = stats.pearsonr(df['even'], df['odd'])
# axs[2, 1].title.set_text('Dorsal horn right')
# Lower part charts
sns.histplot(data=df, x="even", kde=True, bins=10, ax=axs[3, 1], color=dr_color)
axs[3, 1].set_xlim(my_limits)
axs[3, 1].invert_yaxis()
for side in ['top', 'right', 'bottom', 'left']:
    axs[3, 1].spines[side].set_visible(False)
axs[3, 1].get_yaxis().set_visible(False)

axs[3, 0].axis("off")
axs[2, 1].axis("off")
axs[2, 2].axis("off")
axs[2, 3].axis("off")
axs[3, 3].axis("off")

# Left part charts
sns.histplot(data=df, y='odd', kde=True, bins=10, ax=axs[2, 0], color=dr_color)
axs[2, 0].set_ylim(my_limits)
axs[2, 0].invert_xaxis()
for side in ['top', 'right', 'bottom', 'left']:
    axs[2, 0].spines[side].set_visible(False)
axs[2, 0].tick_params(axis='x', which='both', labelbottom=False, bottom=False)

# Linear regression for the middle plot
axs[2, 1].axhline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
axs[2, 1].axvline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
sns.regplot(data=df, marker='o', x='even', y='odd',
            fit_reg=True, ax=axs[2, 1], ci=0,
            scatter_kws={'alpha': 0.3, 's': 8, 'color': dr_color}, line_kws={'alpha': 1, 'color': my_color})
axs[2, 1].set_xlim(my_limits)
axs[2, 1].set_ylim(my_limits)

rect = patches.Rectangle(
    (axs[2, 1].get_position().x0, axs[2, 1].get_position().y0),
    axs[2, 1].get_position().width,
    axs[2, 1].get_position().height,
    linewidth=1.5, edgecolor=my_color, facecolor='none', transform=fig.transFigure
)
fig.patches.append(rect)

# bbox = dict(boxstyle='round', fc='LightBlue', ec='LightBlue', alpha=0.5)
axs[2, 1].text(5, -2, 'r = ' + str(np.round(r, 3)), horizontalalignment='right', color=my_color)

axs[2, 1].text(-1, 5, 'Dorsal horn right')

# 4. Dorsal horn left
df = get_dataframe('dh_left')
r, p = stats.pearsonr(df['even'], df['odd'])
# axs[2, 2].title.set_text('Dorsal horn left')

# Lower part charts
sns.histplot(data=df, x="even", kde=True, bins=10, ax=axs[3, 2], color=dl_color)
axs[3, 2].set_xlim(my_limits)
axs[3, 2].invert_yaxis()
for side in ['top', 'right', 'bottom', 'left']:
    axs[3, 2].spines[side].set_visible(False)
axs[3, 2].get_yaxis().set_visible(False)

# Right part charts
sns.histplot(data=df, y='odd', kde=True, bins=10, ax=axs[2, 3], color=dl_color)
axs[2, 3].set_ylim(my_limits)

# Linear regression for the middle plot
axs[2, 2].axhline(0, color='gray', linestyle='--', linewidth=1, zorder=1)
axs[2, 2].axvline(0, color='gray', linestyle='--', linewidth=1, zorder=1)

sns.regplot(data=df, marker='o', x='even', y='odd',
            fit_reg=True, ax=axs[2, 2], ci=0,
            scatter_kws={'alpha': 0.3, 's': 8, 'color': dl_color}, line_kws={'alpha': 1, 'color': my_color})
axs[2, 2].set_xlim(my_limits)
axs[2, 2].set_ylim(my_limits)

rect = patches.Rectangle(
    (axs[2, 2].get_position().x0, axs[2, 2].get_position().y0),
    axs[2, 2].get_position().width,
    axs[2, 2].get_position().height,
    linewidth=1.5, edgecolor=my_color, facecolor='none', transform=fig.transFigure
)
fig.patches.append(rect)

axs[2, 2].text(5, -2, 'r = ' + str(np.round(r, 3)), horizontalalignment='right', color=my_color)

axs[2, 2].text(-1, 5, 'Dorsal horn left')

fig.savefig(os.path.join(out_dir, 'even_odd_correl_all_regions.svg'), transparent=True, format='svg', bbox_inches='tight')

plt.show()
