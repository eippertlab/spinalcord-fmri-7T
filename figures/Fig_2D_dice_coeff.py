import os
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt

main_path = '/data/pt_02661_raw/Heatpain/derivatives'
screening_path = '/data/pt_02616/data/derivatives'

color_lilac = [117/255, 112/255, 179/255]
color_lilac_bright = [188/255, 189/255, 220/255]


def dice_coefficient(image1, image2):
    intersection = np.sum(image1 * image2)
    return (2. * intersection) / (np.sum(image1) + np.sum(image2))


# dataset 2
subjects = np.arange(1, 17)
selected = pd.DataFrame([])
for isub in subjects:
    sub_id = 'sub-sspr' + '{:02d}'.format(isub)
    # get the EPI segmentation
    bold = os.path.join(main_path, sub_id, 'func', 'second_level', 'MOCO_MEAN_deepseg_manual.nii.gz')
    bold_img = nib.load(bold).get_fdata()
    bold_img_bin = (bold_img > 0).astype(int)
    # get the anatomical segmentation (in EPI space)
    megre = os.path.join(main_path, sub_id, 'anat_t2s', 'MEGRE_RMS_MOCO_MEAN_seg_EPIspace_translation.nii.gz')
    megre_img = nib.load(megre).get_fdata()
    megre_img_bin = (megre_img > 0).astype(int)

    DC_bold = dice_coefficient(bold_img_bin, megre_img_bin)
    df = pd.DataFrame([DC_bold], columns=['dice'])
    df['subject'] = sub_id
    selected = pd.concat([selected, df])

# screening dataset
subjects = np.arange(1, 63)
# two subjects do not have functional images
subjects = subjects[subjects != 29]
subjects = subjects[subjects != 60]
# this one had bad artifacts
subjects = subjects[subjects != 56]
screening = pd.DataFrame([])
for isub in subjects:
    sub_id = 'sub-' + '{:03d}'.format(isub)
    # get the EPI segmentation
    bold = os.path.join(screening_path, sub_id, 'func', 'sosGREp3', 'MOCO_MEAN_deepseg_manual.nii.gz')
    bold_img = nib.load(bold).get_fdata()
    bold_img_bin = (bold_img > 0).astype(int)
    # get the anatomical segmentation (in EPI space)
    megre = os.path.join(screening_path, sub_id, 'anat_t2s', 'MEGRE_RMS_MOCO_MEAN_seg_EPIspace_translation.nii.gz')
    megre_img = nib.load(megre).get_fdata()
    megre_img_bin = (megre_img > 0).astype(int)

    DC_bold = dice_coefficient(bold_img_bin, megre_img_bin)
    df = pd.DataFrame([DC_bold], columns=['dice'])
    df['dataset'] = ['screening']
    df['subject'] = sub_id
    screening = pd.concat([screening, df])

selected.reset_index(drop=True, inplace=True)
screening.reset_index(drop=True, inplace=True)

sns.set_context("poster")
sns.set_style("ticks")
jitter = 0.1

# order: violinplot, boxplot, raw data, raw data, boxplot, violinplot
positions = [0.5, 0.75, 1.25, 2.25, 2.75, 3]
width_violin = 0.8
width_box = 0.2

screening_x_jitter = pd.DataFrame(np.random.normal(loc=positions[2], scale=jitter, size=(len(screening), 1)),
                                  columns=['jitter'])
selected_x_jitter = pd.DataFrame(np.random.normal(loc=positions[3], scale=jitter, size=(len(selected), 1)),
                                 columns=['jitter'])

fig, ax = plt.subplots(1, 1, figsize=(6, 10), sharey=True)
# fig.suptitle("Dice coefficient", fontsize=30)

v1 = ax.violinplot(screening['dice'], positions=[positions[0]], showmeans=False, showextrema=False,
                   showmedians=False, side="low", widths=width_violin)  # `side` param requires matplotlib 3.9+
for pc in v1['bodies']:
    pc.set_facecolor(color_lilac_bright)
    pc.set_edgecolor('black')
    pc.set_linewidth(1.2)
    pc.set_alpha(0.7)

a = sns.boxplot(ax=ax, x=np.repeat(positions[1], len(screening)),
                y=screening['dice'], color=color_lilac_bright, width=width_box, native_scale=True, showfliers=False)
for patch in a.patches:
    r, g, b, al = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.7))

c = sns.boxplot(ax=ax, x=np.repeat(positions[4], len(selected)),
                y=selected['dice'], color=color_lilac, width=width_box, native_scale=True, showfliers=False)
for patch in c.patches:
    r, g, b, al = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.7))

v2 = ax.violinplot(selected['dice'], positions=[positions[5]], showmeans=False, showextrema=False,
                   showmedians=False, side="high", widths=width_violin)
for pc in v2['bodies']:
    pc.set_facecolor(color_lilac)
    pc.set_edgecolor('black')
    pc.set_linewidth(1.2)
    pc.set_alpha(0.7)

df_screening = pd.DataFrame({'jitter': screening_x_jitter['jitter'], 'dice': screening['dice']})
df_selected = pd.DataFrame({'jitter': selected_x_jitter['jitter'], 'dice': selected['dice']})
sns.scatterplot(ax=ax, data=df_screening, x='jitter', y='dice', color=color_lilac_bright, zorder=100,
                edgecolor='black', alpha=0.6, linewidth=1.2)
sns.scatterplot(ax=ax, data=df_selected, x='jitter', y='dice', color=color_lilac, zorder=100,
                edgecolor='black', alpha=0.6, linewidth=1.2)

# these subjects were the ones selected
# sspr01 = screening 001
ax.plot([screening_x_jitter.loc[0, 'jitter'], selected_x_jitter.loc[0, 'jitter']],
        [screening.loc[0, 'dice'], selected.loc[0, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr02 = screening 015
ax.plot([screening_x_jitter.loc[14, 'jitter'], selected_x_jitter.loc[1, 'jitter']],
        [screening.loc[14, 'dice'], selected.loc[1, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr03 = screening 003
ax.plot([screening_x_jitter.loc[2, 'jitter'], selected_x_jitter.loc[2, 'jitter']],
        [screening.loc[2, 'dice'], selected.loc[2, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr04 = screening 032 (=031 because of exclusion of 29)
ax.plot([screening_x_jitter.loc[30, 'jitter'], selected_x_jitter.loc[3, 'jitter']],
        [screening.loc[30, 'dice'], selected.loc[3, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr05 = screening 009
ax.plot([screening_x_jitter.loc[8, 'jitter'], selected_x_jitter.loc[4, 'jitter']],
        [screening.loc[8, 'dice'], selected.loc[4, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr06 = screening 017
ax.plot([screening_x_jitter.loc[16, 'jitter'], selected_x_jitter.loc[5, 'jitter']],
        [screening.loc[16, 'dice'], selected.loc[5, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr07 = screening 022
ax.plot([screening_x_jitter.loc[21, 'jitter'], selected_x_jitter.loc[6, 'jitter']],
        [screening.loc[21, 'dice'], selected.loc[6, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr08 = screening 019
ax.plot([screening_x_jitter.loc[18, 'jitter'], selected_x_jitter.loc[7, 'jitter']],
        [screening.loc[18, 'dice'], selected.loc[7, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr09 = screening 012
ax.plot([screening_x_jitter.loc[11, 'jitter'], selected_x_jitter.loc[8, 'jitter']],
        [screening.loc[11, 'dice'], selected.loc[8, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr10 = screening 048 (=047 because of exclusion of 29)
ax.plot([screening_x_jitter.loc[46, 'jitter'], selected_x_jitter.loc[9, 'jitter']],
        [screening.loc[46, 'dice'], selected.loc[9, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr11 = screening 043 (=042 because of exclusion of 29)
ax.plot([screening_x_jitter.loc[41, 'jitter'], selected_x_jitter.loc[10, 'jitter']],
        [screening.loc[41, 'dice'], selected.loc[10, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr12 = screening 041 (=040 because of exclusion of 29)
ax.plot([screening_x_jitter.loc[39, 'jitter'], selected_x_jitter.loc[11, 'jitter']],
        [screening.loc[39, 'dice'], selected.loc[11, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr13 = screening 046 (=045 because of exclusion of 29)
ax.plot([screening_x_jitter.loc[44, 'jitter'], selected_x_jitter.loc[12, 'jitter']],
        [screening.loc[44, 'dice'], selected.loc[12, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr14 = screening 057 (=055 because of exclusion of 29 and 56)
ax.plot([screening_x_jitter.loc[54, 'jitter'], selected_x_jitter.loc[13, 'jitter']],
        [screening.loc[54, 'dice'], selected.loc[13, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr15 = screening 059 (=057 because of exclusion of 29 and 56)
ax.plot([screening_x_jitter.loc[56, 'jitter'], selected_x_jitter.loc[14, 'jitter']],
        [screening.loc[56, 'dice'], selected.loc[14, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
# sspr16 = screening 058 (=056 because of exclusion of 29 and 56)
ax.plot([screening_x_jitter.loc[55, 'jitter'], selected_x_jitter.loc[15, 'jitter']],
        [screening.loc[55, 'dice'], selected.loc[15, 'dice']],
        color='gray', linewidth=1.0, linestyle='-', zorder=-1)
ax.set_xticks([0.5, 3.5])
ax.set_xticklabels(['screening', 'selection'])
ax.tick_params(axis='x', which='major', labelsize=25)
ax.set_ylim(0.6, 1)
sns.despine()
ax.set_xlabel('Dataset', fontsize=25)
ax.set_ylabel('Dice coefficient', fontsize=25)
plt.subplots_adjust(bottom=0.15, left=0.2)
fig.savefig(os.path.join(main_path, 'results', 'dice_coeff_dataset2_with_connections.svg'), 
            transparent=True, bbox_inches='tight', format='svg')
plt.show()


