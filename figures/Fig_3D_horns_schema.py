import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import nibabel as nib

project_dir = '/data/pt_02661_raw/Heatpain/derivatives'
out_dir = '/data/pt_02661_raw/Heatpain/derivatives/results/spatial_specificity'

mask_dir = os.path.join(project_dir, 'group_analysis', 'both', 'masks')

# get cord mask
cord = os.path.join(mask_dir, 'C6_cut.nii.gz')
imgcord = nib.load(cord)
header = imgcord.header
print(header.get_data_shape())
print(header)
cord_data = imgcord.get_fdata()

slice_cord = cord_data[:, :, 47]
slice_cord_array = pd.DataFrame(columns=["x", "y"])
slice_cord_t = np.transpose(slice_cord)
for x in range(slice_cord.shape[0]):
    for y in range(slice_cord.shape[1]):
        if slice_cord[x, y] == 1:
            slice_cord_array = pd.concat([slice_cord_array,
                                          pd.DataFrame({"x": float(x), "y": float(y)}, index=[0])],
                                         ignore_index=True)
slice_cord_array["x"] = pd.to_numeric(slice_cord_array['x'], downcast='float')
slice_cord_array["y"] = pd.to_numeric(slice_cord_array['y'], downcast='float')

# get horn masks
# dhl
dhl_file = os.path.join(mask_dir, 'dh_left_C6_cut.nii.gz')
imgdhl = nib.load(dhl_file)
header = imgdhl.header
print(header.get_data_shape())
print(header)
dhl_data = imgdhl.get_fdata()

dhl = dhl_data[:, :, 47]
dhl_df = pd.DataFrame(columns=["x", "y"])
for x in range(dhl.shape[0]):
    for y in range(dhl.shape[1]):
        if dhl[x, y] == 1:
            dhl_df = pd.concat([dhl_df, pd.DataFrame({"x": float(x), "y": float(y)}, index=[0])],
                               ignore_index=True)
dhl_df["x"] = pd.to_numeric(dhl_df['x'], downcast='float')
dhl_df["y"] = pd.to_numeric(dhl_df['y'], downcast='float')

# dhr
dhr_file = os.path.join(mask_dir, 'dh_right_C6_cut.nii.gz')
imgdhr = nib.load(dhr_file)
header = imgdhr.header
print(header.get_data_shape())
print(header)
dhr_data = imgdhr.get_fdata()

dhr = dhr_data[:, :, 47]
dhr_df = pd.DataFrame(columns=["x", "y"])
for x in range(dhr.shape[0]):
    for y in range(dhr.shape[1]):
        if dhr[x, y] == 1:
            dhr_df = pd.concat([dhr_df, pd.DataFrame({"x": float(x), "y": float(y)}, index=[0])],
                               ignore_index=True)
dhr_df["x"] = pd.to_numeric(dhr_df['x'], downcast='float')
dhr_df["y"] = pd.to_numeric(dhr_df['y'], downcast='float')

# vhr
vhr_file = os.path.join(mask_dir, 'vh_right_C6_cut.nii.gz')
imgvhr = nib.load(vhr_file)
vhr_data = imgvhr.get_fdata()

vhr = vhr_data[:, :, 47]
vhr_df = pd.DataFrame(columns=["x", "y"])
for x in range(vhr.shape[0]):
    for y in range(vhr.shape[1]):
        if vhr[x, y] == 1:
            vhr_df = pd.concat([vhr_df, pd.DataFrame({"x": float(x), "y": float(y)}, index=[0])],
                               ignore_index=True)
vhr_df["x"] = pd.to_numeric(vhr_df['x'], downcast='float')
vhr_df["y"] = pd.to_numeric(vhr_df['y'], downcast='float')

# vhl
vhl_file = os.path.join(mask_dir, 'vh_left_C6_cut.nii.gz')
imgvhl = nib.load(vhl_file)
vhl_data = imgvhl.get_fdata()

vhl = vhl_data[:, :, 47]
vhl_df = pd.DataFrame(columns=["x", "y"])
for x in range(vhl.shape[0]):
    for y in range(vhl.shape[1]):
        if vhl[x, y] == 1:
            vhl_df = pd.concat([vhl_df, pd.DataFrame({"x": float(x), "y": float(y)}, index=[0])],
                               ignore_index=True)
vhl_df["x"] = pd.to_numeric(vhl_df['x'], downcast='float')
vhl_df["y"] = pd.to_numeric(vhl_df['y'], downcast='float')

mpl.rcParams['pdf.fonttype'] = 42
sns.set(style="white")
color_c6 = "orangered"
with sns.plotting_context('paper', font_scale=2):
    fig, ax = plt.subplots(figsize=(4, 4))
    # add contour
    sns.kdeplot(data=slice_cord_array, x="x", y="y", levels=1, bw_method=0.2,
                color='k', linewidths=2, ax=ax)
    sns.kdeplot(data=dhl_df, x="x", y="y", levels=1, bw_method=0.4, color="k", linewidths=2, ax=ax)
    sns.kdeplot(data=dhr_df, x="x", y="y", levels=1, bw_method=0.4, color="k", linewidths=2, ax=ax)
    sns.kdeplot(data=vhl_df, x="x", y="y", levels=1, bw_method=0.4, color="k", linewidths=2, ax=ax)
    sns.kdeplot(data=vhr_df, x="x", y="y", levels=1, bw_method=0.4, color="k", linewidths=2, ax=ax)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.xaxis.set_tick_params(labelbottom=False)
ax.yaxis.set_tick_params(labelleft=False)
ax.set(xlabel=None)
ax.set(ylabel=None)

mpl.axes.Axes.set_aspect(ax, aspect="equal")

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'schema_horns.pdf'),
            bbox_inches='tight', transparent=True)
plt.show()
