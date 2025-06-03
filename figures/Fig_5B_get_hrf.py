import pandas as pd
import os
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

main_path = '/data/pt_02661_raw/Heatpain/derivatives'

design = 'onset'  # onset or block

my_color = 'black'

# pick an example subject (not the first one, there were some trigger issues)
design_file = os.path.join(main_path, 'sub-sspr02', 'func', 'run-1', 'feat_' + design + '_2mm_susan.feat', 'design.mat')
timing_file = os.path.join(main_path, 'sub-sspr02', 'func', 'run-1', 'design_block_regr.txt')

tr = 1.12

# read the HRF of onset regressor
matrix_rows = []
with open(design_file, 'r') as file:
    # skip the first 5 lines
    for _ in range(5):
        next(file)

    # read the rest of the file
    for line in file:
        row_values = line.strip().split()
        row_values = [float(value) for value in row_values]
        matrix_rows.append(row_values)

df = pd.DataFrame(matrix_rows)

# where are the block starts and ends
onset = pd.read_csv(timing_file, delimiter=' ', names=['onset', 'duration', 'value'])
rectangles = []
for o in onset['onset']:
    my_tuple = (o/tr, (o + 30)/tr)
    rectangles.append(my_tuple)

sns.set_context("poster")
sns.set_style("ticks")

mpl.rcParams['text.color'] = my_color
mpl.rcParams['xtick.color'] = my_color
mpl.rcParams['ytick.color'] = my_color
mpl.rcParams['axes.labelcolor'] = my_color
mpl.rcParams['axes.edgecolor'] = my_color

fig, ax = plt.subplots(1, 1, figsize=(12, 4))

for start, end in rectangles:
    rect = patches.Rectangle((start, 0), end - start, 1, 
                             transform=ax.get_xaxis_transform(),
                             edgecolor=None, facecolor='gray', alpha=0.2)
    ax.add_patch(rect)

# first column is onset regressor
plt.plot(df[0], color=my_color)
ax.set_xlabel('Volume')
ax.set_ylabel('HRF amplitude')
plt.subplots_adjust(bottom=0.2, left=0.1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(os.path.join(main_path, 'results', design + '_design.svg'), 
            transparent=True, format='svg', bbox_inches='tight')

plt.show()
