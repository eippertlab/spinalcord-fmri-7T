#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 16:13:12 2025

@author: uhorn
"""

import numpy as np
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

main_path = '/data/pt_02661_raw/Heatpain_behav/data'
raw_path = os.path.join(main_path, 'raw_data')
output_path = os.path.join(main_path, 'derivatives')
result_path = os.path.join(output_path, 'results')

if not os.path.exists(output_path):
    os.makedirs(output_path)

if not os.path.exists(result_path):
    os.makedirs(result_path)

subjects = np.arange(1, 21)
excluded = []

new_rc_params = {"font.family": 'Arial', "font.size": 12, "font.serif": [],
                 "svg.fonttype": 'none'}
mpl.rcParams.update(new_rc_params)

long_color = [197/255, 27/255, 125/255]
short_color = [222/255, 119/255, 174/255]

df = pd.read_csv(os.path.join(raw_path, 'McGill_Resulttable.csv'))

means_short = df.filter(regex='^short').mean()
means_long = df.filter(regex='^long').mean()

means_short.index = means_short.index.str.replace('^short_', '', regex=True)
means_long.index = means_long.index.str.replace('^long_', '', regex=True)

data = pd.concat([means_short, means_long], axis=1).reset_index()
data.columns = ['Category', 'Short', 'Long']

df_melted = data.melt(id_vars='Category', var_name='Stimulus', value_name='Mean')

sns.set_context("poster")
sns.set_style("ticks")

fig, ax = plt.subplots(1, 1, figsize=(9, 8))

sns.lineplot(data=df_melted, x='Category', y='Mean', 
             hue='Stimulus', style='Stimulus', marker='o', 
             hue_order=['Long', 'Short'],
             palette=[long_color, short_color])
plt.xticks(rotation=90)
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(os.path.join(result_path, 'Questionnaire.svg'), 
             transparent=True, format='svg', bbox_inches='tight')
plt.show()