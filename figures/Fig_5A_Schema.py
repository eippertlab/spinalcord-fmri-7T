#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data were extracted from Treede et al using the WebPlotDigitizer Tool 
see automeris.io
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d

project_dir = '/data/pt_02661_raw/Heatpain/derivatives'

#%% ---------------------------------------------------------------------------
df1 = pd.read_csv(os.path.join(project_dir, 'results', 'Treede', 'Fig3B.csv'), 
                  sep=';', decimal=',', header=None, names=['Bar', 'value'])
df1['x'] = df1.index / 5

df1['y'] = df1['value'].rolling(15, center=True, min_periods=1).sum()/15

sns.set_context("poster")
sns.set_style("ticks")
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 4))
sns.lineplot(data=df1, x='x', y='y', color='black')
ax1.set_ylim([0, 10])
plt.xlabel('')
plt.ylabel('AP rate (Hz)')

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
fig1.savefig(os.path.join(project_dir, 'results', 'Treede',
                          'Fig5A_block.svg'), 
             bbox_inches='tight', transparent=True, format='svg')
plt.show()

#%% ---------------------------------------------------------------------------
df2 = pd.read_csv(os.path.join(project_dir, 'results', 'Treede', 'Fig3C.csv'), 
                  sep=';', decimal=',', header=None, names=['Bar', 'value'])
df2['x'] = df2.index / 5

df2['y'] = df2['value'].rolling(5, center=True, min_periods=1).sum()/5

sns.set_context("poster")
sns.set_style("ticks")
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 4))
sns.lineplot(data=df2, x='x', y='y', color='black')
ax2.set_ylim([0, 10])
plt.xlabel('')
plt.ylabel('AP rate (Hz)')

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig2.savefig(os.path.join(project_dir, 'results', 'Treede',
                          'Fig5A_onset.svg'), 
             bbox_inches='tight', transparent=True, format='svg')

plt.show()