#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 22:52:14 2023

@author: daniel
"""
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(font_scale = 1.5)
sns.set_style("ticks")

vs_widths = []
vs_amps = []
vs = []
    
filename = './results/alpha/noise_spillover_alpha_0.062.dat'
with open(filename, 'r', encoding = 'utf-8') as f:
    fload  = json.load(f)
    res_dict = json.loads(fload)

num_syns = res_dict['num_syns']
independent_dends = res_dict['independent_dends']
trials = res_dict['trials']
vs = res_dict['vs']
vs = np.reshape(vs, (len(num_syns), trials, len(independent_dends), 500))
colors = sns.color_palette("icefire", len(num_syns))
fig_vs = plt.figure()
ax_vs = fig_vs.add_subplot(111)
for i in range(0,len(num_syns)):
    ax_vs.plot(vs[i, 20, 5, :], color = colors[i])

plt.show()