#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:05:58 2023

@author: daniel
"""
from neuron import h
import d1msn as msn
import pickle
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import parameters as p

def mg_block(mg, alpha, v, eta = 1/3.57, x_offset = 0):
    block = 1 / (1 + mg * eta * np.exp(-alpha * (v + x_offset)) )
    return block

def deriv_mg_block(mg, alpha, v, eta = 1/3.57, x_offset = 0):
    block = 1 / (1 + mg * eta * np.exp(-alpha * (v + x_offset)) )
    return block * (1-block) *(alpha)

params  = "./params_dMSN.json"

with open('D1_71bestFit_updRheob.pkl', 'rb') as f:
    model_sets  = pickle.load(f, encoding="latin1")

with open(params) as file:
    par = json.load(file)

cell_index = 34
variables = model_sets[cell_index]['variables'] 
cell = msn.MSN(params = params, variables = variables)

v = np.linspace(-100, 50, 15000);
alphal = np.arange(0.054, 0.076, 0.002)
eta = np.arange(0.02, 0.04, 0.002)
eta = 0.38
mg = 1.0
num_syns = range(1,21)

mbl = []
derivs = []
for alpha in alphal:
    mbl.append(mg_block(mg, alpha, v, eta))
    derivs.append(deriv_mg_block(mg, alpha, v, eta))

maximums = []
for d in derivs:
    maximums.append(max(d))
arctan_maximums = np.arctan(maximums)
norm_arctan_maximums = np.multiply(np.divide(arctan_maximums, arctan_maximums[4]), 100)
norm_arctan_maximums = np.round(norm_arctan_maximums)
print(np.degrees(np.arctan(maximums)))

# Load the array from the file
sns.set(font_scale = 1.5)
sns.set_style("ticks")
vmin = -80
vmax = -50

vs_amps_mean = np.load('./results/vs_amps_mean.npy')
figs = []
axs = []

for dend in range(0,11):
    vs_amps_df = pd.DataFrame(vs_amps_mean[:,:,dend], index = [i for i in norm_arctan_maximums],
                              columns = [i for i in range(1, 21)] )
    figs.append(plt.figure())
    axs.append(figs[-1].add_subplot(111))
    ax =  sns.heatmap(vs_amps_df, cmap = "icefire", cbar_kws = {'label': 'soma Vm (mV)'}, vmin = vmin, vmax = vmax)
    axs.append(ax)
    axs[-1].set_xlabel('cluster size')
    axs[-1].set_ylabel('steepness of Mg block (%)')
    idx = p.independent_dends[dend]
    axs[-1].set_title('Dendrite %d' % idx)
    
    sec = cell.dendlist[idx]; s = p.cluster_start_poss[dend]; e = p.cluster_end_poss[dend]
    start = h.distance(s, sec = sec)
    end = h.distance(e, sec = sec)
    L = end - start
    d = sec.diam
    print('Dendrite %d: Start = %.2f, End = %.2f, L = %.2f, d = %.2f' %( idx, start, end, L, d))

colors = sns.color_palette("icefire", len(num_syns))
plt.show()
