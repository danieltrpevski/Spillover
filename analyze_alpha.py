# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:48:41 2016

@author: daniel
"""
import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def mg_block(mg, alpha, v, eta = 1/3.57, x_offset = 0):
    block = 1 / (1 + mg * eta * np.exp(-alpha * (v + x_offset)) )
    return block

def deriv_mg_block(mg, alpha, v, eta = 1/3.57, x_offset = 0):
    block = 1 / (1 + mg * eta * np.exp(-alpha * (v + x_offset)) )
    return block * (1-block) *(alpha)

v = np.linspace(-100, 50, 15000);
alphal = np.arange(0.054, 0.076, 0.002)
eta = np.arange(0.02, 0.04, 0.002)
eta = 0.38
mg = 1.0

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

vs_widths = []
vs_amps = []
vs = []
alphal = np.arange(0.054, 0.076, 0.002)
for alpha in alphal:
    
#    filename = './results/alpha/data_no_spillover_alpha_%.3f.dat' % alpha
    filename = './results/alpha/noise_spillover_alpha_%.3f.dat' % alpha
    with open(filename, 'r', encoding = 'utf-8') as f:
        fload  = json.load(f)
        res_dict = json.loads(fload)
    num_syns = res_dict['num_syns']
    independent_dends = res_dict['independent_dends']
    trials = res_dict['trials']
    vs.append(res_dict['vs'])

vs_widths_diff = []

vs = np.reshape(vs, (len(alphal), len(num_syns), trials, len(independent_dends), 500))
vs_amps = np.max(vs, axis = (4))
for a in range(0, len(alphal)):
    for b in range(0, len(num_syns)):
        for c in range(0, trials):
            for d in range(0, len(independent_dends)):
                if (vs_amps[a, b, c, d] > -50):
                    vs_amps[a, b, c, d] = -50

vs_amps_mean = np.mean(vs_amps, axis = (2))

vs_amps_diff = vs_amps[:,1:,:,:] - vs_amps[:,0:(len(num_syns)-1),:,:]
max_amps_diff = np.max(vs_amps_diff, axis = 1)
mean_max_amps_diff = np.mean(max_amps_diff, axis = (1,2))
sd_max_amps_diff = np.std(max_amps_diff, axis = (1,2))

vs_amps_df = pd.DataFrame(vs_amps_mean[:,:,0], index = [i for i in norm_arctan_maximums],
                  columns = [i for i in range(1, 21)] )

sns.set(font_scale = 1.5)
sns.set_style("ticks")
vmin = -78
vmax = -50
ax2 = sns.heatmap(vs_amps_df, cmap = "icefire", cbar_kws = {'label': 'soma Vm (mV)'}, vmin = vmin, vmax = vmax)
ax2.set_xlabel('cluster size')
ax2.set_ylabel('steepness of Mg block (%)')

colors = sns.color_palette("icefire", len(num_syns))

fig_vs = plt.figure()
ax_vs = fig_vs.add_subplot(111)
for i in range(0,len(num_syns)):
    ax_vs.plot(vs[4, i, 40, 0, :], color = colors[i])

ax_vs.set_xlabel('t (ms)')
ax_vs.set_yticks([-80, -75, -70, -65, -50])
ax_vs.set_ylim([-84,-50])
ax_vs.set_ylabel('soma Vm (mV)')

fig_bar = plt.figure()
ax_bar = fig_bar.add_subplot(111)
ax_bar.bar(norm_arctan_maximums, mean_max_amps_diff, yerr = sd_max_amps_diff)
plt.show()

# Save the array to a file
np.save('vs_amps_mean.npy', vs_amps_mean)