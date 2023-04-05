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
alpha = 0.086
#etas = np.arange(0.02, 0.4, 0.04)
etas = np.exp([-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1])
v12 = [-48, -40, -32, -24, -16, -8, 0, 8, 16]
eta = 0.38
mg = 1.0

vs_widths = []
vs_amps = []
vs = []

for eta in etas:
    
#    filename = './results/alpha/data_no_spillover_alpha_%.3f.dat' % alpha
    filename = './results/eta/no_spillover_eta_%.3f.dat' % eta
    with open(filename, 'r', encoding = 'utf-8') as f:
        fload  = json.load(f)
        res_dict = json.loads(fload)
    num_syns = res_dict['num_syns']
    independent_dends = res_dict['independent_dends']
    trials = res_dict['trials']
    vswidths = res_dict['vs_widths']
#    vdwidths = res_dict['vd_widths']
    vsamps = res_dict['vs_amps']
#    vdamps = res_dict['vd_amps']
    vs.append(res_dict['vs'])
    vs_widths.append(vswidths)
    vs_amps.append(vsamps)    
#    vd_widths.append(vdwidths)
#    vd_amps.append(vdamps)    

vs_widths_diff = []

vs = np.reshape(vs, (len(etas), len(num_syns), trials, len(independent_dends), 500))
vs_widths = np.reshape(vs_widths, (len(etas), len(num_syns), trials, len(independent_dends)))
vs_amps = np.max(vs, axis = (4))
d1 = vs_amps.shape[0]; d2 = vs_amps.shape[1]; d3 = vs_amps.shape[2]; d4 = vs_amps.shape[3]; 
for i in range(0,d1):
    for j in range(0,d2):
        for k in range(0,d3):
            for l in range(0,d4):
                if vs_amps[i,j,k,l] > -50:
                    vs_amps[i,j,k,l] = -50
vs_widths_mean = np.mean(vs_widths, axis = (2,3))
vs_amps_mean = np.mean(vs_amps, axis = (2))

vs_amps_diff = vs_amps[:,1:,:,:] - vs_amps[:,0:(len(num_syns)-1),:,:]
max_amps_diff = np.max(vs_amps_diff, axis = 1)
mean_max_amps_diff = np.mean(max_amps_diff, axis = (1,2))
sd_max_amps_diff = np.std(max_amps_diff, axis = (1,2))
#for d, a in zip(vs_widths, vs_amps):
#    end = len(d)
#    ddiff = np.asarray(d[1:]) - np.asarray(d[0:(end-1)])
#    adiff = np.asarray(a[1:]) - np.asarray(a[0:(end-1)])
#    vs_widths_diff.append(ddiff.tolist())
#    vs_amps_diff.append(adiff.tolist())

vs_dur_df = pd.DataFrame(vs_widths_mean, index = etas.round(2),
                  columns = [i for i in range(1, 21)] )
vs_amps_df = pd.DataFrame(vs_amps_mean[:,:,0], index = v12,
                  columns = [i for i in range(1, 21)] )

#ax1 = sns.heatmap(vs_amps_df, cmap = "icefire", cbar_kws = {'label': 'plateau width (ms)'})
#ax1.set_xlabel('cluster size')
#ax1.set_ylabel('steepness of Mg block (%)')

sns.set(font_scale = 1.5)
sns.set_style("ticks")
vmin = -78; vmax = -55
ax2 = sns.heatmap(vs_amps_df, cmap = "icefire", cbar_kws = {'label': 'soma Vm (mV)'}, 
                              vmax = vmax, vmin = vmax)
ax2.set_xlabel('cluster size')
ax2.set_ylabel('V$_{1/2}$ (mv)')

colors = sns.color_palette("icefire", len(num_syns))


fig_vs = plt.figure()
ax_vs = fig_vs.add_subplot(111)
for i in range(0,len(num_syns)):
    ax_vs.plot(vs[8, i, 40, 0, :], color = colors[i])

ax_vs.set_xlabel('t (ms)')
ax_vs.set_yticks([-80, -75, -70, -65])
ax_vs.set_ylim([-84,-60])
ax_vs.set_ylabel('soma Vm (mV)')

fig_bar = plt.figure()
ax_bar = fig_bar.add_subplot(111)
ax_bar.bar(etas, mean_max_amps_diff, yerr = sd_max_amps_diff)
plt.show()
