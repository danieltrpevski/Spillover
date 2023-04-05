# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:17:15 2019

@author: daniel
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def sigmoid(x, x_offset, alpha, Mg, eta):
    s = 1 / (1 + Mg*eta* np.exp(-alpha * (x - x_offset)) )
    return s

Erev_AMPA = 0
Erev_NMDA = 15
Erev_GABA = -80
gAMPA = 0.0
gGABA = 0.4
gNMDA = 0.5
gKir = 0.0
Erev_K = -80
alpha = 0.08
eta = 0.25
Mg = 1.0
alphas = np.arange(0.054, 0.096, 0.004)
etas = np.exp([-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]) #np.arange(0.02, 0.4, 0.04)
#etas = np.exp(np.arange(-80, 80, 20))
v = np.arange(-100, 0,0.01)
#kinf = sigmoid(v, -102, -1/13, 1, 1)

sns.set(font_scale = 1.5)
sns.set_style("ticks")
fa = plt.figure()
ax_a = fa.add_subplot(111)

colors = sns.cubehelix_palette(len(alphas))
for a,c in zip(alphas, colors):    
    block = sigmoid(v, 0, a, Mg, eta)
    i = gGABA*(v - Erev_GABA)+ gAMPA*(v- Erev_AMPA) + np.multiply(gNMDA*(v-Erev_NMDA), block) #+ gKir*kinf*(v - Erev_K)
    ax_a.plot(v, i, color = c)
ax_a.plot(v, np.zeros(v.shape), '--', color = 'grey')
ax_a.set_xlabel('V (mV)')
ax_a.set_ylabel('I (pA)')

fe = plt.figure()
ax_e = fe.add_subplot(111)
alpha = 0.086
colors = sns.cubehelix_palette(len(etas), rot = -0.5)
for e,c in zip(etas, colors):    
    block = sigmoid(v, 0, alpha, Mg, e)
    i = gGABA*(v - Erev_GABA)+ gAMPA*(v- Erev_AMPA) + np.multiply(gNMDA*(v-Erev_NMDA), block) #+ gKir*kinf*(v - Erev_K)
    ax_e.plot(v, i, color = c)
ax_e.plot(v, np.zeros(v.shape), '--', color = 'grey')
ax_e.set_xlabel('V (mV)')
ax_e.set_ylabel('I (pA)')
plt.show()
