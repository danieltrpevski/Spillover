# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:48:41 2016

@author: daniel
"""

from neuron import h
import d1msn as msn
import spillover_experiment as se
import pickle
import parameters as p
import json 

# --- 1. Create a cell and other useful stuff

params  = "./params_dMSN.json"

with open('D1_71bestFit_updRheob.pkl', 'rb') as f:
    model_sets  = pickle.load(f, encoding="latin1")

with open(params) as file:
    par = json.load(file)

cell_index = 34
variables = model_sets[cell_index]['variables'] 
cell = msn.MSN(params = params, variables = variables)

for d in p.input_dends:
    cell.dendlist[d].nseg *=5

for sec in cell.dendlist:
    print(sec.name(), "%f, %f, %f, %f, d = %.2f" % (h.distance(1.0, sec = sec), 
                                      h.distance(0, sec = sec), 
                                      h.distance(0.45, sec = sec) - h.distance(0.3, sec = sec),
                                      h.distance(0.05, sec = sec),
                                      sec.diam))

for sec in cell.somalist:
    print(sec.name(), "%f, %f, %f, d = %.2f" % (h.distance(1, sec = sec), 
                                      h.distance(0, sec = sec), 
                                      h.distance(1, sec = sec) - h.distance(0, sec = sec),
                                      sec.diam))

# --- 2. Insert stimulation to cell

dend_record_list = [12]
dend_stim_list = []                    
plateau_cluster_list = [12]
inhibitory_cluster_dict = {'loc': [53], 
                        'pos': [0.85], 
                        'start': [p.inhibitory_burst_start ],
                        'end': [p.inhibitory_burst_end ] }           


#istim = h.IClamp(cell.somalist[0](0.5))
#istim.dur = 25
#istim.amp = 0.9
#istim.delay = 100
#
#istim2 = h.IClamp(cell.somalist[0](0.5))
#istim2.dur = 2000
#istim2.amp = 0.15

ex = se.Spillover_Experiment('record_ca', cell)
#ex.insert_synapses('noise_SPN')
ex.insert_synapses('my_spillover', plateau_cluster_list, deterministic = 0, 
                   num_syns = 30, add_spine = 1, on_spine = 0)
ex.set_up_recording(dend_record_list)
ex.simulate()
ex.plot_results()