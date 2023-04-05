# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:44:21 2017

@author: daniel
"""

import parameters as p
from neuron import h
from random import randint
import time
import sys

class Experiment(object):
    
    def __init__(self):
        self.inc = []; self.istim = [] 
        self.enc = []; self.estim = []
        self.recorder_nc = []; self.recorder_tvec = []
        self.soma_recorder_nc = []; self.soma_recorder_tvec = []
        self.ramp_enc = []; self.ramp_estim = []
        self.presyn = h.Section()
                
    def insert_synapses(self): 
        """Inserts synaptic input."""
        raise NotImplementedError("insert_synapses() is not implemented.")               
    
    def set_up_recording(self):
        """Sets up recording of cell properties."""
        raise NotImplementedError("set_up_recording() is not implemented.")
                        
    def simulate(self):
        """Runs the simulation."""
        raise NotImplementedError("simulate() is not implemented.")
         
    def plot_results(self):
        """Plots results for this experiment."""
        raise NotImplementedError("plot_results() is not implemented.")
        
    def add_input_generator(self, syn, syntype, freq_multiplier = 1.0):
        if syntype in ['expsyn', 'exp2syn']:
            syn.stim.append(h.SpikeGenerator(0.5, sec = self.presyn))
            syn.stim[-1].seed(int(time.time() + randint(1,10**7)))
            syn.stim[-1].noise = 1
            syn.stim[-1].fast_invl = p.e_interval/freq_multiplier
            syn.stim[-1].slow_invl = 0
            
            syn.nc.append(h.NetCon(self.estim[-1], syn.obj))
            syn.nc[-1].delay = 0
            syn.nc[-1].weight[0] = p.g_expsyn_max
        
        elif syntype in ['inhexpsyn', 'inhexp2syn']:
            syn.stim.append(h.SpikeGenerator(0.5, sec = self.presyn)) 
            syn.stim[-1].seed(int(time.time() + randint(1,10**7)))
            syn.stim[-1].fast_invl = p.i_interval/freq_multiplier
            syn.stim[-1].slow_invl = 0
            syn.stim[-1].noise = 1

            syn.nc.append(h.NetCon(syn.stim[-1], syn.obj))      
            syn.nc[-1].delay = 0
            syn.nc[-1].weight[0] = p.g_inhexpsyn_max

    def synapse_distribution(self, celltype, dend):
        if celltype == 'MSN':
            # Synapses according to Cheng et al. Experimental Neurobiology, 147:287-298 (1997)                
            dist = h.distance(1, sec = dend)
            if dist <= 20:
                exc_mean = 7.7; exc_sem = 1.1;
                inh_mean = 1.71; inh_sem = 0.37;
            elif dist <= 40:
                exc_mean = 22.6; exc_sem = 1.5;
                inh_mean = 1.71; inh_sem = 0.37;
            elif dist <= 60:
                exc_mean = 32.1; exc_sem = 1.1;
                inh_mean = 1.71; inh_sem = 0.37;
            elif dist <= 80:
                exc_mean = 31.3; exc_sem = 1.0;
                inh_mean = 1.71; inh_sem = 0.37;
            elif dist <= 100:
                exc_mean = 32.5; exc_sem = 1.0;
                inh_mean = 1.71; inh_sem = 0.37;
            elif dist <= 120:
                exc_mean = 29.4; exc_sem = 0.8;
                inh_mean = 1.71; inh_sem = 0.37;
            else:
                exc_mean = 25.4; exc_sem = 0.8; # Not sure here!
                inh_mean = 1.71; inh_sem = 0.37;
                
        return [exc_mean, exc_sem, inh_mean, inh_sem]
        
    def helper_insert(self, syntype, pos, dend, freq_multiplier):  
        if syntype in ['expsyn', 'exp2syn',
                       'inhexpsyn', 'inhexp2syn', 'glutamate']:
            syn = self.cell.insert_synapse(syntype, dend, pos)
            self.add_input_generator(syn, syntype, freq_multiplier)
            
    def connect_input_generator(self, syn, syntype, gen, delay = 0):
        if syntype == 's_exglu':
            self.enc.append(h.NetCon(gen, syn.obj))
            self.enc[-1].delay = 0
            self.enc[-1].weight[0] = 1.0
        elif syntype == 't_exglu':
            self.exnc.append(h.NetCon(gen, syn))
            self.exnc[-1].delay = delay
            self.exnc[-1].weight[0]= p.weight
        else:
            syn.nc.append(h.NetCon(gen, syn.obj))
            syn.nc[-1].delay = delay
            syn.stim.append(gen)
            if syntype in ['expsyn', 'exp2syn']:
                syn.nc[-1].weight[0] = p.g_expsyn_max
            elif syntype == 'nmda_plateau':
                syn.nc[-1].weight[0] = p.gNMDAmax_plateau
            elif syntype == 'tmGlut':
                syn.nc[-1].weight[0] = p.gAMPAmax
            elif syntype == 'glutamate':
                syn.nc[-1].weight[0] = p.g_expsyn_max
            elif syntype == 'glutamate_plateau':
                syn.nc[-1].weight[0] = p.gAMPAmax_plateau
            elif syntype in ['AMPA', 'NMDA', 'AMPA_pf', 'NMDA_pf', 'NMDA_stp','AMPA_stp']:
                syn.nc[-1].weight[0] = p.weight            
            elif syntype in ['adaptive_glutamate',
                             'adaptive_glutamate_test', 'glutamate_ica_nmda',
                             'adaptive_glutamate_hom', 
                             'glutamate_xor_test', 'adaptive_shom_AMPA_stp', 'adaptive_shom_NMDA_stp',
                             'generalized_rule', 'generalized_rule_dist', 'adaptive_AMPA', 'adaptive_NMDA', 
                             'AMPA_test', 'NMDA_test', 'adaptive_hom_AMPA', 'adaptive_hom_NMDA',
                             'adaptive_NMDAe','adaptive_shom_AMPA', 'adaptive_pf_AMPA', 'adaptive_pf_NMDA',
                             'adaptive_shom_NMDA','adaptive_glutamate_shom','adaptive_cshom_AMPA', 
                             'adaptive_cshom_NMDA','adaptive_glutamate_cshom', 'adaptive_my_shom_NMDA',
                             'adaptive_sAMPA', 'adaptive_sNMDA', 'adaptive_sglutamate', 
                             'NMDAe', 'adaptive_zahra_NMDA', 'adaptive_zahra_AMPA']:
                syn.nc[-1].weight[0] = 1.0
            else:
                print("Synapse model not available in connect_input_genetrator().")
                sys.exit(-1)
                        
        # ADD THE REST OF THE SYNAPSE TYPES HERE !            

#        if 'nmda' in syn.hname():
#            self.enc[-1].weight[0] = 1          

    def event_recorder(self, source):
        self.recorder_nc.append(h.NetCon(source, None))
        self.recorder_tvec.append(h.Vector())
        self.recorder_nc[-1].record(self.recorder_tvec[-1])
        return self.recorder_tvec[-1]
        
    def spike_recorder(self, source, pos, threshold):
        """
        Input arguments:
        
        source - the section of a cell that is to be recorded
        pos - the position on the section (between 0 and 1)
        threshold - the value of the voltage above which signifies a spike has occured

        """
        self.soma_recorder_nc.append(h.NetCon(source(pos)._ref_v, None, sec = source))
        self.soma_recorder_nc[-1].threshold = threshold
            
        self.soma_recorder_tvec.append(h.Vector())
        self.soma_recorder_nc[-1].record(self.soma_recorder_tvec[-1])
        