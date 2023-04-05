# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 16:36:50 2017

@author: daniel
"""

from neuron import h
import experiment as e
import parameters as p
import time
import random as rnd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns

class Spillover_Experiment(e.Experiment):
    
    def __init__(self, exptype, cell, dendstatobj = [], training_mode = p.training_mode):
        super(Spillover_Experiment, self).__init__()
        self.exptype = exptype
        self.cell = cell
        self.dendstatobj = dendstatobj
        self.training_mode = training_mode
        self.spike_flags = []
        self.exglusec = []
        self.exglu = []
        self.exnc = []
        if (str(type(cell))).find('MSN') != -1:
            self.celltype = 'MSN'
                
        elif (str(type(cell))).find('L5PC') != -1:
            self.celltype = 'L5PC'            
            self.insert_synapses('L5PC')
        
    def insert_synapses(self, syntype, syn_loc = [], deterministic = 0, 
                        num_syns = p.distributed_input_size, add_spine = 0, on_spine = 0):
        if syntype in ['expsyn',  'inhexpsyn']:
            if syntype in ['expsyn']:    
                num_syns = p.distributed_input_size
            elif syntype == 'inhexpsyn':
                num_syns = p.inhibitory_cluster_size
            
            if syn_loc == []:
                for i in range (0,num_syns):
                    syn_loc.append(rnd.randint(0,len(self.cell.dendlist)-1))
             
            for loc in syn_loc:
                if syntype == 'expsyn_hom':
                    self.cell.dendlist[loc].insert('Hom_Cai')
                syn = self.cell.insert_synapse(syntype, self.cell.dendlist[loc], p.pos, add_spine = add_spine, on_spine = on_spine)                
                self.add_input_generator(syn, syntype)
                
        elif syntype == 'input_syn':
            syn_loc = []            
            for i in range(0, num_syns):
                syn_loc.append([rnd.randint(0, len(self.cell.dendlist)-1), rnd.uniform(0,1)])
            
            if deterministic == 1:
                spike_time = []
                for i in range(0, len(syn_loc)):
                    spike_time.append(rnd.uniform(p.distributed_input_start, p.distributed_input_end))
            
            counter = 0
            for loc in syn_loc:
                counter += 1
                syn1 = self.cell.insert_synapse('AMPA', self.cell.dendlist[loc[0]], loc[1], 
                           add_spine = 0, on_spine = 0)
                syn2 = self.cell.insert_synapse('NMDA', self.cell.dendlist[loc[0]], loc[1], 
                                           add_spine = 0, on_spine = 0)
#                 syn = self.cell.insert_synapse('glutamate', self.cell.dendlist[loc[0]], loc[1], add_spine = add_spine, on_spine = on_spine)
                if deterministic == 1:                
#                    self.add_input_generator(syn, syntype, deterministic = deterministic, numsyn = counter, tstart = spike_time[counter-1])
                    self.add_input_generator(syn1, 'AMPA', deterministic = deterministic, numsyn = counter, tstart = spike_time[counter-1])
                    self.connect_input_generator(syn2, 'NMDA', syn1.stim[-1]) 
                else:
#                    self.add_input_generator(syn, syntype, deterministic = deterministic, numsyn = counter)
                    self.add_input_generator(syn1, 'AMPA', deterministic = deterministic, numsyn = counter)
                    self.connect_input_generator(syn2, 'NMDA', syn1.stim[-1]) 
                    
        elif syntype == 'ramp':
            for loc in syn_loc:
                syn = self.cell.insert_synapse('expsyn',self.cell.dendlist[loc],  p.pos, add_spine = add_spine, on_spine = on_spine)
                self.add_input_generator(syn, syntype)
            
        elif syntype == 'noise_SPN':
            for dend in self.cell.dendlist:
                # Synapses according to Cheng et al. Experimental Neurobiology, 147:287-298 (1997)                
                unit_length = 20.0 # Values reported in units per 20 microns in the study
                [exc_mean, exc_sem, inh_mean, inh_sem] = self.synapse_distribution(self.celltype, dend)
                # Insert excitatory synapses in this section, 
                # This contains both an exponential and an NMDA synapse.                
                syntype = 'glutamate'
                if dend.nseg <= exc_mean:
                    freq_multiplier = exc_mean/dend.nseg
                    step = 1/dend.nseg
                    for i in range(0, dend.nseg):
                        pos = (i + (i+1))*step/2
                        self.helper_insert(syntype, pos, dend, freq_multiplier)
                else:                        
                    num_exc_syn = int(dend.L/unit_length * rnd.gauss(exc_mean,exc_sem))
                    freq_multiplier = 1.0                     
                    for i in range(0,num_exc_syn):
                        pos = rnd.uniform(0,1)
                        self.helper_insert(syntype, pos, dend, freq_multiplier)                        
                
                # Insert inhibitory synapses in this section                
                syntype = 'inhexp2syn'
                if dend.nseg <= inh_mean:
                    freq_multiplier = inh_mean/dend.nseg
                    step = 1/dend.nseg
                    for i in range(0, dend.nseg):
                        pos = (i + (i+1))*step/2
                        self.helper_insert(syntype, pos, dend, freq_multiplier)
                else:        
                    num_inh_syn = int(dend.L/unit_length * rnd.gauss(inh_mean,inh_sem))                     
                    freq_multiplier = 1.0                    
                    for i in range(0,num_inh_syn):
                        pos = rnd.uniform(0,1)
                        self.helper_insert(syntype, pos, dend, freq_multiplier)  
        
        elif syntype in ['plateau_cluster', 'inhexpsyn_plateau', 'generalized_rule', 
                         'spillover', 'spillover_test', 'no_spillover', 'my_spillover',
                          'no_spillover_stp', 'my_spillover_stp']:
            if syntype == 'plateau_cluster':
                syntype = 'tmGlut'
            elif syntype == 'spillover':
                syntype = 'adaptive_glutamate_test'
            elif syntype == 'my_spillover' or syntype == 'my_spillover_stp':
                self.exglusec.append(h.Section(name = 'exglusec%d' % len(self.exglusec)))
                self.exglu.append(h.IntFire1(self.exglusec[-1](0.5)))
                self.exglu[-1].tau = p.exglu_tau
                self.exglu[-1].refrac = p.session_length
#                h.setpointer(h._ref_stimulus_flag, 'stimulus_flag', self.exglusec[-1].exglu)
            for num,loc in enumerate(syn_loc): 
                syn_step = 1.0/num_syns
                cluster_start_pos = p.cluster_start_poss[p.input_dends.index(loc)]
                cluster_end_pos = p.cluster_end_poss[p.input_dends.index(loc)]
                for i in range(0, num_syns):
#                    pos = cluster_start_pos + (cluster_end_pos - cluster_start_pos)*i*syn_step
                    pos = cluster_end_pos - (cluster_end_pos - cluster_start_pos)*i*syn_step
                    if syntype == 'inhexpsyn_plateau':
                        if deterministic == 0:
                            start = rnd.uniform(p.inhibitory_burst_start, p.inhibitory_burst_end)
                        elif deterministic == 1:
                            start = p.inhibitory_burst_start + i*p.deterministic_interval    
                        pos = p.inh_cluster_start_pos + (p.inh_cluster_end_pos - p.inh_cluster_start_pos)*i*syn_step
                    else:
                        if deterministic == 0:
                            start = rnd.uniform(p.plateau_burst_start, p.plateau_burst_end)
                        elif deterministic == 1:
                            start = p.plateau_burst_start + i*p.deterministic_interval                                
                    if not (syntype in ['spillover_test', 'no_spillover', 'my_spillover', 'no_spillover_stp',
                                        'my_spillover_stp']):
                        syn = self.cell.insert_synapse(syntype, self.cell.dendlist[loc], 
                                                       pos, add_spine = add_spine, on_spine = on_spine)
                        self.add_input_generator(syn, syntype, deterministic = deterministic, numsyn = i, tstart = start)
                        if syntype == 'generalized_rule':
                            h.setpointer(h._ref_dopamine, 'dopamine', syn.obj)
                            h.setpointer(h._ref_stimulus_flag, 'stimulus_flag', syn.obj)
                        elif syntype == 'adaptive_glutamate_test':
                            h.setpointer(h._ref_dopamine, 'dopamine', syn.obj)
                    
                    elif syntype in ['spillover_test', 'my_spillover']:
                        syn1 = self.cell.insert_synapse('AMPA', self.cell.dendlist[loc], pos, 
                                                   add_spine = add_spine, on_spine = on_spine)
                        self.add_input_generator(syn1, 'AMPA', deterministic = deterministic, numsyn = i, tstart = start)
                        if self.cell.spines != []:                         
                            spines = [s for s in self.cell.spines if s.parent == self.cell.dendlist[loc]]
                            spines[i].syn_on = 0            
                        syn2 = self.cell.insert_synapse('NMDA', self.cell.dendlist[loc], pos, 
                                                   add_spine = 0, on_spine = 1)
                        self.connect_input_generator(syn2, 'NMDA', syn1.stim[-1])                               
                        syn3 = self.cell.insert_synapse('NMDAe', self.cell.dendlist[loc], pos, 
                                                               add_spine = 0, on_spine = 0)
                        syn1.clustered_flag = syn2.clustered_flag = syn3.clustered_flag = True
                        
                        if syntype == 'spillover_test':
                            self.connect_input_generator(syn3, 'NMDAe', syn2.stim[-1], delay = p.delay_exnmda)                                
                        elif syntype == 'my_spillover':
                            self.connect_input_generator(self.exglu[-1], 't_exglu', syn1.stim[-1])
                            self.connect_input_generator(syn3, 's_exglu', self.exglu[-1])
 
                    elif syntype in ['my_spillover_stp']:
                        syn1 = self.cell.insert_synapse('AMPA_stp', self.cell.dendlist[loc], pos, 
                                                   add_spine = add_spine, on_spine = on_spine)
                        self.add_input_generator(syn1, 'AMPA_stp', deterministic = deterministic, numsyn = i, tstart = start)
                        if self.cell.spines != []:                         
                            spines = [s for s in self.cell.spines if s.parent == self.cell.dendlist[loc]]
                            spines[i].syn_on = 0            
                        syn2 = self.cell.insert_synapse('NMDA_stp', self.cell.dendlist[loc], pos, 
                                                   add_spine = 0, on_spine = 1)
                        self.connect_input_generator(syn2, 'NMDA_stp', syn1.stim[-1])                               
                        syn3 = self.cell.insert_synapse('NMDAe', self.cell.dendlist[loc], pos, 
                                                               add_spine = 0, on_spine = 0)
                        syn1.clustered_flag = syn2.clustered_flag = syn3.clustered_flag = True
                        
                        self.connect_input_generator(self.exglu[-1], 't_exglu', syn1.stim[-1])
                        self.connect_input_generator(syn3, 's_exglu', self.exglu[-1])                            
                    
                    elif syntype == 'no_spillover':
                        syn1 = self.cell.insert_synapse('AMPA', self.cell.dendlist[loc], pos, 
                                                   add_spine = add_spine, on_spine = on_spine)
                        self.add_input_generator(syn1, 'AMPA', deterministic = deterministic, numsyn = i, tstart = start)
                        if self.cell.spines != []: 
                            spines = [s for s in self.cell.spines if s.parent == self.cell.dendlist[loc]]
                            spines[i].syn_on = 0
                        syn2 = self.cell.insert_synapse('NMDA', self.cell.dendlist[loc], pos, 
                                                   add_spine = 0, on_spine = 1)
                        self.connect_input_generator(syn2, 'NMDA', syn1.stim[-1])                                
                        syn1.clustered_flag = syn2.clustered_flag = True
                    
                    elif syntype == 'no_spillover_stp':
                        syn1 = self.cell.insert_synapse('AMPA_stp', self.cell.dendlist[loc], pos, 
                                                   add_spine = add_spine, on_spine = on_spine)
                        self.add_input_generator(syn1, 'AMPA_stp', deterministic = deterministic, numsyn = i, tstart = start)
                        if self.cell.spines != []: 
                            spines = [s for s in self.cell.spines if s.parent == self.cell.dendlist[loc]]
                            spines[i].syn_on = 0
                        syn2 = self.cell.insert_synapse('NMDA_stp', self.cell.dendlist[loc], pos, 
                                                   add_spine = 0, on_spine = 1)
                        self.connect_input_generator(syn2, 'NMDA_stp', syn1.stim[-1])                                
                        syn1.clustered_flag = syn2.clustered_flag = True
                
                if syntype == 'my_spillover':
                    self.set_exglu_weights()
                        
        elif syntype == 'pf':
            for loc in syn_loc:
                print(loc, p.input_dends.index(loc))
                syn_step = 1.0/num_syns
                cluster_start_pos = p.cluster_start_poss[p.input_dends.index(loc)]
                cluster_end_pos = p.cluster_end_poss[p.input_dends.index(loc)]
                for i in range(0, num_syns):
                    pos = cluster_start_pos + (cluster_end_pos - cluster_start_pos)*i*syn_step
                    start = rnd.uniform(p.plateau_burst_start, p.plateau_burst_end)
                    syn1 = self.cell.insert_synapse('AMPA_pf', self.cell.dendlist[loc],
                                                   pos, add_spine = add_spine, on_spine = on_spine)
                    self.add_input_generator(syn1, 'AMPA_pf', deterministic = deterministic, 
                                             numsyn = i, tstart = start )
                    syn2 = self.cell.insert_synapse('NMDA_pf', self.cell.dendlist[loc],
                                                   pos, add_spine = add_spine, on_spine = on_spine)
                    self.connect_input_generator(syn2, 'NMDA_pf', syn1.stim[-1])
                                           
#                elif (type(syn_loc) == dict):
#                    for ind, loc in enumerate(syn_loc['loc']):
#                        for i in range(0, p.plateau_cluster_size):
#                            syn = self.cell.insert_synapse('glutamate_ica_nmda', self.cell.dendlist[loc], 
#                                                           syn_loc['pos'][ind], add_spine = add_spine, on_spine = on_spine)
#                            self.add_input_generator(syn, 'glutamate_ica_nmda', 1.0 , tstart = syn_loc['start'][ind], 
#                                                     tend = syn_loc['end'][ind], deterministic = deterministic, numsyn = i)

    def add_input_generator(self, syn, syntype, freq_multiplier = 1, 
                            tstart = p.plateau_burst_start, tend = p.plateau_burst_end, deterministic = 0, numsyn = 1):
        
        if deterministic == 1:
            noise = 0
            start = tstart
            number = p.num_spikes
            interval = p.net_con_interval
            weight = p.gAMPAmax_plateau
            
            if syntype in ['input_syn']:
                weight = p.gAMPAmax
            elif syntype in ['NMDA', 'AMPA', 'NMDAe', 'AMPA_stp', 'NMDA_stp']:
                weight = p.weight
            elif syntype in ['inhexpsyn_plateau']:
                weight = p.gGABAmax_plateau
                number = 1
            
            
        elif deterministic == 0:
            noise = 1
            if syntype in ['ramp']:
                start = p.ramp_burst_start
                end = p.ramp_burst_end
                number = (end-start) * p.ramp_syn_rate
                interval = p.ramp_syn_interval
                weight = p.g_ramp_max

            elif syntype in ['inhexpsyn', 'inhexp2syn']:
                start = 0;
                end = p.simtime
                number = (end-start) * p.irate * freq_multiplier
                interval = p.i_interval/freq_multiplier
                weight = p.g_inhexpsyn_max
                
            elif syntype in ['inhexpsyn_plateau']:
                start = tstart
                end = p.inhibitory_burst_end
                number = p.num_spikes
                interval = 0
                weight = p.gGABAmax_plateau

            elif syntype in ['expsyn_plateau', 'plateau_cluster', 'nmda_plateau', 'tmGlut',
                             'glutamate_ica_nmda',
                            'glutamate_plateau', 
                            'adaptive_glutamate_hom', 
                            'glutamate_xor_test',
                            'AMPA', 'NMDA', 'NMDAe', 
                            'AMPA_stp', 'NMDA_stp', 
                            'adaptive_AMPA', 'adaptive_NMDA',
                            'adaptive_sAMPA', 'adaptive_sNMDA',
                            'adaptive_hom_AMPA', 'adaptive_hom_NMDA',
                            'AMPA_test', 'NMDA_test','adaptive_NMDAe',
                            'adaptive_shom_AMPA', 'adaptive_shom_NMDA', 'adaptive_my_shom_NMDA',
                            'adaptive_shom_AMPA_stp', 'adaptive_shom_NMDA_stp',
                            'adaptive_glutamate_shom',
                            'adaptive_cshom_AMPA', 'adaptive_cshom_NMDA',
                            'adaptive_glutamate_cshom', 'adaptive_sglutamate', 'NMDAe',
                            'adaptive_zahra_NMDA','adaptive_pf_AMPA', 'adaptive_pf_NMDA']:
                start = tstart; end = tend
                number = p.num_spikes#(end-start) * p.plateau_syn_rate
                interval = p.deterministic_interval#p.plateau_syn_interval              
                weight = 1.0
                if syntype in ['expsyn_plateau', 'tmGlut', 
                'glutamate_plateau']:
                    weight = p.gAMPAmax_plateau
                elif syntype in ['plateau_cluster', 'nmda_plateau', 'generalized_rule', 
                                 'generalized_rule_dist']:
                    weight = p.gNMDAmax_plateau
                elif syntype in ['nmda']:
                    weight = p.gNMDAmax
                elif syntype in ['NMDA', 'AMPA', 'NMDA_stp', 'AMPA_stp']:
                    weight = p.weight
                
            elif syntype in ['AMPA_pf', 'NMDA_pf', 'adaptive_pf_AMPA', 'adaptive_pf_NMDA']:
                start = p.pf_input_start
                end = p.pf_input_end
                number = p.pf_num_spikes
                interval = p.pf_input_interval
                weight = p.weight
            
            elif syntype in ['input_syn']:
                start = p.distributed_input_start
                end = p.distributed_input_end
                number = (end-start) * p.distributed_input_rate
                interval = p.distributed_input_interval
                weight = p.gAMPAmax
                
            elif syntype in ['adaptive_glutamate']:
                start = tstart; end = tend                
                number = (end-start) * p.low_rate
                interval = p.low_interval
                weight = 1.0

            elif syntype in ['expsyn', 'exp2syn', 'glutamate']:
                start = 0;
                end = p.simtime
                number = (end-start) * p.erate * freq_multiplier
                interval = p.e_interval/freq_multiplier
                weight = p.g_expsyn_max
                
        gen = h.NetStim(0.5, sec = self.presyn)
        gen.seed(int(time.time() + rnd.randint(1,10**7)))
        gen.start = start
        gen.noise = noise
        gen.number = number        
        gen.interval = interval
        
        nc = h.NetCon(gen, syn.obj)
        nc.delay = 0
        nc.weight[0] = weight
        
        if syntype in ['ramp']:
            self.ramp_estim.append(gen)            
            self.ramp_enc.append(nc)
            
        elif syntype in ['inhexpsyn', 'inhexp2syn']:
            syn.stim.append(gen) 
            syn.nc.append(nc)      

        elif syntype in ['inhexpsyn_plateau']:
            syn.stim.append(gen) 
            syn.nc.append(nc)      

        elif syntype in ['expsyn', 'exp2syn', 'plateau_cluster', 'input_syn', 
        'expsyn_plateau', 'nmda_plateau', 'tmGlut', 'glutamate',
        'adaptive_glutamate', 'glutamate_ica_nmda', 'glutamate_plateau',
        'pf',
        'AMPA', 'NMDA', 'NMDAe', 'adaptive_AMPA', 'NMDA_stp', 'AMPA_stp',
        'adaptive_NMDA', 'adaptive_hom_AMPA', 'adaptive_hom_NMDA', 
        'AMPA_test', 'NMDA_test','adaptive_shom_AMPA', 'adaptive_shom_NMDA',
        'adaptive_shom_AMPA_stp', 'adaptive_shom_NMDA_stp', 'adaptive_my_shom_NMDA',
        'adaptive_NMDAe', 'adaptive_glutamate_shom','adaptive_cshom_AMPA', 
        'adaptive_cshom_NMDA','adaptive_glutamate_cshom', 'adaptive_sAMPA', 'adaptive_sNMDA',
        'adaptive_sglutamate', 'adaptive_zahra_NMDA', 'AMPA_pf', 'NMDA_pf']:
            syn.stim.append(gen)
            syn.nc.append(nc)
            
    def set_up_recording(self, dend_record_list = [], record_step = p.record_step):
        self.dend_record_list = dend_record_list        
        self.vdlist = []
        self.t = h.Vector()
        self.t.record(h._ref_t, record_step)
        self.tv = h.Vector()
        self.tv.record(h._ref_t, p.record_step_v)
        self.tthresh = h.Vector()
        self.tthresh.record(h._ref_t, p.record_step_thresh)
        
        self.cali = []
        self.cali_dend = []
        self.cai_nmda = []
        self.cai = []  
        self.cati = []          
        self.cai_nmda_spine = []
        self.cai_spine = []            
        self.cali_spine = []            
        self.cati_spine = []          
        self.vspine = []
        self.kernel = []
        self.kernel_LTD = []
        self.cao = []
        
        self.vs = []
        self.vs = h.Vector()
        self.vs.record(self.cell.somalist[0](0.5)._ref_v, p.record_step_v)
        self.Cdur = []
                
        if self.exptype == 'record_ca':
            self.gk = []
            self.m = []
            self.ica_nmda = []
            
            self.cai_soma = h.Vector()
            self.cai_soma.record(self.cell.somalist[0](0.5)._ref_cai, record_step)
            for intfire in self.exglu:
                self.m.append(h.Vector())
                self.m[-1].record(intfire._ref_m, record_step)
            
            for d in dend_record_list:
                self.vdlist.append(h.Vector())
                self.vdlist[-1].record(self.cell.dendlist[d](p.pos)._ref_v, p.record_step_v)
                self.gk.append(h.Vector())
                self.gk[-1].record(self.cell.dendlist[d](p.pos)._ref_gk_kaf, record_step)                
            
                self.cai.append(h.Vector())
                self.cali.append(h.Vector())
#                self.cati.append(h.Vector())   
                self.cai_nmda.append(h.Vector())
                self.cali_dend.append(h.Vector())
                pos = p.cluster_start_poss[p.input_dends.index(d)]                
                self.cai[-1].record(self.cell.dendlist[d](pos)._ref_cai, record_step)
                self.cali[-1].record(self.cell.dendlist[d](pos)._ref_cali, record_step)
                self.cai_nmda[-1].record(self.cell.dendlist[d](pos)._ref_ca_nmdai, record_step)
                self.cali_dend[-1].record(self.cell.dendlist[d](pos)._ref_cali, record_step)
                    
                if self.cell.spines != []:
                    self.record_spinelist = [s for s in self.cell.spines if s.parent == self.cell.dendlist[d] and s.syn_on == 1]                                        
                    if p.include_empty_spines:
                        self.record_spinelist.extend([s for s in self.cell.spines if s.parent == self.cell.dendlist[d] and s.syn_on == 0])                                        
                    for spine in self.record_spinelist:
#                    spine = self.record_spinelist[0]                    
                        self.vspine.append(h.Vector())
                        self.vspine[-1].record(spine.head(0.5)._ref_v, p.record_step_v)
                        self.cali_spine.append(h.Vector())
                        self.cali_spine[-1].record(spine.head(0.5)._ref_cali, record_step)
                        self.cati_spine.append(h.Vector())   
                        self.cati_spine[-1].record(spine.head(0.5)._ref_cati, record_step)
                        self.cao.append(h.Vector())
                        self.cao[-1].record(spine.head(0.5)._ref_cao, record_step)
                        self.cai_spine.append(h.Vector())
                        self.cai_spine[-1].record(spine.head(0.5)._ref_cai, record_step)
                        self.cai_nmda_spine.append(h.Vector())
                        self.cai_nmda_spine[-1].record(spine.head(0.5)._ref_ca_nmdai, record_step)

#                        self.ical = []
#                        self.ical.append(h.Vector())
#                        self.ical[-1].record(spine.head(0.5)._ref_ical, record_step)
                        self.ica_nmda = []
                        self.ica_nmda.append(h.Vector())
                        self.ica_nmda[-1].record(spine.head(0.5)._ref_ica_nmda, record_step)

                        self.ica = []
                        self.ica.append(h.Vector())
                        self.ica[-1].record(spine.head(0.5)._ref_ica, record_step)
                        
        if self.exptype == 'record_i':
            self.ina = []
            self.ik = []
            self.iampa = []
            self.iNMDA = []
            self.ica_nmda = []

            for d in dend_record_list:
                self.vdlist.append(h.Vector())
                self.vdlist[-1].record(self.cell.dendlist[d](p.pos)._ref_v, p.record_step_v)
                self.ik.append(h.Vector())
                self.ik[-1].record(self.cell.dendlist[d](p.pos)._ref_ik, record_step)                
                self.ina.append(h.Vector())
                self.ina[-1].record(self.cell.dendlist[d](p.pos)._ref_ina, record_step)
                ampa_syns = [s for s in self.cell.esyn if s.sec == self.cell.dendlist[d] and s.type == 'AMPA']
                nmda_syns = [s for s in self.cell.esyn if s.sec == self.cell.dendlist[d] and s.type == 'NMDA']                
                for s in ampa_syns:
                    self.iampa.append(h.Vector())
                    self.iampa[-1].record(s.obj._ref_iAMPA, record_step)                
                for s in nmda_syns:
                    self.iNMDA.append(h.Vector())
                    self.iNMDA[-1].record(s.obj._ref_iNMDA, record_step)                
                    self.ica_nmda.append(h.Vector())
                    self.ica_nmda[-1].record(s.obj._ref_ica_nmda, record_step)                
            
        self.gsyn = []
        self.inmda = []
        self.A = []
        self.B = []

        self.w_ampa = []
        self.w_nmda = []
        self.lthresh_LTP = []
        self.hthresh_LTP = []
        self.lthresh_LTD = []
        self.cali_agh = []
        self.cati_agh = []
        self.cai_agh = []
        self.cai_nmda_agh = []
        self.v_agh = []
        self.kernel_agh = []
        self.kernel_LTD_agh = []
        self.delta_LTP = []      

    def plot_results(self):
        sns.set(font_scale = 1.5)
        sns.set_style("ticks")
                
        if self.exptype == 'record_ca':
            fig_vs = plt.figure()
            ax_vs = fig_vs.add_subplot(111)
            ax_vs.set_ylabel('Vs')
            ax_vs.set_xlabel('t')
            ax_vs.plot(self.tv, self.vs)
            
#            fig_cai_soma = plt.figure()
#            ax_cai_soma = fig_cai_soma.add_subplot(111)
#            ax_cai_soma.set_ylabel('soma cai')
#            ax_cai_soma.set_xlabel('t')
#            ax_cai_soma.plot(self.tout, self.cai_soma)
            
            fig_vd = plt.figure(); 
            ax_vd = fig_vd.add_subplot(111);
            ax_vd.set_ylabel('Vd')
            ax_vd.set_xlabel('t')
            for i in range(0,len(self.vdlist)):
                ax_vd.plot(self.tv, self.vdlist[i])
                
            fig_m = plt.figure(); 
            ax_m = fig_m.add_subplot(111);
            ax_m.set_ylabel('Intfire m')
            ax_m.set_xlabel('t')
            for i in range(0,len(self.m)):
                ax_m.plot(self.tout, self.m[i])
                    
                    
            legend_list = []
            for i in self.dend_record_list:
                legend_list.append('dend[%d](%.2f) = %.2f um' % (i,p.pos, h.distance(p.pos, sec = self.cell.dendlist[i]) ))
            
#            fig_gk = plt.figure()
#            ax_gk = fig_gk.add_subplot(111)
#            ax_gk.set_ylabel('gk')
#            ax_gk.set_xlabel('t')

#            for i in range(0,len(self.gk)):
#                ax_gk.plot(self.tout, self.gk[i])
            
#            fig_cai = plt.figure(); 
#            ax_cai = fig_cai.add_subplot(111);
#            ax_cai.set_ylabel('[Ca]_i')
#            ax_cai.set_xlabel('t')
#            for i in range(0,len(self.cai)):
#                ax_cai.plot(self.tout, self.cai[i])
#            
#            fig_cao = plt.figure(); 
#            ax_cao = fig_cao.add_subplot(111);
#            ax_cao.set_ylabel('[Ca]_o')
#            ax_cao.set_xlabel('t')
#            for i in range(0,len(self.cao)):
#                ax_cao.plot(self.tout, self.cao[i])
#            
#                    
#            fig_cali = plt.figure(); 
#            ax_cali = fig_cali.add_subplot(111);
#            ax_cali.set_ylabel('[Cal]_i')
#            ax_cali.set_xlabel('t')
#            for i in range(0,len(self.cali)):
#                ax_cali.plot(self.tout, self.cali[i])
#
#            fig_cai_nmda = plt.figure(); 
#            ax_cai_nmda = fig_cai_nmda.add_subplot(111);
#            ax_cai_nmda.set_ylabel('[Ca]_NMDA')
#            ax_cai_nmda.set_xlabel('t')
#            for i in range(0,len(self.cai_nmda)):
#                ax_cai_nmda.plot(self.tout, self.cai_nmda[i])


            if self.cell.spines != []:
                fig_vspine = plt.figure()
                ax_vspine = fig_vspine.add_subplot(111)
                ax_vspine.set_ylabel('Vspine')
                ax_vspine.set_xlabel('t')
                for v in self.vspine:
                    ax_vspine.plot(self.tv, v)

                fig_ica = plt.figure(); 
                ax_ica = fig_ica.add_subplot(111);
                ax_ica.set_ylabel('spine Ica')
                ax_ica.set_xlabel('t')
                for i in range(0,len(self.ica)):
                    ax_ica.plot(self.tout, self.ica[i])

#                fig_cai_nmda_spine = plt.figure(); 
#                ax_cai_nmda_spine = fig_cai_nmda_spine.add_subplot(111);
#                ax_cai_nmda_spine.set_ylabel('spine [Ca]_NMDA')
#                ax_cai_nmda_spine.set_xlabel('t')
#                for i in range(0,len(self.cai_nmda_spine)):
#                    ax_cai_nmda_spine.plot(self.tout, self.cai_nmda_spine[i])                    
#
#                fig_cai_spine = plt.figure(); 
#                ax_cai_spine = fig_cai_spine.add_subplot(111);
#                ax_cai_spine.set_ylabel('spine [Ca]')
#                ax_cai_spine.set_xlabel('t')
#                for i in range(0,len(self.cai_spine)):
#                    ax_cai_spine.plot(self.tout, self.cai_spine[i])
#
#                fig_cati_spine = plt.figure(); 
#                ax_cati_spine = fig_cati_spine.add_subplot(111);
#                ax_cati_spine.set_ylabel('spine [Cat]')
#                ax_cati_spine.set_xlabel('t')
#                for i in range(0,len(self.cai_spine)):
#                    ax_cati_spine.plot(self.tout, self.cati_spine[i])            
#
#                fig_cali_spine = plt.figure(); 
#                ax_cali_spine = fig_cali_spine.add_subplot(111);
#                ax_cali_spine.set_ylabel('spine [Cal]')
#                ax_cali_spine.set_xlabel('t')
#                for i in range(0,len(self.cali_spine)):
#                    ax_cali_spine.plot(self.tout, self.cali_spine[i])
#
#                fig_cali_dend = plt.figure()
#                ax_cali_dend = fig_cali_dend.add_subplot(111)
#                ax_cali_dend.set_ylabel('[Cal]_i dend')
#                ax_cali_dend.set_xlabel('t')
#                for c in self.cali_dend:
#                    ax_cali_dend.plot(self.tout, c)                          
            
            legend_list = []
            for i in self.dend_record_list:
                legend_list.append('dend[%d](%.2f) = %.2f um' % (i, p.pos, h.distance(p.pos, sec = self.cell.dendlist[i]) ))
            
#            ax_vd.legend(legend_list); ax_cai.legend(legend_list); 
#            ax_cali.legend(legend_list); ax_cai_nmda.legend(legend_list);
#            if self.cell.spines != []:
#                ax_cali_dend.legend(legend_list)

            plt.show()
#            return fig_vs, fig_vd, fig_cali, fig_cai_nmda

        if self.exptype in ['record_i']:
#            fig_vs = plt.figure()
#            ax_vs = fig_vs.add_subplot(111)
#            ax_vs.set_ylabel('Vs')
#            ax_vs.set_xlabel('t')
#            ax_vs.plot(self.tv, self.vs)
#            
#            fig_vd = plt.figure(); 
#            ax_vd = fig_vd.add_subplot(111);
#            ax_vd.set_ylabel('Vd')
#            ax_vd.set_xlabel('t')

#            for i in range(0,len(self.vdlist)):
#                ax_vd.plot(self.tv, self.vdlist[i])
#                    
#            legend_list = []
#            for i in self.dend_record_list:
#                legend_list.append('dend[%d](%.2f) = %.2f um' % (i,p.pos, h.distance(p.pos, sec = self.cell.dendlist[i]) ))
            
            fig_ik = plt.figure()
            ax_ik = fig_ik.add_subplot(111)
            ax_ik.set_ylabel('ik')
            ax_ik.set_xlabel('t')
            for i in range(0,len(self.ik)):
                ax_ik.plot(self.tout, self.ik[i])

            fig_ica_nmda = plt.figure(); 
            ax_ica_nmda = fig_ica_nmda.add_subplot(111);
            ax_ica_nmda.set_ylabel('ica_NMDA')
            ax_ica_nmda.set_xlabel('t')
            ax_ica_nmda.plot(self.tout, sum(self.ica_nmda))

            fig_ina = plt.figure()
            ax_ina = fig_ina.add_subplot(111)
            ax_ina.set_ylabel('ina')
            ax_ina.set_xlabel('t')
            for i in range(0,len(self.ina)):
                ax_ina.plot(self.tout, self.ina[i])

#            fig_iampa = plt.figure()
#            ax_iampa = fig_iampa.add_subplot(111)
#            ax_iampa.set_ylabel('iampa')
#            ax_iampa.set_xlabel('t')
#            ax_iampa.plot(self.tout, sum(self.iampa))

            fig_inmda = plt.figure()
            ax_inmda = fig_inmda.add_subplot(111)
            ax_inmda.set_ylabel('inmda')
            ax_inmda.set_xlabel('t')
            ax_inmda.plot(self.tout, sum(self.iNMDA))
            
            fig_i = plt.figure()
            ax_i = fig_i.add_subplot(111)
            ax_i.set_ylabel('i')
            ax_i.set_xlabel('t')
            for i in range(0,len(self.ik)):
                ax_i.plot(self.tout, self.ik[i]+self.ina[i]+sum(self.ica_nmda[i])+sum(self.iampa)+sum(self.iNMDA))

            plt.show()
            return fig_i, fig_ik, fig_inmda, fig_ica_nmda, fig_ina

    def set_up_experiment(self):
        if self.exptype in ['xor', 'xor_hom', 'xor_test_set', 'xor_gen', 'xor_spillover',
                            'xor_hom_spillover', 'xor_spillover_test','xor_shom_spillover',
                            'xor_cshom_spillover', 'xor_sspillover', 'xor_zahra_spillover',
                            'xor_shom_my_spillover', 'xor_shom_my_spillover_stp']:
            self.insert_synapses('noise_SPN')            
            self.create_dopamine()
            ts = self.create_training_set(p.training_set_size)
            print("Training set"); print(ts)
            self.training_set = ts; self.training_set_copy = self.training_set
            self.rewards_delivered = []
            self.xor_input_times, self. pf_input_times, self.reward_times = self.convert_ts_to_times(ts)
            self.create_xor_inputs(p.xor_input_size)
            self.connect_xor_inputs()
            self.spike_recorder(self.cell.somalist[0], 0.5, -10)
            
    def simulate(self, simtime = p.simtime, parallel = False):
        start = time.time()
        gmtime = time.gmtime(start)
        print("Starting simulation... %d:%d:%d" %(gmtime.tm_hour + 1, gmtime.tm_min, gmtime.tm_sec))
        
        if not parallel:        
            h.load_file("stdrun.hoc")
            h.init()
            h.tstop = simtime        
            fih3 = h.FInitializeHandler((self.seti_print_status))

            h.run()

        end = time.time()
        print("It took %.2f hours or %.4f seconds to simulate." % ((end-start)/3600 , (end-start)))

        if p.simtime > 1000:
                self.tout = np.multiply(np.asarray(self.t.to_python()), 0.001)
        else:
            self.tout = self.t.to_python()    

    def get_recording_points(self, dend, step):
        points = []        
        if dend.L < step:
            points.append(1.0)
        else:
            q = int(dend.L/step)
            for q in range(1,q+1):
                points.append(q*step/dend.L)
            points.append(1.0)
        return points
    
    def seti_stimulus_indicators(self, xor_input_times):
        times = []        
        for x in xor_input_times:
            times.extend(x)
        times.sort()
        times = (np.unique(times)).tolist()
        for t in times:
            h.cvode.event(t, (self.seti_stimulus_flag, 1))
    
    def seti_stimulus_flag(self, flag):
        h.stimulus_flag = flag
        if flag:
            h.cvode.event(h.t + p.training_input_length + p.time_to_reward + p.reward_length, (self.seti_stimulus_flag, 0))
        else:
            if self.exptype in ['xor_shom_my_spillover', 'xor_shom_my_spillover_stp']:
                self.set_exglu_weights()
                self.reset_exglu()
            
    def seti_print_status(self):
        update_points = np.arange(0, p.simtime, p.simtime/100. )
        for t in update_points:
            h.cvode.event(t, self.print_status)

    def reset_exglu(self):
        for s in self.exglu:
            s.m = 0
    
    def set_exglu_weights(self):
        if self.exptype == 'xor_shom_my_spillover':
            synlist = self.get_synapse_list('adaptive_my_shom_NMDA', clustered_flag = True)
            for s, exglu in zip(synlist, self.exnc):
                print(exglu.weight[0], s.obj.weight)
                exglu.weight[0] = s.obj.weight*p.exglu_norm_factor
                print(exglu.weight[0])
                
        elif self.exptype == 'xor_shom_my_spillover_stp':
            synlist = self.get_synapse_list('adaptive_shom_NMDA_stp', clustered_flag = True)
            for s, exglu in zip(synlist, self.exnc):
                exglu.weight[0] = s.obj.weight*p.exglu_norm_factor
        
        elif self.exptype == 'record_ca':
            for exglu in self.exnc:
                exglu.weight[0] = p.weight*p.exglu_norm_factor
        
        else:
            print("From method set_exglu_weights in %s" % type(self))
            print("Exptype '%s' not supported" % self.exptype)
            sys.exit(-1)
                
    def print_status(self):
        print("At time t %f, total simtime %f." % (h.t, p.simtime))

    def random_dend_list(self, input_dends, p):
        dend_list = []        
        for d in input_dends:
            if rnd.random() < p:
                dend_list.append(d)
        return dend_list
                    
    def get_synapse_list(self, syntype, clustered_flag = False):
        synlist = []
        for syn in self.cell.esyn:
            if syn.type == syntype and syn.clustered_flag == clustered_flag:
                synlist.append(syn)
        
        return synlist
        
    def delete_everything(self):
        self.dend_record_list = None        
        self.vdlist = None
        self.t = None
        self.tout = None
        self.short_tout = None
        self.glu =  None
        self.cali = None
        self.cali_dend = None
        self.cai_nmda = None
        self.cai = None            
        self.vspine = None
        self.vs = None
        self.record_spinelist = None
        self.vspine = None
        self.ical = None
        self.gsyn = None
        self.inmda = None
        self.A = None
        self.B = None
        self.w_ampa = None
        self.w_nmda = None
        self.lthresh_LTP = None
        self.hthresh_LTP = None
        self.lthresh_LTD = None
        self.soma_recorder_nc = None
        self.soma_recorder_tvec = None
        self.recorder_nc = None
        self.recorder_tvec = None
        self.ramp_enc = None
        self.ramp_estim = None
        
        self.distributed_input = None
        self.distributed_input_vectors = None
        self.distributed_input_nc = None
        self.distributed_input_list = None

        self.xor_input = None
        self.xor_input_times = None        
        self.xor_input_nc = None
        self.xor_input_vectors = None

        self.pf_input = None
        self.pf_input_vectors = None

        self.reward_times = None
        self.rewards_delivered = None
        self.cai_nmda_in_syns = None
        self.dopamine_vec = None
        
        self.enc = []
        self.estim = []
        self.cell.esyn = []
        self.inc = []
        self.istim = []
        self.cell.isyn = []
                
        self.g_nmda_pf = None
        self.cai_nmda_in_syns = None
        self.cali_in_syns = None
        self.presyn = h.Section()
        self.synlist = None
        
        self.training_set = None
        self.training_set_copy = None
        self.training_mode = None
        
    