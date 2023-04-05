# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 13:38:21 2017

@author: daniel
"""

from neuron import h
import parameters as p
import synapse as s
import spine as sp
import sys
import re
import random as rnd
import json

filename = p.params_for_learning_functions

class Neuron(object):
    """Generic neuron class."""
    def __init__(self):
        self.create_morphology()
        self.insert_channels()
        self.esyn = []
        self.isyn = []
        self.spines = []
        
    def create_morphology(self):
        """Create the cell morphology, including axial resistance and membrane capacitance."""
        raise NotImplementedError('"create_morphology() is not implemented."')
    
    def insert_channels(self):
        """Insert channels in the cell."""
        raise NotImplementedError("insert_channels() is not implemented.")
    
    def create_sectionlists(self):
        """Build subset lists. This defines 'all', but subclasses may
        want to define others. If overridden, call super() to include 'all'."""
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)
    
    def connect2target(self, target, thresh=10):
        """Make a new NetCon with this cell's membrane
        potential at the soma as the source (i.e. the spike detector)
        onto the target passed in (i.e. a synapse on a cell).
        Subclasses may override with other spike detectors."""
        nc = h.NetCon(self.soma(1)._ref_v, target, sec = self.soma)
        nc.threshold = thresh
        return nc
            
    def insert_synapse(self, syntype, sec, pos, add_spine = 0, on_spine = 0):
        if add_spine and on_spine:
            print("Arguments add_spine and on_spine can't simultaneously be 1")
            sys.exit(-1)

        if add_spine:
            s_ind = [int(si) for si in re.findall("\d+", sec.name())]
            s_ind = s_ind[0]
            self.num_spines_on_dends[s_ind] += 1

            spine_name = 'spine_' + sec.name() + '(' + str(pos) + ')'
            self.spines.append(sp.Spine(sec, spine_name))
            self.spines[-1].attach(sec, pos, 0)
            self.spines[-1].syn_on = 1
            sec = self.spines[-1].head
            s.spinepos = pos
            pos = 0.5
        
        if on_spine:
            empty_spines = [spine for spine in self.spines if (spine.parent == sec and spine.syn_on == 0)]
            
            if empty_spines == []:
                print("There are no empty spines on dendrite %s" % sec.name())
                sys.exit(-1)
            else:
                sec = empty_spines[0].head
                s.spinepos = pos
                pos = 0.5
                empty_spines[0].syn_on = 1
        
        syn = s.Synapse()
        syn.type = syntype
        syn.sec = sec
        syn.pos = pos        
        
        if syntype in ['expsyn', 'expsyn_plateau']:
            syn.obj = h.ExpSyn(sec(pos))            
            syn.obj.tau = p.esyn_tau
            syn.obj.e = p.e_esyn
            self.esyn.append(syn)
            return syn
        
        elif syntype == 'inhexpsyn' or syntype == 'inhexpsyn_plateau':
            syn.obj = h.InhExpSyn(sec(pos))
            if syntype == 'inhexpsyn':       
                syn.obj.tau = p.isyn_tau
            elif syntype == 'inhexpsyn_plateau':
                syn.obj.tau = p.isyn_plateau_tau
            syn.obj.e = p.e_gaba
            self.isyn.append(syn)
            return self.isyn[-1]
            
        elif syntype == 'exp2syn':
            syn.obj = h.Exp2Syn(sec(pos))
            syn.obj.e = p.e_esyn
            syn.obj.tau2 = p.tau2_exp2syn
            syn.obj.tau1 = p.tau1_exp2syn
            self.esyn.append(syn)
            return syn

        elif syntype == 'inhexp2syn':
            syn.obj = h.InhExp2Syn(sec(pos))
            syn.obj.e = p.e_gaba
            syn.obj.tau2 = p.tau2_inhexp2syn
            syn.obj.tau1 = p.tau1_inhexp2syn
            self.isyn.append(syn)
            return syn
            
        elif syntype == 'tmGlut':
            syn.obj = h.tmGlut(sec(pos))
            syn.obj.tau1_nmda = p.tau1_NMDA
            syn.obj.tau2_nmda = p.tau2_NMDA
            syn.obj.nmda_ratio = p.ratio_glutamate_syn
            syn.obj.U = 0.9
            syn.obj.tauF = 5.0
#            syn.obj.tau = 50
#            syn.obj.tauR = 100.0
            self.esyn.append(syn)
            return syn

        elif syntype == 'glutamate' or syntype == 'glutamate_plateau':
            syn.obj = h.glutamate(sec(pos))
            syn.obj.mg = p.Mg
            syn.obj.eta = p.eta
            syn.obj.alpha = p.alpha
            syn.obj.tau1_nmda = p.tau1_NMDA
            syn.obj.tau2_nmda = p.tau2_NMDA
            syn.obj.ratio = p.ratio_glutamate_syn
            self.esyn.append(syn)
            return syn

        elif syntype in ['glutamate_ica_nmda', 'glutamate_xor_test'] :
            syn.obj = h.glutamate_ica_nmda(sec(pos))
            syn.obj.mg = p.Mg            
            syn.obj.eta = p.eta
            syn.obj.alpha = p.alpha
            syn.obj.tau1_nmda = p.tau1_NMDA
            syn.obj.tau2_nmda = p.tau2_NMDA
            
            syn.obj.w_ampa = p.gAMPAmax_plateau 
            syn.obj.w_nmda = p.gNMDAmax_plateau
            
            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
            self.esyn.append(syn)
            return syn

        elif syntype in ['AMPA' ,'AMPA_test', 'AMPA_pf', 'AMPA_stp']:
            if syntype in ['AMPA', 'AMPA_pf']:
                syn.obj = h.AMPA(sec(pos))
            elif syntype in ['AMPA_stp']:
                syn.obj = h.AMPA_stp(sec(pos))
                syn.obj.U = p.U
                syn.obj.u0 = p.u0
            elif syntype == 'AMPA_test':
                syn.obj = h.AMPA_test(sec(pos))
                syn.obj.weight = p.weight
            syn.obj.gmax = p.gmaxAMPA_spillover
            if syntype == 'AMPA_pf':
                syn.obj.gmax = p.gmaxAMPA_pf
            self.esyn.append(syn)
            return syn         
        
        elif syntype in [ 'NMDA', 'NMDA_test', 'NMDAe', 'NMDA_pf', 'NMDA_stp']:
            if syntype in ['NMDA', 'NMDA_pf']:
                syn.obj = h.NMDA(sec(pos))
            elif syntype in ['NMDA_stp']:
                syn.obj = h.NMDA_stp(sec(pos))
                syn.obj.U = p.U
                syn.obj.u0 = p.u0
            elif syntype == 'NMDA_test':
                syn.obj = h.NMDA_test(sec(pos))
            elif syntype ==  'NMDAe':
                syn.obj = h.NMDAe(sec(pos))
            syn.obj.mg = p.Mg
            syn.obj.eta = p.eta
            syn.obj.alpha = p.alpha
            syn.obj.Cdur = p.Cdur
            if syntype in ['NMDA', 'NMDA_test', 'NMDA_stp']:    
                syn.obj.gmax = p.gmaxNMDA_spillover
            elif syntype in ['NMDA_pf']:    
                syn.obj.gmax = p.gmaxNMDA_pf
            elif syntype in ['NMDAe']:
                syn.obj.gmax = p.gmaxNMDAe_spillover
                syn.obj.Cdur_init = p.eCdur_init
                syn.obj.Cdur_factor = p.eCdur_factor
                syn.obj.weight = p.exglu_weight
            syn.obj.Beta = p.Beta
            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
            if syntype == 'NMDA_test':
                syn.obj.weight = p.weight
            
            self.esyn.append(syn)
            return syn        

        elif syntype in ['adaptive_shom_AMPA', 'adaptive_shom_AMPA_stp']:
            if syntype in ['adaptive_shom_AMPA']:
                syn.obj = h.adaptive_shom_AMPA(sec(pos))
            else:
                syn.obj = h.adaptive_shom_AMPA_stp(sec(pos))
                syn.obj.U = p.U
                syn.obj.u0 = p.u0
                
            syn.obj.gmax = p.gmaxAMPA_spillover
            
            self.esyn.append(syn)
            return syn        

        elif syntype == 'adaptive_pf_AMPA':
            syn.obj = h.adaptive_pf_AMPA(sec(pos))
            syn.obj.gmax = p.gmaxAMPA_pf
            
            self.esyn.append(syn)
            return syn               
        
        elif syntype in ['adaptive_shom_NMDA','adaptive_shom_NMDA_stp','adaptive_my_shom_NMDA']:
            if syntype in ['adaptive_shom_NMDA']:
                syn.obj = h.adaptive_shom_NMDA(sec(pos))
            elif syntype in ['adaptive_shom_NMDA_stp']:
                syn.obj = h.adaptive_shom_NMDA_stp(sec(pos))
                syn.obj.U = p.U
                syn.obj.u0 = p.u0
            else:
                syn.obj = h.adaptive_my_shom_NMDA(sec(pos))
            
            syn.obj.mg = p.Mg
            syn.obj.eta = p.eta
            syn.obj.alpha = p.alpha
            syn.obj.gmax = p.gmaxNMDA_spillover
            syn.obj.Beta = p.Beta
            syn.obj.Cdur = p.Cdur
            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
            
            syn.obj.w0 = p.weight            
            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD
            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD
            syn.obj.learning_rate_thresh_KD_LTD = p.learning_rate_thresh_KD_LTD
            syn.obj.KD1 = p.KD1
            syn.obj.n1 = p.n1
            syn.obj.KD2 = p.KD2
            syn.obj.n2 = p.n2
            syn.obj.KD_LTD = p.KD_LTD
            syn.obj.n_LTD = p.n_LTD    

            self.esyn.append(syn)
            return syn  

        elif syntype == 'adaptive_pf_NMDA':
            syn.obj = h.adaptive_pf_NMDA(sec(pos))
            
            syn.obj.mg = p.Mg
            syn.obj.eta = p.eta
            syn.obj.alpha = p.alpha
            syn.obj.gmax = p.gmaxNMDA_pf
            syn.obj.Beta = p.Beta
            syn.obj.Cdur = p.Cdur_pf
            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
            
            syn.obj.w0 = p.weight            
            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD_pf
            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD
            syn.obj.KD_LTD = p.KD_LTD_pf
            syn.obj.n_LTD = p.n_LTD_pf
            
            self.esyn.append(syn)
            return syn  

        elif syntype == 'adaptive_NMDAe':
            syn.obj = h.adaptive_NMDAe(sec(pos))
            syn.obj.mg = p.Mg
            syn.obj.eta = p.eta
            syn.obj.alpha = p.alpha
            syn.obj.Erev = p.erev_NMDA
            syn.obj.gmax = p.gmaxNMDAe_spillover
            syn.obj.Beta = p.Beta
            syn.obj.Cdur = p.eCdur
            syn.obj.Cdur_init = p.eCdur_init
            syn.obj.Cdur_factor = p.eCdur_factor
            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
            
            self.esyn.append(syn)
            return syn           

#        elif syntype == 'adaptive_glutamate':            
#            syn.obj = h.adaptive_glutamate(sec(pos))
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.w0 = p.gAMPAmax_plateau
#            syn.obj.wmax = p.gAMPAmax_plateau*p.LTP_factor            
#            syn.obj.wmin = p.gAMPAmax_plateau*p.LTD_factor
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD
#            syn.obj.thresh_LTD = p.thresh_LTD
#            syn.obj.thresh_LTP = p.thresh_LTP
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            syn.obj.NMDA_AMPA_ratio == p.ratio_distributed_synapses
#
#            syn.obj.tau1_nmda = p.tau1_NMDA
#            syn.obj.tau2_nmda = p.tau2_NMDA
#            
#            #            if add_spine:
##                sec = self.spines[-1].parent
##                pos = self.spines[-1].pos
##            h.setpointer(sec(pos)._ref_cali, 'cali', syn.obj)
#
#            self.esyn.append(syn)
#            return syn

#        elif syntype == 'adaptive_sglutamate':            
#            syn.obj = h.adaptive_sglutamate(sec(pos))
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.w0 = p.gAMPAmax_plateau
#            syn.obj.NMDA_AMPA_ratio = p.ratio_distributed_synapses            
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            
#            syn.obj.tau1_nmda = p.tau1_NMDA
#            syn.obj.tau2_nmda = p.tau2_NMDA
#            
#            with open(filename, 'r') as f:
#                to_read = json.load(f)
#            res_dict = json.loads(to_read)
#            params_LTP = res_dict['cai_nmda_params_by_dend']
#            params_LTD = res_dict['cali_params_by_dend']
#            r = re.findall("\[\d+\]", sec.name())
#            r = [int(num) for elem in r for num in re.findall("\d+", elem)]
#            KD1, n1, KD2, n2, factor_LTP = params_LTP[0]
#            KD_LTD, n_LTD, factor_LTD = params_LTD[0]
#            
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP*0.25
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD*100
#            syn.obj.n1 = p.n1
#            syn.obj.KD2 = p.KD2
#            syn.obj.n2 = p.n2
#            syn.obj.KD_LTD = p.KD_LTD
#            syn.obj.n_LTD = p.n_LTD   
#            
#            self.esyn.append(syn)
#            return syn

#        elif syntype in ['adaptive_glutamate_hom']:            
#            if p.random_weights == True:
#                w_ampa = rnd.uniform(p.gAMPAmax_plateau*p.LTD_factor, p.gAMPAmax_plateau*p.LTP_factor)
#            else:
#                w_ampa = p.gAMPAmax_plateau
#            syn.obj = h.adaptive_glutamate_hom(sec(pos))
#
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            syn.obj.w0 = w_ampa
#            syn.obj.wmax = p.gAMPAmax_plateau*p.LTP_factor_di            
#            syn.obj.wmin = p.gAMPAmax_plateau*p.LTD_factor_di
#            syn.obj.NMDA_AMPA_ratio = p.ratio_distributed_synapses            
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD            
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.thresh_LTP_min = p.thresh_LTP_min                  
#            syn.obj.thresh_LTD_min = p.thresh_LTD_min                  
#            
#            syn.obj.thresh_LTP_max = p.thresh_LTP_max            
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP_di
#            syn.obj.thresh_LTD_max = p.thresh_LTD_max            
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD_di
#
#            syn.obj.tau1_nmda = p.tau1_NMDA
#            syn.obj.tau2_nmda = p.tau2_NMDA
#            
#            self.esyn.append(syn)
#            return syn
#
#        elif syntype in ['adaptive_glutamate_shom']:            
#            if p.random_weights == True:
#                weight = rnd.uniform(p.gAMPAmax_plateau*p.LTD_factor, p.gAMPAmax_plateau*p.LTP_factor)
#            else:
#                weight = p.gAMPAmax_plateau
#            syn.obj = h.adaptive_glutamate_shom(sec(pos))
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            syn.obj.w0 = weight
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD            
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.hthresh_LTP_0 = p.hthresh_LTP                  
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD
#            syn.obj.NMDA_AMPA_ratio = p.ratio_distributed_synapses
#            syn.obj.tau1_nmda = p.tau1_NMDA
#            syn.obj.tau2_nmda = p.tau2_NMDA            
##            if add_spine:
##                sec = self.spines[-1].parent
##                pos = self.spines[-1].pos
##            h.setpointer(sec(pos)._ref_cali, 'cali', syn.obj)
#
#            self.esyn.append(syn)
#            return syn
#            
#        elif syntype in ['adaptive_glutamate_cshom']:            
#            if p.random_weights == True:
#                weight = rnd.uniform(p.gAMPAmax_plateau*p.LTD_factor, p.gAMPAmax_plateau*p.LTP_factor)
#            else:
#                weight = p.gAMPAmax_plateau
#            syn.obj = h.adaptive_glutamate_cshom(sec(pos))
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            syn.obj.w0 = weight
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD            
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.hthresh_LTP_0 = p.hthresh_LTP
#            syn.obj.lthresh_LTP_min = p.lthresh_LTP_min                  
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP*2.5
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD*2.5
#            syn.obj.NMDA_AMPA_ratio = p.ratio_distributed_synapses
#            syn.obj.tau1_nmda = p.tau1_NMDA
#            syn.obj.tau2_nmda = p.tau2_NMDA
#            syn.obj.steepness_LTP = p.steepness_LTP
#            syn.obj.steepness_LTD = p.steepness_LTD
            #            if add_spine:
#                sec = self.spines[-1].parent
#                pos = self.spines[-1].pos
#            h.setpointer(sec(pos)._ref_cali, 'cali', syn.obj)

#            self.esyn.append(syn)
#            return syn   

#        elif syntype == 'adaptive_AMPA':
#            syn.obj = h.adaptive_AMPA(sec(pos))
#            syn.obj.gmax = p.gmaxAMPA_spillover
#
#            syn.obj.w0 = p.weight
#            syn.obj.wmax = p.weight*p.LTP_factor            
#            syn.obj.wmin = p.weight*p.LTD_factor
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w
#            syn.obj.thresh_LTD = p.thresh_LTD
#            syn.obj.thresh_LTP = p.thresh_LTP
#            
#            self.esyn.append(syn)
#            return syn        

#        elif syntype == 'adaptive_sAMPA':
#            syn.obj = h.adaptive_sAMPA(sec(pos))
#            syn.obj.gmax = p.gmaxAMPA_spillover
#                  
#            self.esyn.append(syn)
#            return syn        
#
#        elif syntype == 'adaptive_NMDA':
#            syn.obj = h.adaptive_NMDA(sec(pos))
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.Erev = p.erev_NMDA
#            syn.obj.gmax = p.gmaxNMDA_spillover
#            syn.obj.Beta = p.Beta
#            syn.obj.Cdur = p.eCdur
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            
#            syn.obj.w0 = p.weight
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w
#            syn.obj.thresh_LTD = p.thresh_LTD
#            syn.obj.thresh_LTP = p.thresh_LTP
#
#            syn.obj.Cdur_init = p.eCdur_init
#            syn.obj.Cdur_factor = p.eCdur_factor            
#
#            self.esyn.append(syn)
#            return syn        
            
#        elif syntype == 'adaptive_sNMDA':
#            syn.obj = h.adaptive_sNMDA(sec(pos))
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.Erev = p.erev_NMDA
#            syn.obj.gmax = p.gmaxNMDA_spillover
#            syn.obj.Beta = p.Beta
#            syn.obj.Cdur = p.Cdur
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            
#            syn.obj.w0 = p.weight

#            with open(filename, 'r') as f:
#                to_read = json.load(f)
#            res_dict = json.loads(to_read)
#            params_LTP = res_dict['cai_nmda_params_by_dend']
#            params_LTD = res_dict['cali_params_by_dend']
#            r = re.findall("\[\d+\]", sec.name())
#            r = [int(num) for elem in r for num in re.findall("\d+", elem)]
#            print(r[0])
#            KD1, n1, KD2, n2, factor_LTP = params_LTP[p.independent_dends.index(r[0])]
#            KD_LTD, n_LTD, factor_LTD = params_LTD[p.independent_dends.index(r[0])]
            
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD
#            syn.obj.KD1 = p.KD1
#            syn.obj.n1 = p.n1
#            syn.obj.KD2 = p.KD2
#            syn.obj.n2 = p.n2
#            syn.obj.KD_LTD = p.KD_LTD
#            syn.obj.n_LTD = p.n_LTD           
#
#            self.esyn.append(syn)
#            return syn        

#        elif syntype == 'adaptive_hom_AMPA':
#            syn.obj = h.adaptive_hom_AMPA(sec(pos))
#            syn.obj.gmax = p.gmaxAMPA_spillover
#            
#            syn.obj.w0 = p.weight
#            syn.obj.wmax = p.weight*p.LTP_factor            
#            syn.obj.wmin = p.weight*p.LTD_factor
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD            
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.thresh_LTD_min = p.thresh_LTD_min
#            syn.obj.thresh_LTP_min = p.thresh_LTP_min
#            syn.obj.thresh_LTD_max = p.thresh_LTD_max
#            syn.obj.thresh_LTP_max = p.thresh_LTP_max
#            syn.obj.LTD_thresh_factor = p.LTD_thresh_factor            
#            
#            self.esyn.append(syn)
#            return syn        

 

#        elif syntype == 'adaptive_zahra_AMPA':
#            syn.obj = h.adaptive_zahra_AMPA(sec(pos))
#            syn.obj.gmax = p.gmaxAMPA_spillover
#            
#            self.esyn.append(syn)
#            return syn        
#
#        elif syntype == 'adaptive_cshom_AMPA':
#            syn.obj = h.adaptive_cshom_AMPA(sec(pos))
#            syn.obj.gmax = p.gmaxAMPA_spillover
#
#            syn.obj.width = p.width
#            syn.obj.w0 = p.weight
#            syn.obj.wmax = p.weight*p.LTP_factor            
#            syn.obj.wmin = p.weight*p.LTD_factor
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.hthresh_LTP_0 = p.hthresh_LTP
#            syn.obj.hthresh_LTP_const = p.hthresh_LTP_const
#            syn.obj.n = p.Hill_coefficient
#            syn.obj.LTD_thresh_factor = p.LTD_thresh_factor
#            syn.obj.lthresh_LTP_min = p.lthresh_LTP_min
#            
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD
#            syn.obj.steepness_LTP = p.steepness_LTP
#            syn.obj.steepness_LTD = p.steepness_LTD
#            self.esyn.append(syn)
#            return syn        
        
#        elif syntype in ['adaptive_hom_NMDA']:
#            syn.obj = h.adaptive_hom_NMDA(sec(pos))
#
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.gmax = p.gmaxNMDA_spillover
#            syn.obj.Beta = p.Beta
#            syn.obj.Cdur = p.eCdur
#            syn.obj.Cdur_init = p.eCdur_init
#            syn.obj.Cdur_factor = p.eCdur_factor
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            
#            syn.obj.w0 = p.weight
#            syn.obj.wmax = p.weight*p.LTP_factor            
#            syn.obj.wmin = p.weight*p.LTD_factor
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.thresh_LTD_min = p.thresh_LTD_min
#            syn.obj.thresh_LTP_min = p.thresh_LTP_min
#            
#            syn.obj.thresh_LTD_max = p.thresh_LTD_max
#            syn.obj.thresh_LTP_max = p.thresh_LTP_max
#            syn.obj.LTD_thresh_factor = p.LTD_thresh_factor
#
#            self.esyn.append(syn)
#            return syn        

#        elif syntype == 'adaptive_zahra_NMDA':
#            syn.obj = h.adaptive_zahra_NMDA(sec(pos))
#            
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.gmax = p.gmaxNMDA_spillover
#            syn.obj.Beta = p.Beta
#            syn.obj.Cdur = p.Cdur
#            
#            syn.obj.w0 = p.weight            
#            syn.obj.rate_ltp = 0.00004
#            syn.obj.rate_ltd = 0.0007
#            syn.obj.rate_ltp_tresh = 0.000006
#            syn.obj.rate_ltd_thrsh = 0.00005
#            syn.obj.tremin = 0.03
#            self.esyn.append(syn)
#            return syn  
#        
#        elif syntype == 'adaptive_cshom_NMDA':
#            syn.obj = h.adaptive_cshom_NMDA(sec(pos))
#            
#            syn.obj.mg = p.Mg
#            syn.obj.eta = p.eta
#            syn.obj.alpha = p.alpha
#            syn.obj.gmax = p.gmaxNMDA_spillover
#            syn.obj.Beta = p.Beta
#            syn.obj.Cdur = p.Cdur
#            syn.obj.n = p.Hill_coefficient
#            syn.obj.nmda_ca_fraction = p.nmda_ca_fraction
#            
#            syn.obj.width = p.width
#            syn.obj.w0 = p.weight
#            syn.obj.learning_rate_w_LTP = p.learning_rate_w_LTP
#            syn.obj.learning_rate_w_LTD = p.learning_rate_w_LTD
#            syn.obj.thresh_LTD_0 = p.thresh_LTD
#            syn.obj.thresh_LTP_0 = p.thresh_LTP
#            syn.obj.hthresh_LTP_0 = p.hthresh_LTP
#            syn.obj.hthresh_LTP_const = p.hthresh_LTP_const
#            syn.obj.LTD_thresh_factor = p.LTD_thresh_factor
#            syn.obj.lthresh_LTP_min = p.lthresh_LTP_min
#            
#            syn.obj.learning_rate_thresh_LTP = p.learning_rate_thresh_LTP
#            syn.obj.learning_rate_thresh_LTD = p.learning_rate_thresh_LTD
#            syn.obj.steepness_LTP = p.steepness_LTP
#            syn.obj.steepness_LTD = p.steepness_LTD
#            self.esyn.append(syn)
#            return syn        
        
        else: 
            print("From method cell.insert_synapse")
            print("Syntype '%s' not supported" % syntype)
            sys.exit(-1)

    def max_dist(self, axon_excluding=True):
        if not hasattr(self, 'somalist'):
            raise NotImplementedError("create_sectionlists() is not implemented or attribute somalist not defined")
        
        h.distance(sec=self.somalist[0])
        dmax = 0
        for sec in self.all:
            if axon_excluding and sec.name().find('axon') == 0: 
                continue
            dmax = max(dmax, h.distance(1, sec=sec))
        return dmax
        
    def get_nsegs(self):
        """Returns the number of segments in the neuron model."""
        nsegs = 0
        for sec in self.all: 
            nsegs += sec.nseg
        return nsegs
        
    def set_nsegs(self):
        """Sets the number of segments in each section of the neuron model
        according to n = 2*int(L/40) + 1, where L is the length of the section."""
        for sec in self.all:
            sec.nseg = 2*int(sec.L/40.0)+1
        if hasattr(self, 'axonlist'):
            for sec in self.axonlist:
                sec.nseg = 2  # two segments in axon initial segment

    def total_dend_length(self):
        """Returns the total dendritic length."""
        total_length = 0             
        for dend in self.dendlist:
            total_length += dend.L
        return total_length
        
    def increase_dend_res(self, dend_list, mult):
        for d in dend_list:
            self.dendlist[d].nseg *= mult
            