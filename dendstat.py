# -*- coding: utf-8 -*-
"""
Created on Thu May 14 20:13:55 2020

@author: daniel
"""
import parameters as p

class DendStat(object):
    def __init__(self):
        self.dends = p.input_dends
        self.dend_inputs = [[] for d in self.dends ]
        self.pf_inputs = [[] for d in self.dends ]
        self.dend_syns = [[] for d in self.dends]
        self.dend_flags = [[0,0] for d in self.dends]
        self.dends00 = 0
        self.dends01 = 0
        self.dends10 = 0
        self.dends11 = 0
        self.syns_ampa = [[] for d in self.dends]
        self.syns_nmda = [[] for d in self.dends]
        self.syns_nmda_e = [[] for d in self.dends]
        self.count_01 = 0
        self.count_10 = 0
        self.count_11 = 0
        self.count_00 = 0
        self.distributed_inputs = []
        self.rewards_delivered = []
        self.training_set = []
        self.verror = []
        self.window_error = []
        self.vfull_error = []
        self.window_full_error = []
        self.syns_ampa_agh = [[],[], [], []]
        self.vdlist = []
        self.vs = []
    
    def dend_idx(self, dend):
        return self.dends.index(dend)