# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 18:55:32 2019

@author: daniel
"""

class Synapse(object):
    """Generic synapse class."""
    def __init__(self):
        self.type = None
        self.sec = None
        self.pos = None
        self.spinepos = None
        self.obj = None
        self.source = None
        self.ref_var_ampa = None
        self.ref_var_nmda = None
        self.ref_var_lthresh_LTP = None
        self.ref_var_hthresh_LTP = None
        self.ref_var_lthresh_LTD = None
        self.ref_var_cai_nmda = None
        self.ref_var_cali = None
        self.erec = None
        self.clustered_flag = False
        self.stim = []
        self.nc = []