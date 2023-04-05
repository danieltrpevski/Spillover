# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:18:11 2019

@author: daniel
"""
from neuron import h
import parameters as p

class Spine():
    """
    Spine class. Create a spine with neck and head.
    Based on Mattioni and Le Novere, (2013).
    https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=150284&file=/TimeScales-master/neuronControl/spine.py#tabs-2
    """
    
    def __init__(self, sec, name,         
                       neck_L = p.neck_L,      
                       neck_dia = p.neck_diam,
                       neck_Ra = p.neck_Ra,   
                       head_L = p.head_L,                          
                       head_dia = p.head_diam,    
                       head_Ra = p.head_Ra     ):
        """ Create a spine with geometry given by the arguments"""
        
        self.name         =   name
        self.neck       =   self.create_neck(neck_L, neck_dia, neck_Ra)
        self.head       =   self.create_head(self.neck, head_L, head_dia, head_Ra)
        self.parent     =   None # the parent section connected to the neck
        self.syn_on = 0 # Is there a synapse attached to the spine
        
    
    def create_neck(self, neck_L, neck_dia, Ra):
        """ Create a spine neck"""
        
        name_sec        =   self.name + "_neck"
        neck = h.Section(name = name_sec)
        neck.nseg       =   1
        neck.L          =   neck_L 
        neck.diam       =   neck_dia
        neck.Ra         =   Ra 
        neck.cm         =   1.0
        
        for mech in [   'pas',      \
                        'cav32',    \
                        'cav33',    \
                        'cadyn']:#,    \
                        #'cadyn_nmda']:
            neck.insert(mech)

        neck.g_pas      =   1.25e-5
        neck.e_pas      =   -85 
#        neck.taur_cadyn_nmda = p.tau_cadyn_nmda
        
        return neck
    
        
    def create_head(self, neck, head_L, head_dia, Ra):
        """Create the head of the spine and populate it with channels"""
        
        name_sec        =   self.name + "_head"
        head = h.Section(name = name_sec)        
        head.nseg       =   1
        head.L          =   head_L
        head.diam       =   head_dia
        head.Ra         =   Ra
        head.cm         =   1.0
        
        for mech in [   'pas',      \
                        'kir',      \
#                        'sk',      \
                        'cav32',    \
                        'cav33',    \
                        'car',      \
                        'cal12',    \
                        'cal13',    \
                        'cadyn',    \
                        'caldyn',   \
                        'catdyn',   \
                        'cadyn_nmda']:

            head.insert(mech)
        
        head.g_pas      =   1.25e-5
        head.e_pas      =   -85 
        head.depth_caldyn = 0.1
        head.depth_cadyn_nmda = 0.1
        head.taur_cadyn_nmda = p.tau_cadyn_nmda
        head.taur_caldyn = p.tau_caldyn
        head.taur_catdyn = p.tau_catdyn
#        head.gbar_sk = 0.5e-4
        head.connect(neck(1),0)
        
        return head

            
    def attach(self, parentSec, parentx, childx):
        """Attach a spine to a parentSec and store the parentSec into an attribute.
        Just an handy variation of the connect method"""
        self.neck.connect(parentSec, parentx, childx)
        self.parent = parentSec
        self.pos = parentx

#        self.head.pbar_cal12 = 	9.15e-7
#        self.head.pbar_cal13 = 	1.525e-7
#        self.head.pbar_cav32 = 1e-8
#        self.head.pbar_cav33 = 1e-8
#        self.head.pbar_car = 50e-6
                
#        self.head.pbar_cal12 = 	1.5e-6
#        self.head.pbar_cal13 = 	0.4e-6
#        self.head.pbar_cav32 = 1e-8
#        self.head.pbar_cav33 = 1e-8
#        self.head.pbar_car = 4e-5

        self.head.pbar_cal12 = 	3e-6
        self.head.pbar_cal13 = 	0.75e-6
        self.head.pbar_cav32 = 0.75e-7
        self.head.pbar_cav33 = 0.75e-7
        self.head.pbar_car = 4e-5

#        self.head.pbar_cal12 = 	self.parent(parentx).pbar_cal12
#        self.head.pbar_cal13 = 	self.parent(parentx).pbar_cal13
#        self.head.pbar_cav32 = self.parent(parentx).pbar_cav32
#        self.head.pbar_cav33 = self.parent(parentx).pbar_cav33
#        self.head.pbar_car = self.parent(parentx).pbar_car
#        self.head.gbar_kir = self.parent(parentx).gbar_kir
#        self.head.gbar_sk = self.parent(parentx).gbar_sk

        self.neck.pbar_cav32 = self.parent(parentx).pbar_cav32
        self.neck.pbar_cav33 = self.parent(parentx).pbar_cav33        