TITLE Very simple glutamate accumulation

UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
}

NEURON {
    POINT_PROCESS exglu
    RANGE glu
    POINTER stimulus_flag
}

ASSIGNED { 
    glu
    stimulus_flag
}

BREAKPOINT {
    if (stimulus_flag == 1) {
    	glu = 0
    }
}

INITIAL {  
    glu = 0
}

NET_RECEIVE( weight ) {
    glu = glu + weight
}


