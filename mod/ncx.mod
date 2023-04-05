TITLE Calcium dynamics for NMDA calcium pool

NEURON {
    SUFFIX ncx
    USEION ca_nmda READ ica_nmda, ca_nmdai WRITE ca_nmdai VALENCE 2
    RANGE cainf, taur, kt, kd, drive
}

UNITS {
    (molar) = (1/liter) 
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    (msM) = (ms mM)
    FARADAY = (faraday) (coulomb)
}

PARAMETER {
    cainf = 50e-6 (mM)
    taur = 100 (ms)
    kt = 10 (mM/ms)
    kd = 1e-3 (mM)
    drive = 0.02
}

STATE { ca_nmdai (mM) }

INITIAL { ca_nmdai = cainf }

ASSIGNED {
    ica_nmda (mA/cm2)
}
    
BREAKPOINT {
    SOLVE state METHOD cnexp
}

DERIVATIVE state { 
    ca_nmdai' = -drive * kt*ca_nmdai/(ca_nmdai+kd)
}

COMMENT

Original model by Wolf (2005) and Destexhe (1992).

Ca shell parameters by Evans (2012), with kb but without pump.

NEURON implementation by Alexander Kozlov <akozlov@nada.kth.se>.

ENDCOMMENT
