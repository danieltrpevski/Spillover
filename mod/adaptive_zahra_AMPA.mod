TITLE simple AMPA receptors
NEURON {
	POINT_PROCESS adaptive_zahra_AMPA
	RANGE R, g, ina, Alpha, Beta, iAMPA
	USEION na WRITE ina
	NONSPECIFIC_CURRENT  iAMPA
    POINTER weight
    RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA,mg, Cmax, eta, alpha, treshf,synon,flag
	USEION ca_nmda READ ca_nmdai WRITE ica_nmda VALENCE 2	
	USEION cal READ cali VALENCE 2
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
    Cmax	= 0.1	(mM)		: max transmitter concentration
	Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
:	Alpha	= 0.94	(/ms)	: forward (binding) rate
	Alpha	= 1	(/ms)	: forward (binding) rate
:	Beta	= 0.018	(/ms)		: backward (unbinding) rate
	Beta	= 0.5 (/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		:0 reversal potential
	gmax    = 0 (uS)
  	  
	conc0=0 (mM)
    flagx=0 (1)
    mltype=0 (mM)
    eta=0   :change in the code
    rate_ltd=0 (1)
    rate_ltp=0 (1)
    treshf=0 (mM)    :change in the code
    treshltp=0.07(mM)
    tremin=0.07 (mM)
}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iAMPA 		(nA)		: current = g*(v - Erev)
	g 		(uS)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
	ina
    ca_nmdai        (mM)
    cali            (mM)
    ica_nmda        (nA)
    r0
    weight

}

STATE {Ron Roff  }

INITIAL {
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
    Rtau = 1 / ((Alpha * Cmax) + Beta)
    synon = 0
    

}

BREAKPOINT {
	SOLVE release METHOD cnexp
	g = (Ron + Roff)* gmax
	iAMPA = g*(v - Erev)
	ina = 0.9*iAMPA
	iAMPA = 0.1*iAMPA

	
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
        


}

: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first

NET_RECEIVE(dummy, on, nspike, r0, t0 (ms)) {
	: flag is an implicit argument of NET_RECEIVE and  normally 0
        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
		nspike = nspike + 1
		if (!on) {
			r0 = r0*exp(-Beta*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight
			state_discontinuity(Ron, Ron + r0)
			state_discontinuity(Roff, Roff - r0)
		}
		: come again in Cdur with flag = current value of nspike
		net_send(Cdur, nspike)
        }
	if (flag == nspike) { : if this associated with last spike then turn off
		r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - weight
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
		on = 0
	}

}



