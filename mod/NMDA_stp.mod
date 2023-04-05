

TITLE simple NMDA receptors

COMMENT
-----------------------------------------------------------------------------

Essentially the same as /examples/nrniv/netcon/ampa.mod in the NEURON
distribution - i.e. Alain Destexhe's simple AMPA model - but with
different binding and unbinding rates and with a magnesium block.
Modified by Andrew Davison, The Babraham Institute, May 2000


	Simple model for glutamate AMPA receptors
	=========================================

  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS

    Whole-cell recorded postsynaptic currents mediated by AMPA/Kainate
    receptors (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994) were used
    to estimate the parameters of the present model; the fit was performed
    using a simplex algorithm (see Destexhe et al., J. Computational Neurosci.
    1: 195-230, 1994).

  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)

    The simplified model was obtained from a detailed synaptic model that
    included the release of transmitter in adjacent terminals, its lateral
    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe
    and Sejnowski, 1995).  Short pulses of transmitter with first-order
    kinetics were found to be the best fast alternative to represent the more
    detailed models.

  - ANALYTIC EXPRESSION

    The first-order model can be solved analytically, leading to a very fast
    mechanism for simulating synapses, since no differential equation must be
    solved (see references below).



References

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
   computing synaptic conductances based on a kinetic model of receptor binding
   Neural Computation 6: 10-14, 1994.

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a
   common kinetic formalism, Journal of Computational Neuroscience 1:
   195-230, 1994.

Orignal file by:
Kiki Sidiropoulou
Adjusted Cdur = 1 and Beta= 0.01 for better nmda spikes
PROCEDURE rate: FROM -140 TO 80 WITH 1000

Modified by Penny under the instruction of M.L.Hines on Oct 03, 2017
	Change gmax

-----------------------------------------------------------------------------
ENDCOMMENT



NEURON {
	POINT_PROCESS NMDA_stp
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA, ica_nmda, nmda_ca_fraction
	USEION ca_nmda READ ca_nmdai, ca_nmdao WRITE ica_nmda VALENCE 2
	NONSPECIFIC_CURRENT iNMDA
	RANGE mg, Cmax, eta, alpha, u0, U
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
	FARADAY = (faraday) (coulomb)
        R = (k-mole) (joule/degC)
}

PARAMETER {
	Cmax	= 1	 (mM)           : max transmitter concentration
	Cdur	= 1	 (ms)		: transmitter duration (rising phase)
	Alpha	= 4 (/ms /mM)	: forward (binding) rate (4)
	Beta 	= 0.01   (/ms)   : backward (unbinding) rate
	Erev	= 15	 (mV)		: reversal potential
    	mg   = 1      (mM)           : external magnesium concentration
    	eta = 0.28 (/mV)
    	alpha = 0.072 (/mV)
	gmax = 1   (uS)
	nmda_ca_fraction = 0.175 : (5% of the current is from Ca at 1mM extracellular Ca concentration) 
: previously nmda_ca_fraction = 0.7, from NMDA is 5-10 times more permeable to Ca++ than Na+ or K+, Ascher and Nowak, 1988; 

        tau = 3 (ms)
        tauR = 100 (ms)  : tauR > tau
        tauF = 0 (ms)  : tauF >= 0 (org: 800 ms)
        U = 0.3 (1) <0, 1>
        u0 = 0 (1) <0, 1>

}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iNMDA 		(nA)		: current = g*(v - e)
	g 		(uS)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
        B                       : magnesium block
	ica_nmda
	ca_nmdai
	ca_nmdao
	celsius (degC)
	x
}

STATE {Ron Roff}

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
        B = mgblock(v)
	g = (Ron + Roff)* gmax * B
	iNMDA = g*(v - Erev)
        :ica_nmda = g*nmda_ca_fraction*ghk(v, ca_nmdai, ca_nmdao)
        :iNMDA = iNMDA - ica_nmda
        ica_nmda = nmda_ca_fraction*iNMDA
        iNMDA = (1 - nmda_ca_fraction)*iNMDA
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
}

FUNCTION mgblock(v(mV)) {
        TABLE
        DEPEND mg
        FROM -140 TO 80 WITH 1000

        : from Jahr & Stevens

	 mgblock = 1 / (1 + mg * eta * exp(-alpha * v) )  :was 0.062, changed to 0.072 to get a better voltage-dependence of NMDA currents, july 2008, kiki

}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    if(z == 0) {
        z = z+1e-6
    }
    eco = co*(z)/(exp(z)-1)
    eci = ci*(-z)/(exp(-z)-1)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}

: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first


NET_RECEIVE(weight, on, nspike, r0, y, z, u, t0 (ms)) {

	INITIAL {
		y = 0
		z = 0
		u = u0
		t0 = t
	}
	z = z*exp(-(t-t0)/tauR)
	z = z + (y*(exp(-(t-t0)/tau) - exp(-(t-t0)/tauR)) / (tau/tauR - 1) )
	y = y*exp(-(t-t0)/tau)
	x = 1-y-z
	if (tauF > 0) {
		u = u*exp(-(t-t0)/tauF)
		u = u + U*(1-u)
	} else {
		u = U
	}


	: flag is an implicit argument of NET_RECEIVE and  normally 0
        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
		nspike = nspike + 1
		if (!on) {
			r0 = r0*exp(-Beta*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight*u
			state_discontinuity(Ron, Ron + r0)
			state_discontinuity(Roff, Roff - r0)
		}
:		 come again in Cdur with flag = current value of nspike
		net_send(Cdur, nspike)
       }
	if (flag == nspike) { : if this associated with last spike then turn off
		r0 = weight*Rinf*u + (r0 - weight*Rinf*u)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - weight*u
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
		on = 0
	}
    	
    	y = y + x*u
}
