TITLE simple AMPA receptors

COMMENT
-----------------------------------------------------------------------------

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

Modified by Penny under the instruction of M.L.Hines on Oct 03, 2017
	Change gmax

-----------------------------------------------------------------------------
ENDCOMMENT



NEURON {
	POINT_PROCESS adaptive_hom_AMPA
	RANGE R, gmax, g, ina, Alpha, Beta, iAMPA
	USEION na WRITE ina
	NONSPECIFIC_CURRENT  iAMPA
	RANGE Cdur, Erev, Rinf, Rtau
        POINTER dopamine, stimulus_flag
	RANGE thresh_LTP, thresh_LTD, w0, wmax, wmin, w_nmda
	RANGE learning_rate_w_LTP, learning_rate_w_LTD, thresh_LTP_max, thresh_LTP_min, thresh_LTP_0, learning_rate_thresh_LTP, thresh_LTD_max, thresh_LTD_min, thresh_LTD_0, learning_rate_thresh_LTD, LTD_thresh_factor
	RANGE ca_nmdai_max, cali_max, deriv, active_syn_flag, last_dopamine
        RANGE weight, thresh_LTP, thresh_LTD
	USEION cal READ cali VALENCE 2	
	USEION ca_nmda READ ca_nmdai

}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
    Cmax	= 0.1	(mM)		: max transmitter concentration
:	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
:	Alpha	= 0.94	(/ms)	: forward (binding) rate
	Alpha	= 1	(/ms)	: forward (binding) rate
:	Beta	= 0.018	(/ms)		: backward (unbinding) rate
	Beta	= 0.5 (/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		:0 reversal potential
	gmax    = 1  (uS)
    	
	learning_rate_w_LTP = 0.01
    	learning_rate_w_LTD = 0.01
    	wmax = 0.006 (uS)
    	wmin = 0.001 (uS)
        w0 = 0.00188 (uS)

	ca_nmdai_max = 0
	cali_max = 0
	active_syn_flag = 1e-6

	thresh_LTP_max = 0.5
	thresh_LTP_0 = 0.07
	thresh_LTP_min = 0.05
    	thresh_LTD_max = 0.05
	thresh_LTD_0 = 0.005
	thresh_LTD_min = 0.0005
        LTD_thresh_factor = 0.5

	learning_rate_thresh_LTP = 0.005
	learning_rate_thresh_LTD = 0.005
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iAMPA 		(nA)		: current = g*(v - Erev)
	g 		(uS)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
	ina
        dopamine
        last_dopamine
	stimulus_flag
        ca_nmdai        (mM)
        cali            (mM)

        weight
        thresh_LTP
        thresh_LTD
}

STATE {Ron Roff}

INITIAL {
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
    Rtau = 1 / ((Alpha * Cmax) + Beta)
    synon = 0
    weight = w0
    thresh_LTP = thresh_LTP_0
    thresh_LTD = thresh_LTD_0
    last_dopamine = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	g = (Ron + Roff)* gmax
	iAMPA = g*(v - Erev)
	ina = 0.9*iAMPA
	iAMPA = 0.1*iAMPA
	if (stimulus_flag == 1) {
        	ca_nmdai_max = max(ca_nmdai, ca_nmdai_max)
        	cali_max = max(cali, cali_max)
		last_dopamine = dopamine
        } else {
	  if (last_dopamine == 1 && active_syn_flag == 1) {

		  weight = weight + learning_rate_w_LTP * pind_LTP(ca_nmdai_max) * (wmax-weight)
		  thresh_LTP = thresh_LTP + learning_rate_thresh_LTP * pind_LTP(ca_nmdai_max)*(thresh_LTP_max - thresh_LTP)
		  thresh_LTD = thresh_LTD + learning_rate_thresh_LTD * pind_LTP(ca_nmdai_max)*(thresh_LTD_max - thresh_LTD)		  
          } else if (last_dopamine == -1 && active_syn_flag == 1) {

		  weight = weight - learning_rate_w_LTD * pind_LTD(cali_max) * (weight - wmin)
		  thresh_LTP = thresh_LTP - learning_rate_thresh_LTP * pind_LTD(cali_max)*(thresh_LTP - thresh_LTP_min)
		  thresh_LTD = thresh_LTD - learning_rate_thresh_LTD * pind_LTD(cali_max)*(thresh_LTD - thresh_LTD_min)

          }
          last_dopamine = dopamine		
          reset_max()
        }
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
		active_syn_flag = 1
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

FUNCTION pind_LTP(conc) {
    if (conc > thresh_LTP) {
	pind_LTP = 1
    } else {
	pind_LTP = 0
    }
}

FUNCTION pind_LTD(conc) {
    if (conc > thresh_LTD) {
	pind_LTD = 1
    } else {
	pind_LTD = 0
    }
}

FUNCTION max(current, maximum) {
   if (current>maximum) { 
      max = current
   } else {
      max = maximum
   }
}

PROCEDURE reset_max() {
	ca_nmdai_max = 0
        cali_max = 0
        active_syn_flag = 0
}
