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
	POINT_PROCESS adaptive_cshom_AMPA
	RANGE R, gmax, g, ina, Alpha, Beta, iAMPA
	USEION na WRITE ina
	NONSPECIFIC_CURRENT  iAMPA
	RANGE Cdur, Erev, Rinf, Rtau
        RANGE weight, lthresh_LTP, lthresh_LTD, hthresh_LTP, last_dopamine, lthresh_LTP_min
        POINTER dopamine, stimulus_flag
	RANGE thresh_LTP, thresh_LTD, learning_rate, w0, wmax, wmin, steepness_LTP, steepness_LTD
	RANGE learning_rate_w_LTP, learning_rate_w_LTD, thresh_LTP_0, learning_rate_thresh_LTP, thresh_LTD_0, learning_rate_thresh_LTD 
	RANGE hthresh_LTP_const, hthresh_LTP_0, hthresh_max, n, delta, LTD_thresh_factor, width, delta_LTP
	RANGE ca_nmdai_max, cali_max, active_syn_flag, Cdur_init, Cdur_factor
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
:	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
:	Alpha	= 0.94	(/ms)	: forward (binding) rate
	Alpha	= 1	(/ms)	: forward (binding) rate
:	Beta	= 0.018	(/ms)		: backward (unbinding) rate
	Beta	= 0.5 (/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		:0 reversal potential
	gmax    = 1  (uS)
    	
	learning_rate = 0.01
	learning_rate_w_LTP = 0.01
	learning_rate_w_LTD = 0.01
    	wmax = 0.006 (uS)
    	wmin = 0.001 (uS)
        w0 = 0.00188 (uS)
	
	ca_nmdai_max = 0
	cali_max = 0
	active_syn_flag = 1e-6

	thresh_LTP_0 = 0.07
	thresh_LTD_0 = 0.005
	hthresh_LTP_0 = 0.5
	hthresh_max = 2.0
        lthresh_LTP_min = 0.055

	delta = 0.65
        width = 0.25
        hthresh_LTP_const = 0.05
	learning_rate_thresh_LTP = 0.005
	learning_rate_thresh_LTD = 0.005
	n = 4 : Hill coefficient
        LTD_thresh_factor = 1.0
	steepness_LTP = 0.25
	steepness_LTD = 2.5
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
        deriv
        ica_nmda        (nA)
	weight
	lthresh_LTP
	lthresh_LTD
	hthresh_LTP
	delta_LTP
}

STATE {Ron Roff}

INITIAL {
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
    Rtau = 1 / ((Alpha * Cmax) + Beta)
    synon = 0
    weight = w0
    lthresh_LTP = thresh_LTP_0
    lthresh_LTD = thresh_LTD_0
    hthresh_LTP = hthresh_LTP_0
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
                 
                 delta_LTP = sigmoidal(ca_nmdai_max, lthresh_LTP, steepness_LTP) * (1 - sigmoidal(ca_nmdai_max, lthresh_LTP, steepness_LTP))
		  weight = weight + learning_rate_w_LTP * delta_LTP
		  lthresh_LTP = lthresh_LTP + learning_rate_thresh_LTP * delta_LTP * (1 - 2*sigmoidal(ca_nmdai_max, lthresh_LTP, steepness_LTP))
          
          } else if (last_dopamine == -1 && active_syn_flag == 1) {

		  weight = weight - learning_rate_w_LTD * sigmoidal(cali_max, lthresh_LTD, steepness_LTD) * weight 
		  lthresh_LTP = lthresh_LTP - learning_rate_thresh_LTP * (lthresh_LTP - lthresh_LTP_min)

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

FUNCTION lthresh(conc, KD, steepness) {
:    lthresh = conc^n/(KD^n + conc^n) 
    lthresh = sigmoidal(conc,KD, steepness) 
}

FUNCTION hthresh(conc, KD, steepness) {
:    hthresh = KD^n/(KD^n + conc^n)
    hthresh = 1 - sigmoidal(conc, KD, steepness) 
}

FUNCTION max(current, maximum) {
   if (current>maximum) { 
      max = current
   } else {
      max = maximum
   }
}

FUNCTION min(current, minimum) {
   if (current<minimum) { 
      min = current
   } else {
      min = minimum
   }
}

PROCEDURE reset_max() {
	ca_nmdai_max = 0
        cali_max = 0
        active_syn_flag = 1e-6
}

FUNCTION sigmoidal(x, x_offset, s) {
    sigmoidal = 1/(1+exp(- s *(x - x_offset)))
}
