

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
	POINT_PROCESS adaptive_thalamus
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA
	NONSPECIFIC_CURRENT  iNMDA
	RANGE mg, Cmax, eta, alpha, nmda_ca_fraction
        POINTER dopamine, stimulus_flag
	RANGE w0, weight, lthresh_LTP, lthresh_LTD, hthresh_LTP, last_dopamine, width, lthresh_LTP_min
	RANGE learning_rate_w_LTP, learning_rate_w_LTD 
        RANGE thresh_LTP_0, learning_rate_thresh_LTP 
        RANGE thresh_LTD_0, learning_rate_thresh_LTD
	RANGE hthresh_LTP_const, hthresh_LTP_0, hthresh_max, n, delta, LTD_thresh_factor
	RANGE ca_nmdai_max, cali_max, active_syn_flag, Cdur_init, Cdur_factor, steepness
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
	Cmax	= 1	 (mM)           : max transmitter concentration
	Alpha	= 4 (/ms /mM)	: forward (binding) rate (4)
	Beta 	= 0.01   (/ms)   : backward (unbinding) rate
	Erev	= 0	 (mV)		: reversal potential
        mg   = 1      (mM)           : external magnesium concentration
        eta = 0.28 (/mV)
        alpha = 0.072 (/mV)
	gmax = 1   (uS)

	Cdur_init = 50
	Cdur_factor = 200

	learning_rate_w_LTP = 0.01
	learning_rate_w_LTD = 0.01
    	w0 = 0.00188 (uS)
    	
	ca_nmdai_max = 0
	cali_max = 0
	active_syn_flag = 1e-6

	thresh_LTP_0 = 0.07
	thresh_LTD_0 = 0.005
	hthresh_LTP_0 = 0.5
	hthresh_max = 2.0
        lthresh_LTP_min = 0.06

	delta = 0.65
        width = 0.25
        hthresh_LTP_const = 0.05
	learning_rate_thresh_LTP = 0.005
	learning_rate_thresh_LTD = 0.005
	n = 4 : Hill coefficient
        LTD_thresh_factor = 1.0
	steepness = 50
	nmda_ca_fraction = 0.15
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iNMDA 		(nA)		: current = g*(v - e)
	g 		(uS)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
        B                       : magnesium block
	ica_nmda        (nA)
	dopamine
	last_dopamine
        stimulus_flag
        cali            (mM)
	ca_nmdai        (mM)

	Cdur (ms)		: transmitter duration (rising phase)
	weight
	lthresh_LTP
	lthresh_LTD
	hthresh_LTP
}

STATE {Ron Roff}

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
	weight = w0
        lthresh_LTP = thresh_LTP_0
        lthresh_LTD = thresh_LTD_0
	hthresh_LTP = lthresh_LTP + width
	Cdur	= 1	: transmitter duration (rising phase)
	last_dopamine = 0
}

BREAKPOINT {
        SOLVE release METHOD cnexp
        B = mgblock(v)
	g = (Ron + Roff)* gmax * B
	iNMDA = g*(v - Erev)
        ica_nmda = nmda_ca_fraction*iNMDA
        iNMDA = (1 - nmda_ca_fraction)*iNMDA

	if (stimulus_flag == 1) {
        	ca_nmdai_max = max(ca_nmdai, ca_nmdai_max)
        	cali_max = max(cali, cali_max)
		last_dopamine = dopamine
        } else {
	  if (last_dopamine == 1 && active_syn_flag == 1) {

		  weight = weight + learning_rate_w_LTP * lthresh(ca_nmdai_max, lthresh_LTP) * hthresh(ca_nmdai_max, hthresh_LTP) * (1 - 2*sigmoidal(ca_nmdai_max, lthresh_LTP, steepness))
		  :lthresh_LTP = lthresh_LTP + learning_rate_thresh_LTP * lthresh(ca_nmdai_max, lthresh_LTP)* ((ca_nmdai_max-width)* delta - lthresh_LTP)
		  :lthresh_LTD = lthresh_LTD + learning_rate_thresh_LTD * lthresh(ca_nmdai_max, lthresh_LTP) * (cali_max - lthresh_LTD)
		  :hthresh_LTP = lthresh_LTP + width
          } else if (last_dopamine == -1 && active_syn_flag == 1) {

		  :weight = weight - learning_rate_w_LTD * lthresh(cali_max, lthresh_LTD)*weight
		  :lthresh_LTP = lthresh_LTP - learning_rate_thresh_LTP * lthresh(cali_max, lthresh_LTD)*(lthresh_LTP - max(ca_nmdai_max, lthresh_LTP_min))
		  :lthresh_LTD = lthresh_LTD - learning_rate_thresh_LTD * lthresh(cali_max, lthresh_LTD)*(lthresh_LTD - delta*cali_max)
		  :hthresh_LTP = lthresh_LTP + width

          }
          last_dopamine = dopamine		
          reset_max()
        }
	Cdur = Cdur_init + weight*Cdur_factor
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
:		 come again in Cdur with flag = current value of nspike
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

FUNCTION lthresh(conc, KD) {
:    lthresh = conc^n/(KD^n + conc^n) 
    lthresh = sigmoidal(conc,KD, steepness) 
}

FUNCTION hthresh(conc, KD) {
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
        active_syn_flag = 0
}

FUNCTION sigmoidal(x, x_offset, s) {
    sigmoidal = 1/(1+exp(-s *(x - x_offset)))
}
