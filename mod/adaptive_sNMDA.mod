

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
	POINT_PROCESS adaptive_sNMDA
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA
	NONSPECIFIC_CURRENT  iNMDA
	RANGE mg, Cmax, eta, alpha, nmda_ca_fraction
        POINTER dopamine, stimulus_flag
	RANGE learning_rate_w_LTP, learning_rate_w_LTD, w0, wmax, wmin
	RANGE ca_nmdai_max, cali_max, active_syn_flag, Cdur_init, Cdur_factor
        RANGE last_dopamine, weight, KD1, KD2, n1, n2, KD_LTD, n_LTD, lthresh_LTP, hthresh_LTP
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

	learning_rate_w_LTP = 0.01
    	learning_rate_w_LTD = 0.01
    	wmax = 0.006 (uS)
    	wmin = 0.001 (uS)
        w0 = 0.00188 (uS)
	ca_nmdai_max = 0
	cali_max = 0
	active_syn_flag = 0
        nmda_ca_fraction = 0.15
        Cdur = 1.1
	n1 = 5
	n2 = 5
	KD1 = 1
	KD2 = 1
	KD_LTD = 1
	n_LTD = 1
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
        stimulus_flag
        cali            (mM)
	ca_nmdai        (mM)
        last_dopamine
        weight
        lthresh_LTP
        hthresh_LTP
}

STATE {Ron Roff}

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
	weight = w0
        last_dopamine = 0
        lthresh_LTP = 0
        hthresh_LTP = 0
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

		  lthresh_LTP = sigmoidal(ca_nmdai_max, KD1, n1)
		  hthresh_LTP = 1 - sigmoidal(ca_nmdai_max, KD2, n2)
		  weight = weight + learning_rate_w_LTP * lthresh_LTP * hthresh_LTP
		  		  
          } else if (last_dopamine == -1 && active_syn_flag == 1) {

		  weight = weight - learning_rate_w_LTD * sigmoidal(cali_max, KD_LTD, n_LTD) * weight

          }
          last_dopamine = dopamine		
          reset_max()
        }

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

FUNCTION lthresh(conc, KD, n) {
    lthresh = conc^n/(KD^n + conc^n) 
}

FUNCTION hthresh(conc, KD, n) {
    hthresh = KD^n/(KD^n + conc^n) 
}

FUNCTION sigmoidal(x, x_offset, s) {
    sigmoidal = 1/(1+exp(-s *(x - x_offset)))
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
