COMMENT
Updated Exp2Syn synapse with Mg-blocked nmda channel.

Defaul values of parameters (time constants etc) set to match synaptic channels in 
striatal medium spiny neurons (Du et al., 2017; Chapman et al., 2003; Ding et al., 2008).

Robert . Lindroos @ ki . se

original comment:
________________
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT



NEURON {
	POINT_PROCESS adaptive_glutamate_cshom
	RANGE tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda
	RANGE erev_ampa, erev_nmda, g, i
	NONSPECIFIC_CURRENT i
	
	RANGE i_ampa, i_nmda, g_ampa, g_nmda, I, G, mg, q, alpha, eta
	RANGE w0, NMDA_AMPA_ratio
	RANGE weight, lthresh_LTP, lthresh_LTD, hthresh_LTP, last_dopamine, lthresh_LTP_min
        RANGE hthresh_LTP_const, hthresh_LTP_0, hthresh_max, n, delta, width, steepness_LTP, steepness_LTD
	RANGE learning_rate_w_LTP, learning_rate_w_LTD, thresh_LTP_0,  learning_rate_thresh_LTP, thresh_LTD_0, learning_rate_thresh_LTD, hthresh_LTP_0
	RANGE ca_nmdai_max, cali_max, active_syn_flag, nmda_ca_fraction, delta_LTP
	POINTER dopamine, stimulus_flag
	USEION ca_nmda READ ca_nmdai WRITE ica_nmda VALENCE 2	
	USEION cal READ cali VALENCE 2
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	erev_ampa        = 0.0       (mV)
	erev_nmda 	 = 15.0 (mV)
	tau1_ampa   = 1.9       (ms)
    	tau2_ampa   = 4.8       (ms)  : tau2 > tau1
    	tau1_nmda   = 5.52      (ms)  : old value was 5.63
    	tau2_nmda   = 231       (ms)  : tau2 > tau1
    
    	mg          = 1         (mM)
    	alpha       = 0.062
    	q           = 2
    	eta 	= 18
	NMDA_AMPA_ratio = 1
	w0 = 0.01

    	learning_rate_w_LTP = 0.01
	learning_rate_w_LTD = 0.01
	
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
        steepness_LTP = 0.25
	steepness_LTD = 2.5
        nmda_ca_fraction = 0.15

	ca_nmdai_max = 0
	cali_max = 0
	active_syn_flag = 1e-6
}


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor_nmda
	factor_ampa
	i_ampa
	i_nmda
	g_ampa
	g_nmda
	block
	I
	G

	stimulus_flag
	dopamine
        last_dopamine
        ica_nmda (nA)
	ca_nmdai (mM)
	cali (mM)

	weight
	lthresh_LTP
	lthresh_LTD
	hthresh_LTP
	delta_LTP
}


STATE {
	A (uS)
	B (uS)
	C (uS)
	D (uS)
}



INITIAL {
	LOCAL tp
	if (tau1_nmda/tau2_nmda > .9999) {
		tau1_nmda = .9999*tau2_nmda
	}
	if (tau1_ampa/tau2_ampa > .9999) {
		tau1_ampa = .9999*tau2_ampa
	}
	
	: NMDA
	A           = 0
	B           = 0
	tp          = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)
	factor_nmda = -exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda)
	factor_nmda = 1/factor_nmda
	
	: AMPA
	C           = 0
	D           = 0
	tp          = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)
	factor_ampa = -exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa)
	factor_ampa = 1/factor_ampa
	
	weight = w0 
	lthresh_LTP = thresh_LTP_0
	lthresh_LTD = thresh_LTD_0
	hthresh_LTP = hthresh_LTP_0
	active_syn_flag = 0
        last_dopamine = 0
}




BREAKPOINT {
	SOLVE state METHOD cnexp
	
	: NMDA
	g_nmda = (B - A)*weight*NMDA_AMPA_ratio
	block  = MgBlock()
	i_nmda = g_nmda * (v - erev_nmda) * block
        ica_nmda = nmda_ca_fraction*i_nmda
        i_nmda = (1 - nmda_ca_fraction)*i_nmda
	
	: AMPA
	g_ampa = (D - C)*weight
	i_ampa = g_ampa * (v - erev_ampa)
	
	: total current
	G = g_ampa + g_nmda
	I = i_ampa
        i = I

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

		  weight = weight - learning_rate_w_LTD * cali_max * sigmoidal(cali_max, lthresh_LTD, steepness_LTD) * weight 
		  lthresh_LTP = lthresh_LTP - learning_rate_thresh_LTP * (lthresh_LTP - lthresh_LTP_min)

          }
          last_dopamine = dopamine		
          reset_max()
}

	:        delta_LTP = lthresh(ca_nmdai_max, lthresh_LTP, steepness_LTP) * hthresh(ca_nmdai_max, lthresh_LTP, steepness_LTP)
	:	  weight = weight + learning_rate_w_LTP * lthresh(ca_nmdai_max, lthresh_LTP, steepness_LTP) * hthresh(ca_nmdai_max, hthresh_LTP, steepness_LTP)
	:	  lthresh_LTP = lthresh_LTP + learning_rate_thresh_LTP * delta_LTP * (1 - 2*sigmoid(ca_nmdai_max, lthresh_LTP, steepness_LTP))
	:	  lthresh_LTP = max(lthresh_LTP, lthresh_LTP_min)
		  
		  :lthresh_LTP = lthresh_LTP + learning_rate_thresh_LTP * lthresh(ca_nmdai_max, lthresh_LTP, steepness_LTP)* ((ca_nmdai_max-width)* delta - lthresh_LTP)
		  :lthresh_LTD = lthresh_LTD + learning_rate_thresh_LTD *lthresh(ca_nmdai_max, lthresh_LTP, steepness_LTP) * (cali_max*delta - lthresh_LTD)
		  :hthresh_LTP = lthresh_LTP + width

          
	:	  weight = weight - learning_rate_w_LTD * cali_max * lthresh(cali_max, lthresh_LTD, steepness_LTD) * weight
	:	  lthresh_LTP = lthresh_LTP - learning_rate_thresh_LTP * lthresh(cali_max, lthresh_LTD, steepness_LTD)*(lthresh_LTP - max(ca_nmdai_max, lthresh_LTP_min))
	:	  lthresh_LTP = max(lthresh_LTP, lthresh_LTP_min)
		  
		  :lthresh_LTP = lthresh_LTP - learning_rate_thresh_LTP *lthresh(cali_max, lthresh_LTD, steepness_LTD)*(lthresh_LTP - max(ca_nmdai_max, lthresh_LTP_min))
		  :lthresh_LTD = lthresh_LTD - learning_rate_thresh_LTD *lthresh(cali_max, lthresh_LTD, steepness_LTD)*(lthresh_LTD - delta*cali_max)
		  :hthresh_LTP = lthresh_LTP + width

         
}



DERIVATIVE state {
	A' = -A/tau1_nmda*q
	B' = -B/tau2_nmda*q
	C' = -C/tau1_ampa
	D' = -D/tau2_ampa
}



NET_RECEIVE(dummy (uS)) {
	active_syn_flag = 1
	
	A = A + factor_nmda
	B = B + factor_nmda
	C = C + factor_ampa
	D = D + factor_ampa
}


FUNCTION MgBlock() {
    
    MgBlock = 1 / (1 + mg * eta * exp(-alpha * v)  )
    
}

FUNCTION lthresh(conc, KD, steepness) {
:    lthresh = conc^n/(KD^n + conc^n)
    lthresh = sigmoidal(conc,KD, steepness)  
}

FUNCTION hthresh(conc, KD, steepness) {
:    hthresh = KD^n/(KD^n + conc^n)
    hthresh = 1 - sigmoidal(conc, KD, steepness)  
}


FUNCTION reset_max() {
	ca_nmdai_max = 0
        cali_max = 0
	active_syn_flag = 0
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

FUNCTION sigmoidal(x, x_offset, s) {
    sigmoidal = 1/(1+exp(-s *(x - x_offset)))
}
