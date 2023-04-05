TITLE simple NMDA receptors
NEURON {
	POINT_PROCESS adaptive_shom_NMDA
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA,mg, Cmax, eta, alpha, nmda_ca_fraction
	RANGE w0, threshf, synon, flagx, conc0, mltypeMin, mltype, kernel, kernel_LTD, lthresh_LTP, hthresh_LTP, lthresh_LTD
	RANGE ca_nmdai_max, cali_max, cati_max, cai_max
	RANGE learning_rate_w_LTD, learning_rate_w_LTP, learning_rate_thresh_LTD, learning_rate_thresh_LTP, learning_rate_thresh_KD_LTD, KD1, KD2, KD_LTD, n1, n2, n_LTD
	NONSPECIFIC_CURRENT iNMDA
	POINTER dopamine, stimulus_flag
	USEION ca_nmda READ ca_nmdai VALENCE 2	
	USEION cal READ cali VALENCE 2
	USEION cat READ cati VALENCE 2
	USEION ca READ cai WRITE ica VALENCE 2

}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	Cmax	= 1	 (mM)           : max transmitter concentration
	Cdur      = 1  (ms)		: transmitter duration (rising phase)
	Alpha	= 4 (/ms /mM)	: forward (binding) rate (4)
	Beta 	= 0.01   (/ms)   : backward (unbinding) rate
	Erev	= 15	 (mV)		: reversal potential
        mg   = 0      (mM)           : external magnesium concentration
        eta = 0  :change in the code
        alpha = 0 (/mV)
	gmax = 0   (uS)
        w0 = 0
        conc0=0 (mM)
        flagx=0 (1)
    	mltype=0 (mM)
    	mltypeMin=0.00015 (mM)
    	learning_rate_w_LTD=0 (1)
    	learning_rate_w_LTP=0 (1)
    	learning_rate_thresh_LTD=0 (1)
    	learning_rate_thresh_LTP=0 (1)
    	learning_rate_thresh_KD_LTD = 0
    	threshf=0 (mM)    :change in the code
    	threshltp=0.02(mM)
    	f=0   (ms)
    	conc   (mM)
	n1 = 5
	n2 = 5
	KD1 = 1
	KD2 = 1
	KD_LTD = 1
	n_LTD = 1
	nmda_ca_fraction = 0.175
	kernel_LTD = 0
        cai_max = 0
	ca_nmdai_max = 0
	cali_max = 0
	cati_max = 0


}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iNMDA 		(nA)		: current = g*(v - e)
	g 		(uS)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
        B                       : magnesium block
	ica        (nA)
	dopamine
        cali            (mM)
	ca_nmdai        (mM)
	cai  (mM)
	cati
	stimulus_flag
	kernel
	lthresh_LTP
	hthresh_LTP
	lthresh_LTD
	
}

STATE {Ron Roff weight thresh}

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
	weight = w0
	threshf = KD1
        thresh = threshf
   

}

BREAKPOINT {
	SOLVE release METHOD cnexp
    B = mgblock(v)
	g = (Ron + Roff)* gmax*B
	iNMDA = g*(v - Erev)
    ica = 0.175*iNMDA   :(5-10 times more permeable to Ca++ than Na+ or K+, Ascher and Nowak, 1988)
    iNMDA = 0.825*iNMDA
    kernel = supra(dopamine,weight,ca_nmdai,cali,cai,g,thresh)
    weight=weight+supra(dopamine,weight,ca_nmdai,cali,cai,g,thresh)
    thresh=thresh+sTresh(dopamine,thresh,conc0)
    lthresh_LTP = thresh
    hthresh_LTP = thresh
    lthresh_LTD = thresh
    ca_nmdai_max = conc0
    cali_max = mltype
   
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


FUNCTION supra(dopam,we(uS),conc(mM),calic(mM),cai(mM),g_p(uS),tresh(mM))(uS/ms) {
    UNITSOFF
    
    
    if (g_p < 1e-7) {
       conc0=0
       mltype=0
    } else {
    
       if (conc>conc0 ){
          conc0=conc
       }
       if (calic>mltype ){
          mltype=calic
       }
    }
    
	supra =((learning_rate_w_LTD*we*mltype*funcCalMin(mltype)*(dopam-1))+(learning_rate_w_LTP*funcCal(conc0,tresh)*(dopam+1)))*trap(g_p)*dopam*dopam
	
	UNITSON    
}
FUNCTION sTresh(dopam,tresh(mM),conc0)(uM/ms) {
    UNITSOFF
    sTresh =((-learning_rate_thresh_LTD*funcCalt(conc0,tresh)*(dopam-1))-learning_rate_thresh_LTP*(funcCalt(conc0,tresh)*(dopam+1)))*trap(g)*dopam*dopam
	UNITSON    
}


FUNCTION trap(g_p(uS))() {
    UNITSOFF
	if (g_p < 1e-7) {
		trap = 0
	} else {
	    trap = 1
        }
    UNITSON
}


FUNCTION funcCal(calnm(mM),tresh(mM))() {
    UNITSOFF
    funcCal= (1-(1 / (1 + exp((-1000*calnm + (1000*tresh))/1))))*(1 / (1 + exp((-1000*calnm+(1000*tresh))/1))) 
    UNITSON
}
FUNCTION funcCalt(calnm(mM),tresh(mM))() {
    UNITSOFF
    funcCalt= (1-(1 / (1 + exp((-1000*calnm + (1000*tresh))/3))))*(1 / (1 + exp((-1000*calnm+(1000*tresh))/3)))
    UNITSON
}
FUNCTION funcCalteg(calnm(mM),tresh(mM))() {
    UNITSOFF
    funcCalteg= (-1 / (1 + exp((-1000*calnm+(1000*tresh))/2.5))) 
    UNITSON
}
FUNCTION funcCalMin(calnm(mM))() {
    UNITSOFF
    funcCalMin=1 / (1 + exp((-Norm(calnm)+7)/1))
    UNITSON
}


FUNCTION Norm(calnm(mM))() {
    UNITSOFF
    Norm= 100000*calnm
    UNITSON
}
FUNCTION funcCalMin1(calnm(mM))() {
    UNITSOFF
    funcCalMin1= 1 / (1 + exp((-100000*calnm+(100000*mltypeMin))/1)) 
    UNITSON
}
