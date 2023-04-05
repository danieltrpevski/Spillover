TITLE simple NMDA receptors
NEURON {
	POINT_PROCESS adaptive_zahra_NMDA
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA, mg, Cmax, eta, alpha, w0, weight
	RANGE treshf, synon, conc0, mltypeMin, mltype, rate_ltd, rate_ltp, rate_ltd_thrsh, rate_ltp_tresh, tremin
	NONSPECIFIC_CURRENT  iNMDA
	POINTER dopamine
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
	Cdur      =1  (ms)		: transmitter duration (rising phase)
	Alpha	= 4 (/ms /mM)	: forward (binding) rate (4)
	Beta 	= 0.01   (/ms)   : backward (unbinding) rate
	Erev	= 15	 (mV)		: reversal potential
    mg   = 0      (mM)           : external magnesium concentration
    eta = 0  :change in the code
    alpha = 0 (/mV)
	gmax = 0   (uS)
    w0 = 0
    conc0=0 (mM)
    mltype=0 (mM)
    mltypeMin=0.00015 (mM)
    rate_ltd=0 (1)
    rate_ltp=0 (1)
    rate_ltd_thrsh=0 (1)
    rate_ltp_tresh=0 (1)
    treshf=0 (mM)    :change in the code
    treshltp=0.02(mM)
    tremin=0 (mM)
    f=0   (ms)
    



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
    cali            (mM)
	ca_nmdai        (mM)
	
}

STATE {Ron Roff weight tresh}

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
	weight = w0
    tresh = treshf
   

}

BREAKPOINT {
	SOLVE release METHOD cnexp
    B = mgblock(v)
	g = (Ron + Roff)* gmax*B
	iNMDA = g*(v - Erev)
    ica_nmda = 0.175*iNMDA   :(5-10 times more permeable to Ca++ than Na+ or K+, Ascher and Nowak, 1988)
    iNMDA = 0.825*iNMDA
    weight=weight+supra(dopamine,weight,ca_nmdai,cali,g,tresh)
    tresh=tresh+sTresh(dopamine,tresh,conc0)

   
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


FUNCTION supra(dopam,we(uS),conc(mM),calic(mM),g_p(uS),tresh(mM))(uS/ms) {
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
    
	supra =((rate_ltd*we*mltype*funcCalMin(mltype)*(dopam-1))+(rate_ltp*funcCal(conc0,tresh)*(dopam+1)))*trap(g_p)*dopam*dopam
	
	UNITSON    
}
FUNCTION sTresh(dopam,tresh(mM),conc0)(uM/ms) {
    UNITSOFF
    sTresh =((rate_ltd_thrsh*(tresh-tremin)*(dopam-1))+rate_ltp_tresh*(funcCalsecd(conc0,tresh)*(dopam+1)))*trap(g)*dopam*dopam
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
    funcCal= (1-(1 / (1 + exp((-1000*calnm + (1000*tresh))/2.5))))*(1 / (1 + exp((-1000*calnm+(1000*tresh))/2.5))) 
    UNITSON
}
FUNCTION funcCalsecd(calnm(mM),tresh(mM))() {
    UNITSOFF
    funcCalsecd= (1-(1 / (1 +exp((-1000*calnm+ (1000*tresh))/2.5))))*(1 / (1 + exp((-1000*calnm+(1000*tresh))/2.5)))* (1-(2 / (1 + exp((-1000*calnm+ (1000*tresh))/2.5))))
    UNITSON
}
FUNCTION funcCalMin(calnm(mM))() {
    UNITSOFF
    funcCalMin= 1 / (1 + exp((-100000*calnm+(100000*mltypeMin))/1)) 
    UNITSON
}



