TITLE Slow Ca-dependent potassium current
:
:   Ca++ dependent K+ current responsible for slow AHP

NEURON {
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gbar, po, ik
	GLOBAL m_inf, tau_m
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

ASSIGNED {       : parameters needed to solve DE
	v               (mV)
	celsius         (degC)
	ek              (mV)
	cai             (mM)           : initial [Ca]i
	ik              (mA/cm2)
	po
	m_inf
	tau_m           (ms)
}

PARAMETER {
	gbar    = 0   (mho/cm2)
        
: 	taumin  = 180    (ms)            : minimal value of the time cst
 	taumin  = 150    (ms)            : minimal value of the time cst
	b = 0.004		(/ms)
}


STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	po = m*m
	ik = gbar*po*(v - ek)    : potassium current induced by this channel
}

DERIVATIVE states {
	rates(cai)
	m' = (m_inf - m) / tau_m
}


INITIAL {
	rates(cai)
	m = m_inf
}


PROCEDURE rates(cai(mM)) {  LOCAL a 
	a = cai/b
 	m_inf = a/(a+1)
	tau_m = taumin+ 1(ms)*1(mM)*b/(cai+b)

}
