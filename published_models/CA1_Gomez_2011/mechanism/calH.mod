TITLE Ca L-type channel with high treshold of activation
: inserted in distal dendrites to account for distally
: restricted initiation of Ca++ spikes
: uses channel conductance (not permeability)
: written by Yiota Poirazi, 1/8/00 poirazi@LNC.usc.edu
:
: Updated to use Cvode  - Carl Gold, 08/10/03
: Updated by Yiota Poirazi  1/3/2005


NEURON {
	SUFFIX calH
	USEION ca READ cai, cao WRITE ica
        RANGE gcalbar, m, h, ica
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER { 
        ki     = 0.1  (mM)            : middle point of inactivation fct
        gcalbar = 0     (mho/cm2) : initialized conductance
}


ASSIGNED {                        : parameters needed to solve DE
	ica 		(mA/cm2)
        inf[2]				: inf/tau[0] = m
	tau[2]    	(ms)		: inf/tau[1] = h     
        v               (mV)
        celsius	        (degC)
	ecan             (mV)            : Ca++ reversal potential
	cai             (mM)      : initial internal Ca++ concentration
	cao             (mM)      : initial external Ca++ concentration
}

STATE {	
	m 
	h 
}                  

INITIAL {
        rates(v)
        m = 0    : initial activation parameter value
	h = 1    : initial inactivation parameter value
}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}



BREAKPOINT {
	SOLVE states METHOD cnexp
	ecan = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcalbar*m*m*m*h*h2(cai)*(v - ecan)       
}

DERIVATIVE states {
	rates(v)
	m' = (inf[0]-m)/tau[0]
	h' = (inf[1]-h)/tau[1]
}


PROCEDURE rates(v(mV)) {LOCAL a, b :rest = -70
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
	}
}


FUNCTION varss(v(mV), i) {
	if (i==0) { 
:             varss = 1 / (1 + exp((v+37)/(-1(mV))))  : Ca activation 
             varss = 1 / (1 + exp((v+37.7)/(-1(mV))))  : Ca activation 
	}
	else if (i==1) { 
             varss = 1 / (1 + exp((v+41)/(0.5(mV)))) : Ca inactivation 
	}
}

FUNCTION vartau(v(mV), i) (ms){
	if (i==0) {
:           vartau = 3.6(ms)  : activation variable time constant
           vartau = 3.5(ms)  : activation variable time constant
        }
	else if (i==1) {
:           vartau = 29(ms)   : inactivation variable time constant
           vartau = 20(ms)   : inactivation variable time constant
        }
}
