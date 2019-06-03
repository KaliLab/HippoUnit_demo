TITLE Ca R-type channel with medium threshold for activation
: used in somatic regions. It has lower threshold for activation/inactivation
: and slower activation time constant
: than the same mechanism in dendritic regions
: uses channel conductance (not permeability)
: written by Yiota Poirazi on 3/12/01 poirazi@LNC.usc.edu

NEURON {
	SUFFIX somacar
	USEION ca READ cai, cao WRITE ica
        RANGE gcabar, m, h, ica
:	RANGE inf, fac, tau, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {     
        gcabar = 0      (mho/cm2) : initialized conductance
}


ASSIGNED {      : parameters needed to solve DE
        v               (mV)
 	celsius         (degC)
	ica             (mA/cm2)
	ecar             (mV)      : Ca++ reversal potential
        inf[2]
	fac[2]
	tau[2]
        cai		(mM)
        cao		(mM)
}

STATE {	
	m 
	h 
} 


INITIAL {
	m = 0    : initial activation parameter value
	h = 1    : initial inactivation parameter value
	rates(v)

}

BREAKPOINT {
:	rates(v)
	SOLVE states METHOD cnexp
        ecar = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar*m*m*m*h*(v - ecar)
}


DERIVATIVE states {
        rates(v)
	m' = (inf[0]-m)/tau[0]
	h' = (inf[1]-h)/tau[1]
}


PROCEDURE rates(v(mV)) {
	FROM i=0 TO 1 {
		tau[i] = vartau(i)
		inf[i] = varss(v,i)
	}
}



FUNCTION varss(v(mV), i) {
	if (i==0) {
	   varss = 1 / (1 + exp((v+60)/(-3))) :Ca activation
	}
	else if (i==1) {
           varss = 1/ (1 + exp((v+62)/(1)))   :Ca inactivation
	}
}

FUNCTION vartau(i) {
	if (i==0) {
           vartau = 100  : activation variable time constant
        }
	else if (i==1) {
           vartau = 5    : inactivation variable time constant
       }
	
}	
