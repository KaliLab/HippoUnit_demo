TITLE HH channel
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)

NEURON {
	SUFFIX hh3
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el
	RANGE gna,gk
	GLOBAL inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	celsius = 37	(degC)
	dt (ms)
:	gnabar=.20 (mho/cm2)
:	gkbar=.12 (mho/cm2)
:	gl=.0001 (mho/cm2)
	gnabar=0 (mho/cm2)
	gkbar=0(mho/cm2)
	gl=0 (mho/cm2)
	ena = 40 (mV)
	ek = -100 (mV)
	el = -70.0 (mV)	: steady state at v = -65 mV
}
STATE {
	m h n
}
ASSIGNED {
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	inf[3]
	gna
	gk
}
LOCAL	fac[3]


BREAKPOINT {
	SOLVE states
	gna = gnabar*m*m*h
	ina = gna*(v - ena)
	gk = gkbar*n*n
	ik =  gk*(v - ek)
	il = gl*(v - el)
}

INITIAL {
	mhn(v*1(/mV))		
	m = inf[0]
	h = inf[1]
	n = inf[2]
}

PROCEDURE states() {	: exact when v held constant
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	n = n + fac[2]*(inf[2] - n)
	VERBATIM
	return 0;
	ENDVERBATIM
}

UNITSOFF
FUNCTION expM1(x,y) {
	if (fabs(x/y) < 1e-6) {
		expM1 = y*(1 - x/2)
	}else{
		expM1 = x/(exp(x/y) - 1)
	}
}

FUNCTION varss(v, i) {
	if (i==0) {
		varss = 1 / (1 + exp((v + 40)/(-3))) :Na activation
	}
	else if (i==1) {
		varss = 1 / (1 + exp((v + 45)/(3))) :Na inactivation
	}
	else {
		:varss = 0
		varss = 1 / (1 + exp((v + 40)/(-3))) :K activation
	}
}

FUNCTION vartau(i) {
	if (i==0) {
		vartau = 0.05  :Na activation tau
	}
	else if (i==1) {
		vartau = 0.5   :Na inactivation tau
	}
	else {
		vartau = 2     :K activation
	}
}	


PROCEDURE mhn(v) {LOCAL a, b, tau :rest = -70
	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 2 {
		tau = vartau(i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau))
	}
}
UNITSON
