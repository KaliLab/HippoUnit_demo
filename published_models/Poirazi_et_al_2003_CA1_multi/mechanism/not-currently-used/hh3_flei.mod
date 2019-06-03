TITLE HH channel
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)

: $Header: /disks/redondo/brannon/rs/bp/mechanism/hh3_flei.mod,v 1.6 1999/09/29 04:16:17 brannon Exp $

NEURON {
	SUFFIX hh3_flei
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el
	RANGE b_i
	GLOBAL W
	RANGE inf, tau
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
	gnabar=.20 (mho/cm2)
	gkbar=.12 (mho/cm2)
	gl=.0001 (mho/cm2)
	ena = 40 (mV)
	ek = -100 (mV)
	el = -70.0 (mV)	: steady state at v = -65 mV
	W = 0.016 (/mV)    : this 1/61.5 mV. see the paper
	b_i = 0.2
}

STATE {
	m h n s
}

ASSIGNED {
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	inf[4]
	tau[4]
}

LOCAL	fac[4]


BREAKPOINT {
	SOLVE states
	ina = gnabar*m*m*h*s*(v - ena)
	ik = gkbar*n*n*(v - ek)
	il = gl*(v - el)
}

INITIAL {
	states()
	m = inf[0]
	h = inf[1]
	n = inf[2]
	s = inf[3]
}

PROCEDURE states() {	: exact when v held constant
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	n = n + fac[2]*(inf[2] - n)
	s = s + fac[3]*(inf[3] - s)
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
	else if (i==2) {
		:varss = 0
		varss = 1 / (1 + exp((v + 40)/(-3))) :K activation
	} else {
                :"s" activation system - Migliore 96 model
:		varss = ( 1 + b_i*exp((v-58) / 2) ) / ( 1 + exp((v-58) / 2) )
		varss = 0.5
       }
}

FUNCTION vartau(v, i) {
	LOCAL alpha, beta, sum
	if (i==0) {
		vartau = 0.05  :Na activation tau
	}
	else if (i==1) {
		vartau = 0.5   :Na inactivation tau
	}
	else if (i==2) {
		vartau = 2     :K activation
	} else {
                :"s" activation system - Migliore 96 model
        alpha = 4e-4 * exp( 9.6 * (v+60) * W)
        beta  = 4e-4 * exp(-2.4 * (v+60) * W)
        sum = alpha + beta
	if ((1/sum) < 3.00) { 
	   vartau=3.00 
	} else {
           vartau =    1 / sum
	}
        printf ("%f %f\n", v, vartau)
      }
}	


PROCEDURE mhn(v) {LOCAL a, b :rest = -70
	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 3 {
		tau[i] = vartau(v, i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}
UNITSON
