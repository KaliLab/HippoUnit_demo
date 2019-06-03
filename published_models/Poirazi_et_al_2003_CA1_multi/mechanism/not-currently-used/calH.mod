TITLE L-calcium channel
: L-type calcium channel


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt  (ms)
	v (mV)
	celsius = 34	(degC)
	gcalbar = 0 (mho/cm2)
	ki=.001 (mM)
	cai=5.e-5 (mM)
	cao = 2  (mM)
:        tfa=2
        tfa=2
        eca = 140
}


NEURON {
	SUFFIX calH
	USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf,taum
}

STATE {
	m  
}

ASSIGNED {
	ica (mA/cm2)
        gcal  (mho/cm2) 
        minf
        taum
        }

INITIAL {
          rates(v)
          m = minf
	  gcal = gcalbar*m*m*h2(cai)
  }

BREAKPOINT {
	SOLVE states
	gcal = gcalbar*m*m*h2(cai)
	ica = gcal*ghk(v,cai,cao)

}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = -0.1967*(v-50)/(exp(-(v-50)/8.0)-1.0)
:	alpm = -0.1967*(v-30)/(exp(-(v-40)/8.0)-1.0)
}


FUNCTION betm(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
	betm = 0.046*exp(-v/20.73)
}


UNITSON
LOCAL facm

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        m = m + facm*(minf - m)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a
        a = alpm(v)
        taum = 1/(tfa*(a+betm(v)))
        minf = a/(a+betm(v))
        facm = (1 - exp(-dt/taum))

}



