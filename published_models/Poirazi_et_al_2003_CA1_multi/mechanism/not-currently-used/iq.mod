
TITLE Q current

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt (ms)
	v (mV)
	erevq=-55   (mV)
	celsius = 30	(degC)
	gqbar=.003 (mho/cm2)
        vhalf=-93   (mV)
        a0=0.0009      (/ms)
        b0=0.0009     (/ms)
        zeta=5    (1)
        gq=0.4   (1)
}


NEURON {
	SUFFIX q
	NONSPECIFIC_CURRENT Iq
        RANGE Iq,gqbar
        GLOBAL inf,tau
}

STATE {
        q
}

ASSIGNED {
	Iq (mA/cm2)
        inf
        tau
}

INITIAL {
	rate(v)
	q = inf
}

BREAKPOINT {
	SOLVE state
	Iq = gqbar*q*(v-erevq)
}

FUNCTION alp(v(mV)) {
  alp = exp( 1.e-3*zeta*(v-vhalf)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bet(v(mV)) {
  bet = exp(1.e-3*zeta*gq*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
}

LOCAL facq

:if state_borgka is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {     : exact when v held constant; integrates over dt step
        rate(v)
        q = q + facq*(inf - q)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10=5^((celsius-23)/10)
        a = alp(v)
        inf = 1/(1 + a)
        tau = bet(v)/(q10*(a0+b0*a))
        facq = (1 - exp(-dt/tau))
}














