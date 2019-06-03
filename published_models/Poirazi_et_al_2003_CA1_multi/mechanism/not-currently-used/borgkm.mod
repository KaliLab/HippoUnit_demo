
TITLE Borg-Graham K-M channel

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt (ms)
        cai (mM)
	v (mV)
        ek (mV)
	tm=10
	celsius = 30	(degC)
	gkmbar=.003 (mho/cm2)
        vhalf=-55   (mV)
        a0=0.0002      (/ms)
        b0=0.0002     (/ms)
        zeta=-7    (1)
        gm=0.8   (1)
        st=1
}


NEURON {
	SUFFIX borgkm
	USEION k READ ek WRITE ik
        RANGE gkmbar
        GLOBAL inf,tau
}

STATE {
        m
}

ASSIGNED {
	ik (mA/cm2)
        inf
        tau
}

INITIAL {
	rate(v)
	m = inf
}


BREAKPOINT {
	SOLVE state
	ik = gkmbar*m^st*(v-ek)

}

FUNCTION alp(v(mV)) {
  alp = exp( 1.e-3*zeta*(v-vhalf)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bet(v(mV)) {
  bet = exp(1.e-3*zeta*gm*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
}

LOCAL facm

:if state_borgka is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {     : exact when v held constant; integrates over dt step
        rate(v)
        m = m + facm*(inf - m)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10=5^((celsius-22)/10)
        a = alp(v)
        inf = 1/(1 + a)
        tau = bet(v)/(q10*(a0+b0*a))
	if (tau<tm) {tau=tm}
        facm = (1 - exp(-dt/tau))
}














