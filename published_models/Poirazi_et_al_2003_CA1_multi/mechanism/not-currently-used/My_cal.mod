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

NEURON {
	SUFFIX mycal
	USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf, taum, facm, ica
}

PARAMETER {
        dt  (ms)
	v (mV)
	celsius = 22.7	(degC)
	gcalbar = 0 (mho/cm2)
	cai = 5.e-5 (mM)
	cao = 2  (mM)
        tadjm
        eca = 140
}


STATE {
	m h 
}

ASSIGNED {
	ica   (mA/cm2)
        gcal  (mho/cm2) 
        minf
        taum
        facm
}

INITIAL {
         tadjm= 3.55^((celsius-23.5)/10)
         rates(v)
         m = minf
         gcal = gcalbar*m*m
}

BREAKPOINT {
	SOLVE states
	gcal = gcalbar*m*m
	ica = gcal*ghk(v,cai,cao,2)

}

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
:        taum = (1/(exp((v+131.6)/-16.7)+exp((v+16.8)/18.2)) + 0.612) / tadjm
        taum = 15*exp((v+4)/-200) / tadjm
        minf = 1./(1+exp(-(v+4)/6.0))
        facm = (1 - exp(-dt/taum))
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) {
        LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        e = w / (exp(w)-1)
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else
        : denominator is small -> Taylor series
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}

