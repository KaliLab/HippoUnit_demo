TITLE n-calcium channel
: n-type calcium channel

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX mycan
	USEION ca READ cai,cao WRITE ica
        RANGE gcanbar, gcan, hinf, minf, taum, tauh, facm, fach, ica
}

PARAMETER {
        dt  (ms)
	v (mV)
	celsius = 6.3	(degC)
	gcanbar = 0 (mho/cm2)
	cai=5.e-5 (mM)
	cao = 2  (mM)
        tadjm
        tadjh
        eca = 140
}


STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        gcan  (mho/cm2) 
        minf
        hinf
        taum
        tauh
        facm
        fach
}

INITIAL {
         tadjm= 0.3*3.55^((celsius-23.5)/10)
         tadjh= 0.4*2.8^((celsius-23.5)/10) 
	 rates(v)
         m = minf
         h = hinf
         gcan = gcanbar*m*m*h
}

BREAKPOINT {
	SOLVE states
	gcan = gcanbar*m*m*h
:	gcan = gcanbar*m*h
	ica = gcan*ghk(v,cai,cao,2)

}

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        m = m + facm*(minf - m)
        h = h + fach*(hinf - h)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { :callable from hoc
:        taum = (1/(exp((v+131.6)/-16.7)+exp((v+16.8)/18.2)) + 0.612) / tadjm
:        taum = 0.1*((exp((v+11)/-200)+exp((v+11)/200))) / tadjm
:        taum = 3*exp((v+11)/-200)/tadjm
:        taum = 100*(1- 0.95/(1+exp((v+11)/-200)))
        taum = 15*(0.055*(25.01 - v)/(exp((25.01-v)/10) - 1)+9.4*exp((-63.01-v)/20))
       minf = 1/(1+exp(-(v+11)/8.3))
        facm = (1 - exp(-dt/taum))
:        if (v<-80) 
:          { tauh = exp((v+467)/66.6) / tadjh }
:        else  
:          { tauh = (exp((v+21.88)/-10.52)+28) / tadjh }
:        tauh = 15*exp((v+11)/-200) / tadjh
:        tauh = 10*(1+ 2.5/(1+exp((v+11)/-10)))
         tauh = 10/(0.01*(10.01 - v)/(exp((10.01-v)/10) - 1)+20*exp((-110.01-v)/20))
        hinf = 1/(1+exp((v+44)/9.2))
        fach = (1 - exp(-dt/tauh))
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

