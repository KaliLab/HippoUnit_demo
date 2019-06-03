TITLE na3
: Na current for axon. No slow inact.
: M.Migliore Jul. 1997

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nax
	USEION na READ ena WRITE ina
	RANGE  gbar
	GLOBAL minf, hinf, mtau, htau,thinf, qinf
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	
								
	tha  =  -30	(mV)		: v 1/2 for act	
	qa   = 7.2	(mV)		: act slope (4.5)		
	Ra   = 0.4	(/ms)		: open (v)		
	Rb   = 0.124 	(/ms)		: close (v)		

	thi1  = -45	(mV)		: v 1/2 for inact 	
	thi2  = -45 	(mV)		: v 1/2 for inact 	
	qd   = 1.5	(mV)	        : inact tau slope
	qg   = 1.5      (mV)
	mmin=0.02	
	hmin=0.5			
	q10=2
	Rg   = 0.01 	(/ms)		: inact recov (v) 	
	Rd   = .03 	(/ms)		: inact (v)	

	thinf  = -50 	(mV)		: inact inf slope	
	qinf  = 4 	(mV)		: inact inf slope 

	ena		(mV)            : must be explicitly def. in hoc

	v 		(mV)
	dt		(ms)
	celsius=24	(degC)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
}
 

STATE { m h}

INITIAL {
trates(v)
m=minf  
h=hinf
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)

}

BREAKPOINT {
        SOLVE states
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
} 

LOCAL mexp, hexp, sexp

PROCEDURE states() {   
        trates(v)      
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(vm) {  
        LOCAL  a, b, qt
        qt=q10^((celsius-24)/10)
	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)
        mexp = 1 - exp(-dt/mtau)

	a = trap0(vm,thi1,Rd,qd)
	b = trap0(-vm,-thi2,Rg,qg)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf)/qinf))
        hexp = 1 - exp(-dt/htau)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	

        

