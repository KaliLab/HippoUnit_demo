TITLE Borg-Graham-like K-A channel
: INACTIVATING, M.Migliore, BJ, 1996

NEURON {
	SUFFIX borgkdr
	USEION k READ ek WRITE ik
        RANGE  gkdrbar,gkdr,vhalfl,vhalfn
	GLOBAL ninf,linf,taun,taul
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt (ms)
	v (mV)
        ek= -91 (mV)
	celsius = 30	(degC)
	gkdrbar=.003 (mho/cm2)
        vhalfn=-40   (mV)
        vhalfl=-60   (mV)
        a0l=0.001      (/ms)
        a0n=0.03      (/ms)
        b0n=0.03      (/ms)
        b0l=0.001      (/ms)
        zetan=-5    (1)
        zetal=2    (1)
        gmn=0.7   (1)
        gml=1.0   (1)
	nmax=0.3  (1)
}



STATE {
	n
        l
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf      
        gkdr
        taun
        taul
}

INITIAL {
 	rates(v,vhalfn,vhalfl)
        n=ninf
        l=linf
	gkdr = gkdrbar*n^3*l
	ik = gkdr*(v-ek)

}


BREAKPOINT {
	SOLVE states
	gkdr = gkdrbar*n^3*l
	ik = gkdr*(v-ek)

}

FUNCTION alpn(v(mV),vn) {
  alpn = exp(1.e-3*zetan*(v-vn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV),vn) {
  betn = exp(1.e-3*zetan*gmn*(v-vn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV),vl) {
  alpl = exp(1.e-3*zetal*(v-vl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV),vl) {
  betl = exp(1.e-3*zetal*gml*(v-vl)*9.648e4/(8.315*(273.16+celsius))) 
}
LOCAL facn,facl

:if state_borgka is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v,vhalfn,vhalfl)
        n = n + facn*(ninf - n)
        l = l + facl*(linf - l)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV),vn,vl) { :callable from hoc
        LOCAL a,q10
        q10=3^((celsius-30)/10)
        a = alpn(v,vn)
        ninf = 1/(1+a)
        taun = betn(v,vn)/(q10*(a0n+b0n*a))
	if (taun<nmax) {taun=nmax}
        facn = (1 - exp(-dt/taun))
        a = alpl(v,vl)
        linf = 1/(1+a)
        taul = betl(v,vl)/(q10*(a0l + b0l*a))
        facl = (1 - exp(-dt/taul))
}














