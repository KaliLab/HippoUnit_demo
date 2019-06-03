TITLE K-A channel from Klee Ficker and Heinemann
: modified to account for Dax A Current ----------
: M.Migliore Jun 1997
: modified by Poirazi on 10/2/00 according to Hoffman_etal97 
: to account for I_A proximal (<100microns)
: (n) activation, (l) inactivation


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt (ms)
	v (mV)
        ek (mV)              : must be explicitely def. in hoc
	gkabar=.007 (mho/cm2)
        vhalfn= 11   (mV)
        vhalfl=-56   (mV)
        a0l=0.05      (/ms)
        a0n=0.05    (/ms)
}

NEURON {
	SUFFIX kap
	USEION k READ ek WRITE ik
        RANGE gkabar,gka
        GLOBAL ninf,linf,taul,taun
}

STATE {
	n
        l
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf
        taul
        taun
        gka
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
	gka = gkabar*n^4*l
	ik = gka*(v-ek)
}

BREAKPOINT {
	SOLVE states
	gka = gkabar*n^4*l
	ik = gka*(v-ek)
}

FUNCTION alpn(v(mV)) {
  alpn = -0.01*(v+21.3)/(exp((v+21.3)/-35)-1)
}


FUNCTION betn(v(mV)) {
  betn = 0.01*(v+21.3)/(exp((v+21.3)/35)-1)
}

FUNCTION alpl(v(mV)) {
  alpl = -0.01*(v+58)/(exp((v+58)/8.2)-1)
}

FUNCTION betl(v(mV)) {
  betl = 0.01*(v+58)/(exp((v+58)/-8.2)-1)
}

LOCAL facn,facl
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        n = n + facn*(ninf - n)
        l = l + facl*(linf - l)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,b
        a = alpn(v)
        b = betn(v)
        ninf = a/(a + b)
        taun = 0.2
        facn = (1 - exp(-dt/taun))
        a = alpl(v)
        b = betl(v)
        linf = a/(a + b)
        
        if (v > -20) {
	   taul = 5 + 2.6*(v+20)/10
        } else {
	   taul = 5
        }
        facl = (1 - exp(-dt/taul))
}
