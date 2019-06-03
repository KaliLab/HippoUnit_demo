TITLE  Na channel
: with recovery from inactivation and ir. M.Migliore BJ 1996
NEURON {
	SUFFIX namr
	USEION na READ ena WRITE ina
	NONSPECIFIC_CURRENT Ir
        RANGE  gnabar,vhalf,vhalfn,vhalfl,vvh,rmax, r, b
        GLOBAL ninf,linf,taul,taun,rinf,taur,inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {
        dt (ms)
	v (mV)
        ena=50 (mV)
	celsius = 22	(degC)
	gnabar=.01 (mho/cm2)
	vhalf=-50    (mV)
	zeta=20       (1)
        vhalfn=-30   (mV)
        vhalfr=-60   (mV)
        vhalfl=-57   (mV)
        a0l=0.1      (/ms)
        b0l=0.1      (/ms)
        a0n=2    (/ms)
        b0n=2    (/ms)
        a0r=0.0003    (ms)
        b0r=0.0003    (ms)
        zetan=-4    (1)
        zetar=12    (1)
        zetal=4    (1)
        gmn=0.9   (1)
        gmr=0.2   (1)
        gml=0.65   (1)
	lmax=1   (1)
	nmax=0.01   (1)
	rmax=3 (1)
	vvh=-58     (mv)
	vvs=2 	(1)
	b=0
}


STATE {
	n
        l
	r
	fr
}

ASSIGNED {
	Ir  (mA/cm2)
	ina (mA/cm2)
	inf
        ninf
        linf      
	rinf
        taul
        taun
	taur
}

INITIAL {
	rates(v,vhalf,vhalfn,vhalfl,vvh,rmax)
	n=ninf
	l=linf
	r=rinf
	fr=inf
}


BREAKPOINT {
	SOLVE states
	ina = gnabar*n*l*r*(v-ena)
	Ir = -ina*fr
:	Ir = 0
}

FUNCTION alp(v(mV),vf) { 
  alp = exp( 1.e-3*zeta*(v-vf)*9.648e4/(8.315*(273.16+celsius)))
}       

FUNCTION alpn(v(mV),vn) {
  alpn = exp(1.e-3*zetan*(v-vn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV),vn) {
  betn = exp(1.e-3*zetan*gmn*(v-vn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpv(v(mV),vh) {
         alpv = (1+b*exp((v-vh)/vvs))/(1+exp((v-vh)/vvs))
}

FUNCTION alpr(v(mV)) {
  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betr(v(mV)) {
  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV),vl) {
  alpl = exp(1.e-3*zetal*(v-vl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV),vl) {
  betl = exp(1.e-3*zetal*gml*(v-vl)*9.648e4/(8.315*(273.16+celsius))) 
}

LOCAL facn,facl,facr

:if state_borgka is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v,vhalf,vhalfn,vhalfl,vvh,rmax)
        n = n + facn*(ninf - n)
        l = l + facl*(linf - l)
        r = r + facr*(rinf - r)
	fr = inf
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV),vf,vn,vl,vh,rm) { :callable from hoc
        LOCAL a,q10
        q10=1.5^((celsius-22)/10)
        inf = 1/(1 + alp(v,vf))
        a = alpn(v,vn)
        ninf = 1/(1 + a)
        taun = betn(v,vn)/(q10*(a0n+b0n*a))
	if (taun<nmax) {taun=nmax}
        facn = (1 - exp(-dt/taun))
        a = alpl(v,vl)
        linf = 1/(1+ a)
        taul = betl(v,vl)/(q10*(a0l + b0l*a))
	if (taul<lmax) {taul=lmax}
        facl = (1 - exp(-dt/taul))
        rinf = alpv(v,vh)
        taur = betr(v)/(a0r+b0r*alpr(v))
	if (taur<rm) {taur=rm}
        facr = (1 - exp(-dt/taur))
}














