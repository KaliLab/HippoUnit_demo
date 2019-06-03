TITLE HH channel (activity dependant attenuation)
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)
: Brannon-added attenuation 
: Poirazi-modified Kdr and Na threshold/time constants to make it more stable





NEURON {
	SUFFIX hha
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el
	RANGE ar2, vhalfs
	RANGE inf, fac, tau
	RANGE taus
	RANGE W
	GLOBAL taumin
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        a0r=0.0003    (ms)
        b0r=0.0003    (ms)

        zetar=12    (1)
	zetas=12   (1)
	gms=0.2    (1)
        gmr=0.2   (1)

	a0s=0.0003 (ms)
	vvs  = 2 (mV)
        vhalfr=-60   (mV)

	v (mV)
	W = 0.016 (/mV)    : this 1/61.5 mV. see the paper

	celsius = 34	(degC)
	dt (ms)
:	gnabar=.20 (mho/cm2)
:	gkbar=.12 (mho/cm2)
:	gl=.0001 (mho/cm2)
        gnabar=0 (mho/cm2)
	gkbar=0 (mho/cm2)
	gl=0 (mho/cm2)
	ena = 60 (mV)
	ek = -77 (mV)
	el = -70.0 (mV)	: steady state at v = -65 mV
	ar2 = 1.0

	b=0

	taumin=3

}

STATE {
	m h n l s
}

ASSIGNED {
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	inf[5]
	fac[5]
	tau[5]
}




BREAKPOINT {
	SOLVE states
:	ina = gnabar*m*m*h*s*(v - ena)
	ina = gnabar*m*m*h*s*(v - ena)
	ik = gkbar*n*n*l*(v - ek)
	il = gl*(v - el)
}

INITIAL {
	states()
	s=1
        l=1
        h=1
:	ina = gnabar*m*m*h*s*(v - ena)
	ina = gnabar*m*m*h*s*(v - ena)
	ik = gkbar*n*n*l*(v - ek)
	il = gl*(v - el)
}

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	n = n + fac[2]*(inf[2] - n)
        l = l + fac[3]*(inf[3] - l)
	s = s + fac[4]*(inf[4] - s)
}	

PROCEDURE states() {	: exact when v held constant
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}




FUNCTION varss(v, i) {
	LOCAL max, min,vhalf,smooth
	if (i==0) {
		varss = 1 / (1 + exp((v + 40)/(-3))) :Na activation
:		varss = 1 / (1 + exp((v + 39)/(-2.7))) :Na activation
	}
	else if (i==1) {
:		varss = 1 / (1 + exp((v + 45)/(3))) :Na inactivation
		varss = 1 / (1 + exp((v + 45)/(2.5))) :Na inactivation
	}
	else if (i==2) {	
:		varss = 1 / (1 + exp((v + 42)/(-2))) :K activation
		varss = 1 / (1 + exp((v + 40)/(-2))) :K activation
         }
	else if (i==3) {	
		varss = (1 + 1.6/(1+exp((v + 65)/12)))/2.52

	} else {
                :"s" activation system - Migliore 96 model
		max=1
		min=0
		smooth=2
		vhalf=-60
		varss =     alpv(v,vhalfr)
       }
}

FUNCTION alpr(v(mV)) { 
  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}


FUNCTION alpv(v(mV),vh) {
         alpv = (1+ar2*exp((v-vh)/vvs))/(1+exp((v-vh)/vvs))

}

FUNCTION betr(v(mV)) {
  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) }


FUNCTION vartau(v, i) {
	LOCAL alpha, beta, sum,tmp

	if (i==0) {
		vartau = 0.05  :Na activation tau
:                if ( v < -51) {
:                    vartau = (0.4/25.0)*(v +75)+0.1
:                } else if (v < 0) {
:                    vartau = (-0.3/50.0)*(v +50)+0.5
:                } else {
:                  vartau = 0.2
:                } 
	}
	else if (i==1) {
              vartau = 0.5 :Na inactivation tau
:	        vartau = (-7.0/50.0)*(v +75)+8
:                if (v > -25) {
:                  vartau = 1
:                } 
        }
	else if (i==2) {
 		vartau = 1.5     :K activation
        }
	else if (i==3) {
                vartau = 350   :K inactivation
      	} else {
	        tmp = betr(v)/(a0r+b0r*alpr(v))
	        if (tmp<taumin) {tmp=taumin}
	VERBATIM
/*	printf("%g %g\n", _lv, _ltmp);*/
	ENDVERBATIM
		vartau=tmp
       }
}	


PROCEDURE mhn(v) {LOCAL a, b :rest = -70
:	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 4 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
:	printf("v: %g inf[3]: %g\n", v, inf[3])
}















