TITLE HH-type Ca++ channel 
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)


NEURON {
	SUFFIX myca
	USEION ca READ eca WRITE ica
	RANGE gbar
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
	gbar=.12 (mho/cm2)
        eca = 140.0 (mV)	
	ar2 = 1.0

	b=0

	taumin=3

}

STATE {
	m h n s
}

ASSIGNED {
	ica (mA/cm2)
	inf[4]
	fac[4]
	tau[4]
}




BREAKPOINT {
	SOLVE states
	ica = gbar*m*m*h*s*(v - eca)
	}

INITIAL {
	states()
	s=1
	ica = gbar*m*m*h*s*(v - eca)
	}

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	n = n + fac[2]*(inf[2] - n)
	s = s + fac[3]*(inf[3] - s)
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
		varss = 1 / (1 + exp((v + 40)/(-3))) :Ca activation
	}
	else if (i==1) {
		varss = 1 / (1 + exp((v + 45)/(3))) :Ca inactivation
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
  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}


FUNCTION vartau(v, i) {
	LOCAL alpha, beta, sum,tmp

	if (i==0) {
		vartau = 0.05  :Ca activation tau
	}
	else if (i==1) {
		vartau = 0.5   :Ca inactivation tau
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
	FROM i=0 TO 2 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
:	printf("v: %g inf[2]: %g\n", v, inf[2])
}





