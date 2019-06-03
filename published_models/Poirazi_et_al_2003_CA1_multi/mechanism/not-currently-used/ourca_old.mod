TITLE Ca R-type channel 
: Yiota Poirazi, 11/13/00


NEURON {
	SUFFIX car
	USEION ca READ eca WRITE ica
        RANGE gcabar, m, h
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 34	(degC)
	dt (ms)
        gcabar=0 (mho/cm2)
	eca = 140 (mV)
        }

STATE {
	m h
}

ASSIGNED {
	ica (mA/cm2)
        inf[2]
	fac[2]
	tau[2]
}

BREAKPOINT {
	SOLVE states
	ica = gcabar*m*m*m*h*(v - eca)
        :ica = gcabar*m*m*m*(v - eca )
	}

INITIAL {
        h=1
        m=0
	states()
	ica = gcabar*m*m*m*h*(v - eca)
	:ica = gcabar*m*m*m*(v - eca)
        }

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	}	

PROCEDURE states() {	: exact when v held constant
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}




FUNCTION varss(v, i) {
	if (i==0) {
		varss = 1 / (1 + exp((v+49)/(-3))) :Ca activation
	}
	else if (i==1) {
                 varss = 1/ (1 + exp((v+52)/(1))) :Ca inactivation
	}
	
}

FUNCTION vartau(v, i) {
	
	if (i==0) {
        : vartau = 5+9/(exp((v+35)/5) + exp(-(v+35)/5))
          vartau = 50
        }
	else if (i==1) {
         : vartau = 7+ 50/(exp((v+35)/10) + exp(-(v+35)/20))
          vartau = 10
       }
	
}	

PROCEDURE mhn(v) {LOCAL a, b :rest = -70
:	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}















