TITLE Kdr channel with slow inactivation
:Yiota Poirazi, 11/22/00


NEURON {
	SUFFIX kdr_inac
        USEION k READ ek WRITE ik
        RANGE  gkbar
	RANGE inf, fac, tau
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        
	celsius = 34	(degC)
	dt (ms)
	gkbar=0 (mho/cm2)
	ek = -77 (mV)
        v (mV)
}

STATE {
	 n l 
}

ASSIGNED {
        ik (mA/cm2)
        inf[2]
	fac[2]
	tau[2]
}


BREAKPOINT {
	SOLVE states
	ik = gkbar*n*n*l*(v - ek)
}

INITIAL {
	states()
        l=1
	ik = gkbar*n*n*l*(v - ek)
}

PROCEDURE calcg() {
	mhn(v*1(/mV))
	n = n + fac[0]*(inf[0] - n)
        l = l + fac[1]*(inf[1] - l)
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
		:varss = 1 / (1 + exp((v + 40)/(-2))) :K activation
                 varss = 1 / (1 + exp((v + 40)/(-3))) :K activation
         } else {	
		:varss = (1 + 1.6/(1+exp((v + 65)/12)))/2.52
                :varss = (1 + 1.8/(1+exp((v + 65)/10)))/2.75
                varss = (1 + 2.2/(1+exp((v + 53)/3)))/3.2
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
 		:vartau = 1.5     :K activation
                 vartau = 3     :K activation
        } else {
                vartau = 100   :K inactivation      		       
	VERBATIM
/*	printf("%g %g\n", _lv, _ltmp);*/
	ENDVERBATIM
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















