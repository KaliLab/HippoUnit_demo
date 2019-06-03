TITLE K channel which is both voltage and Ca++ dependent
: fast activation/ slow inactivation
: Yiota Poirazi, 3/2/01


NEURON {
	SUFFIX vkca
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE gk, gbar, m, h, c, c_inf
      	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (molar) = (1/liter)
        (mM) = (millimolar)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}


PARAMETER {
        v                (mV)
        celsius = 34	(degC)
	dt               (ms)
        gbar = 0        (mho/cm2)
        ek = -80         (mV)
        cai              (mM)
        cac = 0.025      (mM)  
        gk             
               
}

STATE {
	m h c
}

ASSIGNED {
	ik
        inf[3]
	fac[3]
	tau[3]
        c_inf
}

BREAKPOINT {
	SOLVE states 
        gk = gbar*m*m*m*h*c*c
	ik = gk*(v - ek)       
	}

INITIAL {
        h = 1
        m = 0
        c = 0
	states()
        gk = gbar*m*m*m*h*c*c
	ik = gk*(v - ek)
        }


PROCEDURE calcg() {
	mhn(v*1(/mV), cai)
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	c = c + fac[2]*(inf[2] - c)
	}	

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car

         car = (cai/cac)^2
         c_inf = 1
:car / ( 1 + car )         
}
                            


PROCEDURE states() {	: exact when v held constant
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}


FUNCTION varss(v, i) { 

	if (i==0) {
          varss = 1 / (1 + exp((v+60)/(-1))) : activation
	}
	else if (i==1) {
          varss = 1 / (1 + exp((v+57)/(0.5))) : inactivation
       	}
        	
}

FUNCTION vartau(v, i) {
	
	if (i==0) {
         vartau = 2
        }
	else if (i==1) {
           vartau = 45
        }
        else if (i==2) {
         :  if (v < -55) {
         :     vartau = 10
         :  } else {
              vartau = 2
         :  }
        }
	
}	

PROCEDURE mhn(v, cai) {LOCAL a, b :rest = -70
:	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	
        FROM i=0 TO 1 {
           if (cai < 0.001 || cai > 0.2) {
              : inf[i] = 1
           } else {
               inf[i] = varss(v,i) 
           }
		tau[i] = vartau(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
        evaluate_fct(v,cai)
        inf[2] = c_inf
        tau[2] = vartau(v,2)
        fac[2] = (1 - exp(-dt/tau[2]))
}














