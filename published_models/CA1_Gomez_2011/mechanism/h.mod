TITLE  H-current that uses Na ions
: Updated to use Cvode by Yiota Poirazi 12/1/2005

NEURON {
	SUFFIX h
        RANGE  gbar,vhalf, K, taun, ninf, g, ihi
:	USEION na READ ena WRITE ina      
	USEION hi READ ehi WRITE ihi VALENCE 1      
:	NONSPECIFIC_CURRENT i
}

UNITS {
	(um) = (micrometer)
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(pmho) = (picomho)
	(mmho) = (millimho)
}

:INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {              : parameters that can be entered when function is called in cell-setup
:        dt             (ms)
:	v              (mV)
        ena    = 50    (mV)
:        eh     = -10   (mV)
:        ehi     = -30   (mV)
        ehi     = -10   (mV)
	K      = 8.5   (mV)
	gbar   = 0     (mho/cm2)  : initialize conductance to zero
	vhalf  = -90   (mV)       : half potential
}	


STATE {                : the unknown parameters to be solved in the DEs
	n
}

ASSIGNED {             : parameters needed to solve DE
        v 
:	ina (mA/cm2)
	ihi (mA/cm2)
	ninf
	taun (ms)
	g
}

        


INITIAL {               : initialize the following parameter using states()
	rates()	
	n = ninf
	g = gbar*n
:	ina = g*(v-eh)
	ihi = g*(v-ehi)
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n
:	ina = g*(v-eh)  
	ihi = g*(v-ehi)  
}

DERIVATIVE states {
	rates()
        n' = (ninf - n)/taun
}

PROCEDURE rates() {  
 
 	if (v > -30) {
	   taun = 1
	} else {
           taun = 2*(1/(exp((v+145)/-17.5)+exp((v+16.8)/16.5)) + 5) :h activation tau

	}  
         ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))                  :steady state value
}



