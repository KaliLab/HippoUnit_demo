TITLE  H-current

NEURON {
	SUFFIX h
        RANGE  gbar,vhalf, K, taun, ninf, g  :range variables
	USEION na READ ena WRITE ina        : which ions to use and how to treat them ena = reversal pot for Na
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

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {                    : parameters that can be entered when function is called in cell-setup
        dt             (ms)
	v              (mV)
        ena=0            (mV)

	K      = 8.5   (mV)
	gbar   = 0.1   (mmho/cm2) : this is the somatic value. dendritic value is 6x higher
	vhalf  =-90    (mV)
}	


STATE {                 : the unknown parameters to be solved in the DEs
	n
}

ASSIGNED {             : parameters neede to solve DE
	ina (mA/cm2)
	ninf
	taun (ms)
	g
}

        


INITIAL {      : initialize the following parameter using states()
	states()	
	n = ninf
	g = gbar*n
	ina = g*(v-ena)*(0.001)
}


BREAKPOINT {
	SOLVE h METHOD derivimplicit
	g = gbar*n
	ina = g*(v-ena)*(0.001)
}

DERIVATIVE h {
	states()
        n' =(ninf - n)/taun
}

PROCEDURE states() {  
 
 	if (v > -30) {
	   taun = 1
	} else {
:           taun = -0.5(ms/mV)*v - 10(ms)
            taun = -0.5(ms/mV)*v - 10(ms)
	}
 
        ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))

}



