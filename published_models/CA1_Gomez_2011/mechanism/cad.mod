TITLE decay of submembrane calcium concentration
:
: Internal calcium concentration due to calcium currents and decay.
: (decay can be viewed as simplified buffering)
:
:  This is a simple pool model of [Ca++]. 
:  cai' = drive_channel + (cainf-cai)/taur,
:  where the first term
:  drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
:  describes the change caused by Ca++ inflow into a compartment
:  with volume u (u is restricted to the volume of a submembrane shell).
: (Units checked using "modlunit" -> factor 10000 needed in ca entry.)
:
:  The second is a decay term that causes [Ca++] to decay exponentially 
:  (with a time constant taur) to the baseline concentration cainf
:  Simple first-order decay or buffering:
:
:       Cai + B <-> ...
:
:   which can be written as:
:
:       dCai/dt = (cainf - Cai) / taur
:
:   where cainf is the equilibrium intracellular calcium value (usually
:   in the range of 200-300 nM) and taur is the time constant of calcium 
:   removal.  The dynamics of submembranal calcium is usually thought to
:   be relatively fast, in the 1-10 millisecond range (see Blaustein, 
:   TINS, 11: 438, 1988).
:   Or, taur >= 0.1ms (De Schutter and Bower 1994),
:       taur <= 50 ms (Traub and Llinas 1977).
:
: Written by Alain Destexhe, Salk Institute, Nov 12, 1992
:
:

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ca
       	GLOBAL depth,cainf,taur
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulomb)
}


PARAMETER {
	depth	= .1	(um)		: depth of shell
	taur	= 200	(ms)		: rate of calcium removal
	cainf	= 100e-6(mM)
        cai	(mM)
}

STATE {
	ca		(mM) 
}

INITIAL {
	ca = cainf
}

ASSIGNED {
	ica		(mA/cm2)
        drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
:	SOLVE state METHOD euler
}

DERIVATIVE state { 

	drive_channel =  - (10000) * (ica) / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   : cannot pump inward          
	ca' = drive_channel/18 + (cainf-ca)/taur*7
	cai = ca
}







