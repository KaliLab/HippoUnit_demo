
: km.mod
: Potassium channel, Hodgkin-Huxley style kinetics
: Based on I-M (muscarinic K channel)
: Slow, noninactivating
: Author: Zach Mainen, Salk Institute, 1995, zach@salk.edu
:
: modified to use CVode --Carl Gold 08/12/03



NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE n, gbar,ik
	GLOBAL Ra, Rb, ninf, ntau
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
:	(pS) = (picosiemens)
:	(um) = (micron)
} 

PARAMETER {
	gbar = 0.03     (mho/cm2)
:	gbar = 10   	(pS/um2)	: 0.03 mho/cm2
	tha  = -30	(mV)		: v 1/2 for inf
	qa   = 9	(mV)		: inf slope		
	Ra   = 0.001	(/ms)		: max act rate  (slow)
	Rb   = 0.001	(/ms)		: max deact rate  (slow)
	temp = 23	(degC)		: original temp 	
	q10  = 2.3			: temperature sensitivity
	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	celsius		(degC)
	v 		(mV)
	ik 		(mA/cm2)
	ek		(mV)
	ninf
	ntau 		(ms)
	tadj
}
 

STATE { 
	n 
}

INITIAL {
        tadj = q10^((celsius - temp)/10(degC))  :temperature adjustment
	rates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
: 	ik = tadj* gbar*n * (v - ek)
	ik = (1e-4) *tadj* gbar*n * (v - ek)
}


DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/ntau
}


PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL a,b
        a = Ra * (v - tha)*1(/mV) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha)*1(/mV) / (1 - exp((v - tha)/qa))
        ntau = 1/(a+b)/tadj
	ninf = a*ntau
}

