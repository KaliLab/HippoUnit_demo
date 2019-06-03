COMMENT

This .mod file MUST use the euler method of differentiation. Other
methods make numerous calls to the check() block within a single
timestep and hence do unintuitive things to the state variables there.

ENDCOMMENT

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
	(mM) = (milli/liter)
      (umho) = (micromho)
        (nA) = (nanoamp)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
	POINT_PROCESS abbott_nmda
	POINTER vpre
	NONSPECIFIC_CURRENT i
	RANGE d, s
	RANGE tauD, tauS, taug, G
	RANGE erev, thresh
}

STATE {
	D 
	S
	g 		(umho)		: conductance
	h 		(umho)		: conductance
}

PARAMETER {
	v    (mV)
	erev = 0 (mV)
	thresh = 0.5 (mV) : useful for 0-1 files
	d = 0.4
	s = 1.0
	tauD = 300 (ms)
	tauS = 20e3 (ms)
	taug = 2 (ms)
	G    = 10e-12 (umho)
	eta     = 0.33  (/mM)
	mag     = 1     (mM)
	gamma   = 0.06  (/mV)
}	

ASSIGNED {

	i     (nA)
	firing
	vpre  (mV)
}

INITIAL {
	firing = 0
	D = 1
	S = 1
	g = 0
}

BREAKPOINT {

   SOLVE depression METHOD euler : so now we're clinical psychiatrists
   aux()
}

DERIVATIVE depression {
   check()
   D' = (1/tauD) * ( 1 - D ) 
   S' = (1/tauS) * ( 1 - S )
   g' = (1/taug) * (-g) 
:   printf("D: %f, S: %f, g: %f\n", D, S, g)
}

PROCEDURE aux() {
	h = 1/(1 + eta * mag * exp( - (gamma * v)))
	i = (g * h * (v - erev))
}

PROCEDURE check() {
:	  printf ("--------------------------------------------------\n")
:	  printf ("t: %f\n", t)


	if (firing && (vpre < thresh)) {
		firing = 0
		}
	if ((vpre >= thresh) && !firing) {
		firing = 1
		D = d * D
		S = s * S
		g = g + G * D * S
:		printf(" **** firing: %f \n", firing)
		}

:	  printf("vpre: %f, firing: %f\n", vpre, firing)
}

