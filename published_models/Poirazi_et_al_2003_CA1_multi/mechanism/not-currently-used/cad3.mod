TITLE decay of internal calcium concentration
:
: Internal calcium concentration due to calcium currents and pump.
: Differential equations.
:
: Simple model of ATPase pump with 3 kinetic constants (Destexhe 92)
:     Cai + P <-> CaP -> Cao + P  (k1,k2,k3)
: A Michaelis-Menten approximation is assumed, which reduces the complexity
: of the system to 2 parameters: 
:       kt = <tot enzyme concentration> * k3  -> TIME CONSTANT OF THE PUMP
:       kd = k2/k1 (dissociation constant)    -> EQUILIBRIUM CALCIUM VALUE
: The values of these parameters are chosen assuming a high affinity of 
: the pump to calcium and a low transport capacity (cfr. Blaustein, 
: TINS, 11: 438, 1988, and references therein).
:
: Units checked using "modlunit" -> factor 10000 needed in ca entry
:
: VERSION OF PUMP + DECAY (decay can be viewed as simplified buffering)
:
: All variables are range variables
:
: Written by Alain Destexhe, Salk Institute, Nov 12, 1992
:


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX cad3
        USEION ca READ ica, cai WRITE cai
        RANGE m
        GLOBAL depth,kt,kd,cainf,taur
}

UNITS {
        (molar) = (1/liter)                     : moles do not appear in units
        (mM)    = (millimolar)
        (um)    = (micron)
        (mA)    = (milliamp)
        (msM)   = (ms mM)
        FARADAY = (faraday) (coulomb)
}

:CONSTANT {
:        FARADAY = 96489         (coul)          : moles do not appear in units
:       FARADAY = 96.489        (k-coul)        
:}

PARAMETER {
        depth   = .1    (um)            : depth of shell
        :taur    = 700   (ms)            : rate of calcium removal
        taur    = 200   (ms)            : rate of calcium removal
        cainf   = 1e-8  (mM)
        cai		(mM)
        kt      = 1     (mM/ms)         : estimated from k3=.5, tot=.001
       : kt      = 5     (mM/ms)         : estimated from k3=.5, tot=.001
        kd      = 5e-4  (mM)            : estimated from k2=250, k1=5e5
}

STATE {
        m             (mM) 
}

INITIAL {
        m =  cainf 
}

ASSIGNED {
        ica             (mA/cm2)
        drive_channel   (mM/ms)
        drive_pump      (mM/ms)
}
        
BREAKPOINT {
        SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

        drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

        if (drive_channel <= 0. ) { drive_channel = 0. } : cannot pump inward
:        if (drive_pump <= -1e-04) { drive_pump = 0. } : cannot pump if cai too low

:       drive_pump = -tot * k3 * m / (m + ((k2+k3)/k1) )    : quasistat
        drive_pump = -kt * m / (m + kd )            : Michaelis-Menten

         m' = drive_channel + drive_pump + (cainf-m)/taur
         cai = m
         
}

