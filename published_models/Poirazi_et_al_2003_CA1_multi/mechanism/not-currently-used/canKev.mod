TITLE fast HVA calcium current

COMMENT
 fast high-voltage-activated calcium channel model
 Nif/AgTx/CgTx resistant VSCC from rat sensorimotor pyramidal cells
 Based on Lorenzon and Foehring (1995), J. Neurophysiol. 73(4):1430-1442

 Written by Kevin A. Archie, karchie@lnc.usc.edu
 (GHK code taken from Arthur Houweling's MyFirstNEURON models)

$Log: hvaccf.mod,v $
Revision 1.2  2000/09/27 22:45:41  karchie
Incorporated a few minor changes.

Revision 1.1  2000/09/21 17:50:33  karchie
Initial revision

ENDCOMMENT

VERBATIM

extern double nrn_ghk(double, double, double, double);

static const char rcsid[]="$Id: hvaccf.mod,v 1.2 2000/09/27 22:45:41 karchie Exp karchie $";

ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX ca
        USEION ca READ cai,cao WRITE ica
        RANGE pcabar, ica, m_inf, h_inf
}

UNITS {
        (mA)    = (milliamp)
        (mV)    = (millivolt)
        (mM)    = (milli/liter)
        FARADAY = 96480 (coul)
        R       = 8.314 (volt-coul/degC)
}

PARAMETER {
        v                       (mV)
        celsius                 (degC)
        dt                      (ms)
        cai             = 5.e-05(mM)
        cao             = 2.5   (mM)
        pcabar                  (cm/s)  
        tauM            = 5     (ms)
:        vHalfM          = -22   (mV)
:        slopeM          = 12    (mV)
        vHalfM          = 3   (mV)
        slopeM          = 8.3    (mV)
        tauH            = 0.8   (ms)
:        vHalfH          = -24   (mV)  : given slopeH, 80% inactivation @ -10mV
:        slopeH          = 10    (mV)  : close to slopeM (no data for this)
        vHalfH          = -39   (mV)  : given slopeH, 80% inactivation @ -10mV
        slopeH          = 9.2    (mV)  : close to slopeM (no data for this)
        tBase           = 23.5  (degC) : temperature for which tau is correct
:        mpow            = 2           : power of m in state equation
}

STATE {
        m
        h
}

ASSIGNED {
        ica             (mA/cm2)
        m_inf
        h_inf
}

INITIAL {
        LOCAL tadj

        : adjust rate constants based on temperature.
        : original experiments performed at room temperature
        : assumes that temperature remains constant through the sim
        tadj = 3^((celsius-tBase)/10)   : assume Q10 of 3
        tauM = tauM / tadj
        tauH = tauH / tadj

        : set initial values of state variables.
        rates(v)
        m = m_inf
        h = h_inf
}

BREAKPOINT {
        SOLVE states

:        ica = pcabar * pow(m,mpow) * h * nrn_ghk((v),(cai),(cao),2);
VERBATIM
        ica = pcabar * m * m * h * nrn_ghk((v),(cai),(cao),2);
ENDVERBATIM
}

PROCEDURE states() {
        rates(v)
        m = m + (1-exp(-dt/tauM))*(m_inf-m)
        h = h + (1-exp(-dt/tauH))*(h_inf-h)
}


PROCEDURE rates(v(mV)) {
        m_inf = 1/(1+exp(-(v-vHalfM)/slopeM))
        h_inf = 1/(1+exp((v-vHalfH)/slopeH))
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) {
        LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        e = w / (exp(w)-1)
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else
        : denominator is small -> Taylor series
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}







