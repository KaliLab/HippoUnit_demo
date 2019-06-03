TITLE CaGk
: Calcium activated K channel.
: From Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

UNITS {
        (molar) = (1/liter)
}

UNITS {
        (mV) =  (millivolt)
        (mA) =  (milliamp)
        (mM) =  (millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

NEURON {
        SUFFIX cagk
        USEION ca READ cai
        USEION k READ ek WRITE ik
        RANGE gkbar
        GLOBAL oinf, tau
}

UNITS {
        FARADAY = (faraday)  (kilocoulombs)
        R = 8.313424 (joule/degC)
}

PARAMETER {
        celsius=20      (degC)
        v               (mV)
        gkbar=.01       (mho/cm2)       : Maximum Permeability
        cai = 1e-3      (mM)
        ek              (mV)
        dt              (ms)
        

        d1 = .84
        d2 = 1.
        k1 = .18        (mM)
        k2 = .011       (mM)
        abar = .28      (/ms)
        bbar = .48      (/ms)
}

ASSIGNED {
        ik              (mA/cm2)
        oinf
        tau             (ms)
}

STATE { o }             : fraction of open channels

BREAKPOINT {
        SOLVE state
        ik = gkbar*o*(v - ek)
}

LOCAL fac

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {     : exact when v held constant; integrates over dt step
        rate(v, cai)
        o = o + fac*(oinf - o)
        VERBATIM
        return 0;
        ENDVERBATIM
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { :callable from hoc
        alp = abar/(1 + exp1(k1,d1,v)/ca)
}

FUNCTION bet(v (mV), ca (mM)) (1/ms) { :callable from hoc
        bet = bbar/(1 + ca/exp1(k2,d2,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { :callable from hoc
        exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius))
}

PROCEDURE rate(v (mV), ca (mM)) { :callable from hoc
        LOCAL a
        a = alp(v,ca)
        tau = 1/(a + bet(v, ca))
        oinf = a*tau
        fac = (1 - exp(-dt/tau))
}
