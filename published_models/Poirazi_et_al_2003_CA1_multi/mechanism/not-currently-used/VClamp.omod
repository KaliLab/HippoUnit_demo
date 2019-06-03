TITLE svclmp.mod

COMMENT

Single electrode Voltage clamp with three levels
------------------------------------------------

Series Resistance added; backards compatible, except parameters 
e0,vo0,vi0,gain,rstim,tau1,tau2 that no longer exist

Clamp is on at time 0, and off at time dur[0]+dur[1]+dur[2]. When clamp is off
the injected current is 0.  The clamp levels are amp[0], amp[1], amp[2].  i is
the injected current, vc measures the control voltage) Do not insert several
instances of this model at the same location in order to make level changes.
That is equivalent to independent clamps and they will have incompatible
internal state values.

The electrical circuit for the clamp is exceedingly simple:

        rs           Rin
vc ---'\/\/`---o---'\/\/`---o
               |            |
               |____| |_____|
                    | |
                     Cm

Note that since this is an electrode current model v refers to the internal
potential which is equivalent to the membrane potential v when there is no
extracellular membrane mechanism present but is v+vext when one is present. 
Also since i is an electrode current, positive values of i depolarize the
cell. (Normally, positive membrane currents are outward and thus hyperpolarize
the cell)

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 3

NEURON {
        POINT_PROCESS SEVClamp
        ELECTRODE_CURRENT i
        RANGE dur, amp, rs, vc, i
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (micromho)
}


PARAMETER {
        v (mV)
        rs = 1 (megohm)		: series resistance
}

ASSIGNED {
        i (nA)
        vc (mV)
        ic (nA)
        tc2 (ms)
        tc3 (ms)
	dur[NSTEP] (ms)
	amp[NSTEP] (mV)
        on
}

INITIAL {
        tc2 = dur[0] + dur[1]
        tc3 = tc2 + dur[2]
        on = 0
}

BREAKPOINT {
        SOLVE vstim
        if (on) {
                i = (vc - v)/rs
        }else{
                i = 0
        }
}

PROCEDURE vstim() {
        on = 1
        if (t < dur[0]) {
                vc = amp[0]
        }else if (t < tc2) {
                vc = amp[1]
        }else if (t < tc3) {
                vc = amp[2]
        }else {
                vc = 0
                on = 0
        }
        if (on) {
        }else{
                ic = 0
        }
        VERBATIM
        return 0;
        ENDVERBATIM
}
