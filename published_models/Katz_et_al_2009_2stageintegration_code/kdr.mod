TITLE K-DR channel
: from Klee Ficker and Heinemann
: modified to account for Dax et al.
: M.Migliore 1997

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (mol) = (1)

}

NEURON {
        SUFFIX kdr
        USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar,ik
        RANGE ninf,taun
        GLOBAL nscale
}

PARAMETER {
        dt                      (ms)
        v                       (mV)
        ek                      (mV)    : must be explicitely def. in hoc
        celsius                 (degC)

        temp    = 24            (degC)

        gkdrbar = 0.003         (mho/cm2)

        vhalfn  = 13            (mV)
        a0n     = 0.02          (/ms)
        zetan   = -3            (1)
        gmn     = 0.7           (1)

        nmin    = 1             (ms)
        q10     = 1
        nscale  = 1
}

STATE {
        n
}

ASSIGNED {
        ik                      (mA/cm2)
        ninf
        gkdr                    (mho/cm2)
        taun                    (ms)
}

INITIAL {
        rates(v)
        n=ninf
        gkdr = gkdrbar*n
        ik = gkdr*(v-ek)
}        

BREAKPOINT {
        SOLVE states METHOD cnexp
        gkdr = gkdrbar*n
        ik = gkdr*(v-ek)
}

DERIVATIVE states {
        rates(v)
        n' = (ninf-n)/taun
}

FUNCTION alpn(v(mV)) {
        alpn = exp(zetan*(v-vhalfn)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/degC/mol)*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) {
        betn = exp(zetan*gmn*(v-vhalfn)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/degC/mol)*(273.16(degC)+celsius))) 
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-temp)/10(degC))
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
        if (taun<nmin) {taun=nmin}
        taun=taun/nscale
}














