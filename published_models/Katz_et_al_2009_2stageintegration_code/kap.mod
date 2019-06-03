TITLE K-A channel from Klee Ficker and Heinemann
: modified to account for Dax A Current ----------
: M.Migliore Jun 1997

NEURON {
        SUFFIX kap
        USEION k READ ek WRITE ik
        RANGE gkabar,gka,ik
        RANGE ninf,linf,taul,taun
        RANGE vhalfn,vhalfl
        GLOBAL lmin,nscale,lscale
}

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (mol) = (1)
}

PARAMETER {
        dt                              (ms)
        v                               (mV)
        ek                              (mV)
        celsius                         (degC)

        temp    = 24                    (degC)

        gkabar                          (mho/cm2)

        vhalfn  = 11                    (mV)
        a0n     = 0.05                  (/ms)
        zetan   = -1.5                  (1)
        gmn     = 0.55                  (1)
        pw      = -1                    (1)
        tq      = -40                   (mV)
        qq      = 5                     (mV)
        nmin    = 0.1                   (ms)
        nscale  = 1

        vhalfl  = -56                   (mV)
        a0l     = 0.05                  (/ms)
        zetal   = 3                     (1)
        lmin    = 2                     (ms)
        lscale  = 1

        q10     = 5
}

STATE {
        n
        l
}

ASSIGNED {
        ik (mA/cm2)
        ninf
        linf      
        taul  (ms)
        taun  (ms)
        gka   (mho/cm2)
        qt
}

INITIAL {
        rates(v)
        n=ninf
        l=linf
        gka = gkabar*n*l
        ik = gka*(v-ek)
}        

BREAKPOINT {
        SOLVE states METHOD cnexp
        gka = gkabar*n*l
        ik = gka*(v-ek)
}

DERIVATIVE states {
        rates(v)
        n' = (ninf-n)/taun
        l' = (linf-l)/taul
}

FUNCTION alpn(v(mV)) {
LOCAL zeta
        zeta=zetan+pw/(1+exp((v-tq)/qq))
        alpn = exp(zeta*(v-vhalfn)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/mol/degC)*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) {
LOCAL zeta
        zeta=zetan+pw/(1+exp((v-tq)/qq))
        betn = exp(zeta*gmn*(v-vhalfn)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/mol/degC)*(273.16(degC)+celsius))) 
}


FUNCTION alpl(v(mV)) {
        alpl = exp(zetal*(v-vhalfl)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/mol/degC)*(273.16(degC)+celsius))) 
}

FUNCTION betl(v(mV)) {
        betl = exp(zetal*(v-vhalfl)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/mol/degC)*(273.16(degC)+celsius))) 
}


PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10(degC))
        a = alpn(v)
        ninf = 1/(1 + a)
        taun = betn(v)/(qt*a0n*(1+a))
        if (taun<nmin) {taun=nmin}
        taun=taun/nscale

        a = alpl(v)
        linf = 1/(1 + a)
        taul = 0.26(ms/mV)*(v+50)
        if (taul<lmin) {taul=lmin}
        taul=taul/lscale
}
