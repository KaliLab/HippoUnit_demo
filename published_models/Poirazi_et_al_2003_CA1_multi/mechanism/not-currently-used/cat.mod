TITLE transient and low threshold calcium current (T-current)

COMMENT
        *********************************************
        reference:      Huguenard & McCormick (1992) 
                        J.Neurophysiology 68(4), 1373-1383
        found in:       thalamic relay neurons
        *********************************************
        Assembled for MyFirstNEURON by Arthur Houweling
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX cat
        USEION ca READ cai,cao
        USEION Ca WRITE iCa VALENCE 2
        : The T-current does not activate calcium-dependent currents.
        : The construction with dummy ion Ca prevents the updating of the 
        : internal calcium concentration. 
        RANGE pcabar, m_inf, h_inf, tau_m, tau_h, iCa
}

UNITS {
        (mA)    = (milliamp)
        (mV)    = (millivolt)
        (mM)    = (milli/liter)
        FARADAY = 96480 (coul/mol)
        R       = 8.314 (volt-coul/degC)
}

PARAMETER {
        v               (mV)
        celsius         (degC)
        dt              (ms)
        cai = 5.e-05     (mM)
        cao = 2         (mM)
: maximum permiability!!! 
        pcabar= 0.0001  (cm/s)          
}

STATE {
        m h
}

ASSIGNED {
        iCa             (mA/cm2)
        tau_m           (ms)
        tau_h           (ms)
        m_inf 
        h_inf
        tadjm
        tadjh
}

BREAKPOINT { 
        SOLVE state :METHOD euler
        iCa = pcabar * m*m*h * ghk(v,cai,cao,2)
}

:DERIVATIVE state {
:       rates(v)
:   
:       m'= (m_inf-m) / tau_m
:       h'= (h_inf-h) / tau_h
:}
 
PROCEDURE state() {
        rates(v)

        m= m + (1-exp(-dt/tau_m))*(m_inf-m)
        h= h + (1-exp(-dt/tau_h))*(h_inf-h)
}

UNITSOFF
INITIAL {
        tadjm= 3.55^((celsius-23.5)/10)
        tadjh= 2.8^((celsius-23.5)/10)
        rates(v)
        m = m_inf
        h = h_inf
:        iCa = pcabar * m*m*h * ghk(v,cai,cao,2)
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else : denominator is small -> Taylor series
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}
UNITSOFF

PROCEDURE rates(v(mV)) { 
        tau_m = (1/(exp((v+131.6)/-16.7)+exp((v+16.8)/18.2)) + 0.612) / tadjm 
:        m_inf = 1 / (1+exp((v+60.5)/-6.2))
        m_inf = 1 / (1+exp((v+32)/-7.0))
        if (v<-80) 
          { tau_h = exp((v+467)/66.6) / tadjh }
        else  
          { tau_h = (exp((v+21.88)/-10.52)+28) / tadjh }
:        h_inf = 1 / (1+exp((v+84)/4.03)) 
         h_inf = 1 / (1+exp((v+67)/6.5)) 
}

UNITSON






