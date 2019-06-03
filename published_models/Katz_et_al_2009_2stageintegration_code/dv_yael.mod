COMMENT
	deltav, deltavmax and vrest
ENDCOMMENT

NEURON {
	SUFFIX dv
	RANGE vrest, deltav, vmax, vmaxt, vmin, vmint
	RANGE dvmax, dvmaxt, dvmin, dvmint
	RANGE vfall, vfallt, dvfall, dvfallt
}

ASSIGNED {
	v (millivolt)
	vrest (millivolt)
	deltav (millivolt)
	vmax (millivolt)
	vmaxt (ms)
	vmin (millivolt)
	vmint (ms)
	dvmax (millivolt)
	dvmaxt (ms)
	dvmin (millivolt)
	dvmint (ms)
	vfall (millivolt)
	vfallt (ms)
	dvfall (millivolt)
	dvfallt (ms)
}

INITIAL {
	vrest = v
	deltav = 0
	vmax = v
	vmaxt = 0
	vmin = v
	vmint = 0
	dvmax = 0
	dvmaxt = 0
	dvmin = 0
	dvmint = 0
	vfall = 0
	vfallt = 0
	dvfall = 0
	dvfallt = 0
}

BREAKPOINT {

	if (t<dt) {
		vrest=v
	}
	deltav=v-vrest
	if (v>vmax) {
		vmax=v
		vmaxt=t
	}
	if (v<vmin) {
		vmin=v
		vmint=t
	}
	if (deltav>dvmax) {
		dvmax=deltav
		dvmaxt=t
	}
	if (deltav<dvmin) {
		dvmin=deltav
		dvmint=t
	}
	if (t>vmaxt+1-dt && t<vmaxt+1+dt) {
		vfall=v
		vfallt=t
		dvfall=deltav
		dvfallt=t
	}
}
