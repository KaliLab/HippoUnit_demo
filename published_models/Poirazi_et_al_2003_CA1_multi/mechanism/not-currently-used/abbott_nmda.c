/* Created by Language version: 4.1  of 8/16/98 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#define exp hoc_Exp
extern double hoc_Exp();
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define delta_t dt
#define erev _p[0]
#define thresh _p[1]
#define d _p[2]
#define s _p[3]
#define tauD _p[4]
#define tauS _p[5]
#define taug _p[6]
#define G _p[7]
#define i _p[8]
#define D _p[9]
#define S _p[10]
#define g _p[11]
#define h _p[12]
#define DD _p[13]
#define DS _p[14]
#define Dg _p[15]
#define Dh _p[16]
#define firing _p[17]
#define _nd_area  *_ppvar[0].pval
#define vpre	*_ppvar[2].pval
#define _p_vpre	_ppvar[2].pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  2;
 /* external NEURON variables */
 extern double dt;
 extern double t;
 /* declaration of user functions */
 static double _hoc_aux();
 static double _hoc_check();
 static int _mechtype;
extern int nrn_get_mechtype();
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(_vptr) void* _vptr; { destroy_point_process(_vptr);}
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 static _hoc_setdata(_vptr) void* _vptr; { Prop* _prop;
 _prop = ((Point_process*)_vptr)->prop;
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 0,0
};
 static struct Member_func {
	char* _name; double (*_member)();} _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "aux", _hoc_aux,
 "check", _hoc_check,
 0, 0
};
 /* declare global and static user variables */
#define eta eta_abbott_nmda
 double eta = 0.33;
#define gamma gamma_abbott_nmda
 double gamma = 0.06;
#define mag mag_abbott_nmda
 double mag = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "eta_abbott_nmda", "/mM",
 "mag_abbott_nmda", "mM",
 "gamma_abbott_nmda", "/mV",
 "erev", "mV",
 "thresh", "mV",
 "tauD", "ms",
 "tauS", "ms",
 "taug", "ms",
 "G", "umho",
 "g", "umho",
 "h", "umho",
 "i", "nA",
 "vpre", "mV",
 0,0
};
 static double D0 = 0;
 static double S0 = 0;
 static double g0 = 0;
 static double h0 = 0;
 static double t0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "eta_abbott_nmda", &eta,
 "mag_abbott_nmda", &mag,
 "gamma_abbott_nmda", &gamma,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static nrn_alloc(), nrn_init(), nrn_cur(), nrn_state();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
extern int nrn_cvode_;
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "abbott_nmda",
 "erev",
 "thresh",
 "d",
 "s",
 "tauD",
 "tauS",
 "taug",
 "G",
 0,
 "i",
 0,
 "D",
 "S",
 "g",
 "h",
 0,
 "vpre",
 0};
 
static nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = (double *)ecalloc(18, sizeof(double));
 	/*initialize range parameters*/
 	erev = 0;
 	thresh = 0.5;
 	d = 0.4;
 	s = 1;
 	tauD = 300;
 	tauS = 20000;
 	taug = 2;
 	G = 1e-11;
  }
 	_prop->param = _p;
 	_prop->param_size = 18;
  if (!nrn_point_prop_) {
 	_ppvar = (Datum *)ecalloc(4, sizeof(Datum));
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 _abbott_nmda_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_state, nrn_init,
	 hoc_nrnpointerindex,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func,
	 _vectorized);
 _mechtype = nrn_get_mechtype(_mechanism[0]);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 abbott_nmda /disks/redondo/brannon/rs/bp/mechanism/abbott_nmda.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static aux();
static check();
 
static int _ode_spec1(), _ode_matsol1();
 static double *_temp1;
 static int _slist1[3], _dlist1[3];
 static int depression();
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   check (  ) ;
   DD = ( 1.0 / tauD ) * ( 1.0 - D ) ;
   DS = ( 1.0 / tauS ) * ( 1.0 - S ) ;
   Dg = ( 1.0 / taug ) * ( - g ) ;
   }
 return _reset;
}
 static int _ode_matsol1() {
 check (  ) ;
 DD = DD  / (1. - dt*( (( 1.0 / tauD ))*(( ( - 1.0 ) )) )) ;
 DS = DS  / (1. - dt*( (( 1.0 / tauS ))*(( ( - 1.0 ) )) )) ;
 Dg = Dg  / (1. - dt*( (( 1.0 / taug ))*(( - 1.0 )) )) ;
}
 /*END CVODE*/
 
static int depression () {_reset=0;
 {
   check (  ) ;
   DD = ( 1.0 / tauD ) * ( 1.0 - D ) ;
   DS = ( 1.0 / tauS ) * ( 1.0 - S ) ;
   Dg = ( 1.0 / taug ) * ( - g ) ;
   }
 return _reset;}
 
static int  aux (  )  {
   h = 1.0 / ( 1.0 + eta * mag * exp ( - ( gamma * v ) ) ) ;
   i = ( g * h * ( v - erev ) ) ;
    return 0; }
 static double _hoc_aux(_vptr) void* _vptr; {
 double _r;
 	_hoc_setdata(_vptr);
 _r = 1.;
 aux (  ) ;
 return(_r);
}
 
static int  check (  )  {
   if ( firing  && ( vpre < thresh ) ) {
     firing = 0.0 ;
     }
   if ( ( vpre >= thresh )  &&  ! firing ) {
     firing = 1.0 ;
     D = d * D ;
     S = s * S ;
     g = g + G * D * S ;
     }
    return 0; }
 static double _hoc_check(_vptr) void* _vptr; {
 double _r;
 	_hoc_setdata(_vptr);
 _r = 1.;
 check (  ) ;
 return(_r);
}
 
static int _ode_count() { return 3;}
 
static int _ode_spec(_nd, _pp, _ppd) Node* _nd; double* _pp; Datum* _ppd; {
	_p = _pp; _ppvar = _ppd; v = _nd->_v;
  _ode_spec1();
 }
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol) int _ieq; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; {
	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static int _ode_matsol(_nd, _pp, _ppd) Node* _nd; double* _pp; Datum* _ppd; {
	_p = _pp; _ppvar = _ppd; v = _nd->_v;
 _ode_matsol1();
 }

static initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = t0;
{
  D = D0;
  S = S0;
  g = g0;
  h = h0;
 {
   firing = 0.0 ;
   D = 1.0 ;
   S = 1.0 ;
   g = 0.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static nrn_init(_nd, _pp, _ppd) Node *_nd; double *_pp; Datum* _ppd; {
 double _v;
 _p = _pp; _ppvar = _ppd;
 _v = _nd->_v;
 v = _v;
 initmodel();
}

static double _nrn_current(_v) double _v;{double _current=0.;v=_v;{ {
   aux (  ) ;
   }
 _current += i;

} return _current;
}

static nrn_cur(_nd, _pp, _ppd) Node *_nd;double *_pp; Datum* _ppd;{
 double _g, _rhs, _v, *_pdiag;
 _p = _pp; _ppvar = _ppd;
_pdiag = &_nd->_d;
 _v = _nd->_v;
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
 _rhs -= _g*_v;
 *_pdiag += _g;
 _nd->_rhs -= _rhs;
 
}

static nrn_state(_nd, _pp, _ppd) Node *_nd; double *_pp; Datum* _ppd; {
 double _break, _save;
 double _v;
 _p = _pp; _ppvar = _ppd;
 _v = _nd->_v;
 _break = t + .5*dt; _save = t; delta_t = dt;
 v=_v;
{
 { {
 for (; t < _break; t += delta_t) {
 error =  euler(_ninits, 3, _slist1, _dlist1, _p, &t, delta_t, depression, &_temp1);
 if(error){fprintf(stderr,"at line 66 in file abbott_nmda.mod:\n   SOLVE depression METHOD euler : so now we're clinical psychiatrists\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
  depression();
 }
}
}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(D) - _p;  _dlist1[0] = &(DD) - _p;
 _slist1[1] = &(S) - _p;  _dlist1[1] = &(DS) - _p;
 _slist1[2] = &(g) - _p;  _dlist1[2] = &(Dg) - _p;
_first = 0;
}
