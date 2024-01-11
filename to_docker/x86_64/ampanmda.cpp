/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format off
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
// clang-format on
#include "neuron/cache/mechanism_range.hpp"
static constexpr auto number_of_datum_variables = 2;
static constexpr auto number_of_floating_point_variables = 28;
namespace {
template <typename T>
using _nrn_mechanism_std_vector = std::vector<T>;
using _nrn_model_sorted_token = neuron::model_sorted_token;
using _nrn_mechanism_cache_range = neuron::cache::MechanismRange<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_mechanism_cache_instance = neuron::cache::MechanismInstance<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_non_owning_id_without_container = neuron::container::non_owning_identifier_without_container;
template <typename T>
using _nrn_mechanism_field = neuron::mechanism::field<T>;
template <typename... Args>
void _nrn_mechanism_register_data_fields(Args&&... args) {
  neuron::mechanism::register_data_fields(std::forward<Args>(args)...);
}
}
 
#if !NRNGPU
#undef exp
#define exp hoc_Exp
#endif
 
#define nrn_init _nrn_init__AmpaNmda
#define _nrn_initial _nrn_initial__AmpaNmda
#define nrn_cur _nrn_cur__AmpaNmda
#define _nrn_current _nrn_current__AmpaNmda
#define nrn_jacob _nrn_jacob__AmpaNmda
#define nrn_state _nrn_state__AmpaNmda
#define _net_receive _net_receive__AmpaNmda 
#define release release__AmpaNmda 
 
#define _threadargscomma_ _ml, _iml, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _internalthreadargsprotocomma_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _ml, _iml, _ppvar, _thread, _nt
#define _threadargsproto_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, NrnThread* _nt
#define _internalthreadargsproto_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t _nt->_t
#define dt _nt->_dt
#define mg _ml->template fpfield<0>(_iml)
#define mg_columnindex 0
#define gmax _ml->template fpfield<1>(_iml)
#define gmax_columnindex 1
#define ltdinvl _ml->template fpfield<2>(_iml)
#define ltdinvl_columnindex 2
#define ltpinvl _ml->template fpfield<3>(_iml)
#define ltpinvl_columnindex 3
#define sighalf _ml->template fpfield<4>(_iml)
#define sighalf_columnindex 4
#define sigslope _ml->template fpfield<5>(_iml)
#define sigslope_columnindex 5
#define sigexp _ml->template fpfield<6>(_iml)
#define sigexp_columnindex 6
#define x _ml->template fpfield<7>(_iml)
#define x_columnindex 7
#define mgid _ml->template fpfield<8>(_iml)
#define mgid_columnindex 8
#define ggid _ml->template fpfield<9>(_iml)
#define ggid_columnindex 9
#define srcgid _ml->template fpfield<10>(_iml)
#define srcgid_columnindex 10
#define training _ml->template fpfield<11>(_iml)
#define training_columnindex 11
#define i _ml->template fpfield<12>(_iml)
#define i_columnindex 12
#define inmda _ml->template fpfield<13>(_iml)
#define inmda_columnindex 13
#define iampa _ml->template fpfield<14>(_iml)
#define iampa_columnindex 14
#define gnmda _ml->template fpfield<15>(_iml)
#define gnmda_columnindex 15
#define Rinf _ml->template fpfield<16>(_iml)
#define Rinf_columnindex 16
#define Rtau _ml->template fpfield<17>(_iml)
#define Rtau_columnindex 17
#define Ron _ml->template fpfield<18>(_iml)
#define Ron_columnindex 18
#define Roff _ml->template fpfield<19>(_iml)
#define Roff_columnindex 19
#define gampa _ml->template fpfield<20>(_iml)
#define gampa_columnindex 20
#define synon _ml->template fpfield<21>(_iml)
#define synon_columnindex 21
#define DRon _ml->template fpfield<22>(_iml)
#define DRon_columnindex 22
#define DRoff _ml->template fpfield<23>(_iml)
#define DRoff_columnindex 23
#define Dgampa _ml->template fpfield<24>(_iml)
#define Dgampa_columnindex 24
#define v _ml->template fpfield<25>(_iml)
#define v_columnindex 25
#define _g _ml->template fpfield<26>(_iml)
#define _g_columnindex 26
#define _tsav _ml->template fpfield<27>(_iml)
#define _tsav_columnindex 27
#define _nd_area *_ml->dptr_field<0>(_iml)
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_mgblock(void*);
 static double _hoc_plast(void*);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {0, 0}
};
 static Member_func _member_func[] = {
 {"loc", _hoc_loc_pnt},
 {"has_loc", _hoc_has_loc},
 {"get_loc", _hoc_get_loc_pnt},
 {"mgblock", _hoc_mgblock},
 {"plast", _hoc_plast},
 {0, 0}
};
#define mgblock mgblock_AmpaNmda
#define plast plast_AmpaNmda
 extern double mgblock( _internalthreadargsprotocomma_ double );
 extern double plast( _internalthreadargsprotocomma_ double );
 /* declare global and static user variables */
#define Alpha Alpha_AmpaNmda
 double Alpha = 0.35;
#define Beta Beta_AmpaNmda
 double Beta = 0.012;
#define Cdur Cdur_AmpaNmda
 double Cdur = 1;
#define E E_AmpaNmda
 double E = 0;
#define ampatau ampatau_AmpaNmda
 double ampatau = 3;
#define gampafactor gampafactor_AmpaNmda
 double gampafactor = 0.001;
#define nmdafactor nmdafactor_AmpaNmda
 double nmdafactor = 0.07;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"Cdur_AmpaNmda", "ms"},
 {"Alpha_AmpaNmda", "/ms"},
 {"Beta_AmpaNmda", "/ms"},
 {"E_AmpaNmda", "mV"},
 {"gampafactor_AmpaNmda", "1"},
 {"nmdafactor_AmpaNmda", "1"},
 {"ampatau_AmpaNmda", "ms"},
 {"mg", "mM"},
 {"gmax", "umho"},
 {"ltdinvl", "ms"},
 {"ltpinvl", "ms"},
 {"sighalf", "1"},
 {"sigslope", "1"},
 {"x", "um"},
 {"gampa", "umho"},
 {"i", "nA"},
 {"inmda", "nA"},
 {"iampa", "nA"},
 {"gnmda", "umho"},
 {"Rtau", "ms"},
 {0, 0}
};
 static double Roff0 = 0;
 static double Ron0 = 0;
 static double delta_t = 0.01;
 static double gampa0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"Cdur_AmpaNmda", &Cdur_AmpaNmda},
 {"Alpha_AmpaNmda", &Alpha_AmpaNmda},
 {"Beta_AmpaNmda", &Beta_AmpaNmda},
 {"E_AmpaNmda", &E_AmpaNmda},
 {"gampafactor_AmpaNmda", &gampafactor_AmpaNmda},
 {"nmdafactor_AmpaNmda", &nmdafactor_AmpaNmda},
 {"ampatau_AmpaNmda", &ampatau_AmpaNmda},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void nrn_cur(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_jacob(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[3].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"AmpaNmda",
 "mg",
 "gmax",
 "ltdinvl",
 "ltpinvl",
 "sighalf",
 "sigslope",
 "sigexp",
 "x",
 "mgid",
 "ggid",
 "srcgid",
 "training",
 0,
 "i",
 "inmda",
 "iampa",
 "gnmda",
 "Rinf",
 "Rtau",
 0,
 "Ron",
 "Roff",
 "gampa",
 0,
 0};
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
  if (nrn_point_prop_) {
    _nrn_mechanism_access_alloc_seq(_prop) = _nrn_mechanism_access_alloc_seq(nrn_point_prop_);
    _ppvar = _nrn_mechanism_access_dparam(nrn_point_prop_);
  } else {
   _ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 28);
 	/*initialize range parameters*/
 	mg = 1;
 	gmax = 2;
 	ltdinvl = 250;
 	ltpinvl = 33.33;
 	sighalf = 25;
 	sigslope = 5;
 	sigexp = 4;
 	x = 0;
 	mgid = -1;
 	ggid = -1;
 	srcgid = -1;
 	training = 1;
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 28);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 
#define _tqitem &(_ppvar[2])
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _ampanmda_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"mg"} /* 0 */,
                                       _nrn_mechanism_field<double>{"gmax"} /* 1 */,
                                       _nrn_mechanism_field<double>{"ltdinvl"} /* 2 */,
                                       _nrn_mechanism_field<double>{"ltpinvl"} /* 3 */,
                                       _nrn_mechanism_field<double>{"sighalf"} /* 4 */,
                                       _nrn_mechanism_field<double>{"sigslope"} /* 5 */,
                                       _nrn_mechanism_field<double>{"sigexp"} /* 6 */,
                                       _nrn_mechanism_field<double>{"x"} /* 7 */,
                                       _nrn_mechanism_field<double>{"mgid"} /* 8 */,
                                       _nrn_mechanism_field<double>{"ggid"} /* 9 */,
                                       _nrn_mechanism_field<double>{"srcgid"} /* 10 */,
                                       _nrn_mechanism_field<double>{"training"} /* 11 */,
                                       _nrn_mechanism_field<double>{"i"} /* 12 */,
                                       _nrn_mechanism_field<double>{"inmda"} /* 13 */,
                                       _nrn_mechanism_field<double>{"iampa"} /* 14 */,
                                       _nrn_mechanism_field<double>{"gnmda"} /* 15 */,
                                       _nrn_mechanism_field<double>{"Rinf"} /* 16 */,
                                       _nrn_mechanism_field<double>{"Rtau"} /* 17 */,
                                       _nrn_mechanism_field<double>{"Ron"} /* 18 */,
                                       _nrn_mechanism_field<double>{"Roff"} /* 19 */,
                                       _nrn_mechanism_field<double>{"gampa"} /* 20 */,
                                       _nrn_mechanism_field<double>{"synon"} /* 21 */,
                                       _nrn_mechanism_field<double>{"DRon"} /* 22 */,
                                       _nrn_mechanism_field<double>{"DRoff"} /* 23 */,
                                       _nrn_mechanism_field<double>{"Dgampa"} /* 24 */,
                                       _nrn_mechanism_field<double>{"v"} /* 25 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 26 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 27 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<void*>{"_tqitem", "netsend"} /* 2 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 3 */);
  hoc_register_prop_size(_mechtype, 28, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 6;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 AmpaNmda /mnt/mydata/cerebellum/mod/ampanmda.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "simple NMDA receptors";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[3], _dlist1[3];
 static int release(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   DRon = ( synon * Rinf - Ron ) / Rtau ;
   DRoff = - Beta * Roff ;
   Dgampa = - gampa / ampatau ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 DRon = DRon  / (1. - dt*( ( ( ( - 1.0 ) ) ) / Rtau )) ;
 DRoff = DRoff  / (1. - dt*( ( - Beta )*( 1.0 ) )) ;
 Dgampa = Dgampa  / (1. - dt*( ( - 1.0 ) / ampatau )) ;
  return 0;
}
 /*END CVODE*/
 static int release (_internalthreadargsproto_) { {
    Ron = Ron + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / Rtau)))*(- ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - Ron) ;
    Roff = Roff + (1. - exp(dt*(( - Beta )*( 1.0 ))))*(- ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - Roff) ;
    gampa = gampa + (1. - exp(dt*(( - 1.0 ) / ampatau)))*(- ( 0.0 ) / ( ( - 1.0 ) / ampatau ) - gampa) ;
   }
  return 0;
}
 
double mgblock ( _internalthreadargsprotocomma_ double _lv ) {
   double _lmgblock;
 _lmgblock = 1.0 / ( 1.0 + exp ( 0.062 * - _lv ) * ( mg / 3.57 ) ) ;
   
return _lmgblock;
 }
 
static double _hoc_mgblock(void* _vptr) {
 double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _nrn_mechanism_cache_instance _ml_real{_p};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  _ppvar = _nrn_mechanism_access_dparam(_p);
  _thread = _extcall_thread.data();
  _nt = static_cast<NrnThread*>(_pnt->_vnt);
 _r =  mgblock ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 
double plast ( _internalthreadargsprotocomma_ double _lstep ) {
   double _lplast;
 _lplast = pow( ( 1.0 - 1.0 / ( 1.0 + exp ( ( _lstep - sighalf ) / sigslope ) ) ) , sigexp ) ;
   
return _lplast;
 }
 
static double _hoc_plast(void* _vptr) {
 double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _nrn_mechanism_cache_instance _ml_real{_p};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  _ppvar = _nrn_mechanism_access_dparam(_p);
  _thread = _extcall_thread.data();
  _nt = static_cast<NrnThread*>(_pnt->_vnt);
 _r =  plast ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  Prop* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
   _thread = nullptr; _nt = (NrnThread*)_pnt->_vnt;   _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = nullptr;}
 {
   if ( _lflag  == 0.0 ) {
     if ( t - _args[3] < ltpinvl ) {
       if ( training ) {
         _args[1] = _args[1] + 1.0 ;
         if ( _args[1] > 75.0 ) {
           _args[1] = 75.0 ;
           }
         }
       }
     else if ( t - _args[3] > ltdinvl ) {
       }
     else {
       if ( training ) {
         _args[1] = _args[1] - 1.0 ;
         if ( _args[1] < 0.0 ) {
           _args[1] = 0.0 ;
           }
         }
       }
     _args[3] = t ;
     if ( _args[1]  == 0.0 ) {
       _args[2] = 0.0 ;
       }
     else {
       _args[2] = _args[0] * plast ( _threadargscomma_ _args[1] ) ;
       }
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = gampa;
    double __primary = (gampa + _args[2] * gmax * gampafactor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / ampatau ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / ampatau ) - __primary );
    gampa += __primary;
  } else {
 gampa = gampa + _args[2] * gmax * gampafactor ;
       }
 _args[4] = _args[4] * exp ( - Beta * ( t - _args[5] ) ) ;
     _args[5] = t ;
     synon = synon + _args[2] ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron + _args[4]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / Rtau ) ) )*( - ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron + _args[4] ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff - _args[4]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - Beta )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff - _args[4] ;
       }
 net_send ( _tqitem, _args, _pnt, t +  Cdur , _args[2] + 1.0 ) ;
     }
   else {
     _args[4] = ( _lflag - 1.0 ) * Rinf + ( _args[4] - ( _lflag - 1.0 ) * Rinf ) * exp ( - ( t - _args[5] ) / Rtau ) ;
     _args[5] = t ;
     synon = synon - ( _lflag - 1.0 ) ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron - _args[4]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / Rtau ) ) )*( - ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron - _args[4] ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff + _args[4]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - Beta )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff + _args[4] ;
       }
 }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
     _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  Datum* _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  Datum* _thread = (Datum*)0;
  NrnThread* _nt = (NrnThread*)_pnt->_vnt;
 if ( _args[1]  == 0.0 ) {
     _args[2] = 0.0 ;
     }
   else {
     _args[2] = _args[0] * plast ( _threadargscomma_ _args[1] ) ;
     }
   _args[3] = - 1e9 ;
   _args[4] = 0.0 ;
   _args[5] = - 1e9 ;
   }
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_threadargs_);
 }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 3; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _ode_matsol1 (_threadargs_);
 }
 
static void _ode_matsol(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  Roff = Roff0;
  Ron = Ron0;
  gampa = gampa0;
 {
   Rinf = Alpha / ( Alpha + Beta ) ;
   Rtau = 1.0 / ( Alpha + Beta ) ;
   synon = 0.0 ;
   gampa = 0.0 ;
   }
 
}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _tsav = -1e20;
   _v = _vec_v[_ni[_iml]];
 v = _v;
 initmodel(_threadargs_);
}
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   gnmda = mgblock ( _threadargscomma_ v ) * ( Ron + Roff ) * gmax * nmdafactor ;
   inmda = gnmda * ( v - E ) ;
   iampa = gampa * ( v - E ) ;
   i = iampa + inmda ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_rhs = _nt->node_rhs_storage();
auto const _vec_sav_rhs = _nt->node_sav_rhs_storage();
auto const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ _rhs = _nrn_current(_threadargscomma_ _v);
 	}
 _g = (_g_local - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
	 _vec_rhs[_ni[_iml]] -= _rhs;
 
}
 
}

static void nrn_jacob(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_d = _nt->node_d_storage();
auto const _vec_sav_d = _nt->node_sav_d_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
  _vec_d[_ni[_iml]] += _g;
 
}
 
}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni;
_ni = _ml_arg->_nodeindices;
size_t _cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
for (size_t _iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
 {   release(_threadargs_);
  }}}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {Ron_columnindex, 0};  _dlist1[0] = {DRon_columnindex, 0};
 _slist1[1] = {Roff_columnindex, 0};  _dlist1[1] = {DRoff_columnindex, 0};
 _slist1[2] = {gampa_columnindex, 0};  _dlist1[2] = {Dgampa_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/ampanmda.mod";
    const char* nmodl_file_text = 
  "TITLE simple NMDA receptors\n"
  "\n"
  ": Hines combined AMPA and NMDA and spike dependent plasticity\n"
  "\n"
  ": Modified from the original AMPA.mod, M.Migliore Jan 2003\n"
  ": A weight of 0.0035 gives a peak conductance of 1nS in 0Mg\n"
  "\n"
  "COMMENT\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "	Simple model for glutamate AMPA receptors\n"
  "	=========================================\n"
  "\n"
  "  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS\n"
  "\n"
  "    Whole-cell recorded postsynaptic currents mediated by AMPA/Kainate\n"
  "    receptors (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994) were used\n"
  "    to estimate the parameters of the present model; the fit was performed\n"
  "    using a simplex algorithm (see Destexhe et al., J. Computational Neurosci.\n"
  "    1: 195-230, 1994).\n"
  "\n"
  "  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)\n"
  "\n"
  "    The simplified model was obtained from a detailed synaptic model that \n"
  "    included the release of transmitter in adjacent terminals, its lateral \n"
  "    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe\n"
  "    and Sejnowski, 1995).  Short pulses of transmitter with first-order\n"
  "    kinetics were found to be the best fast alternative to represent the more\n"
  "    detailed models.\n"
  "\n"
  "  - ANALYTIC EXPRESSION\n"
  "\n"
  "    The first-order model can be solved analytically, leading to a very fast\n"
  "    mechanism for simulating synapses, since no differential equation must be\n"
  "    solved (see references below).\n"
  "\n"
  "\n"
  "\n"
  "References\n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for\n"
  "   computing synaptic conductances based on a kinetic model of receptor binding\n"
  "   Neural Computation 6: 10-14, 1994.  \n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for\n"
  "   excitable membranes, synaptic transmission and neuromodulation using a \n"
  "   common kinetic formalism, Journal of Computational Neuroscience 1: \n"
  "   195-230, 1994.\n"
  "\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS AmpaNmda\n"
  "	RANGE R, g, mg, inmda, iampa, gnmda, gampa\n"
  "	RANGE x, mgid, ggid, srcgid, gmax, training\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL Cdur, Alpha, Beta, E, ampatau\n"
  "    RANGE Rinf, Rtau\n"
  "	GLOBAL gampafactor, nmdafactor\n"
  "	RANGE ltdinvl, ltpinvl, sighalf, sigslope, sigexp\n"
  "    THREADSAFE\n"
  "}\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "\n"
  "	Cdur	= 1		(ms)	: transmitter duration (rising phase)\n"
  "	Alpha	= 0.35		(/ms)	: forward (binding) rate\n"
  "	Beta	= 0.012		(/ms)	: backward (unbinding) rate\n"
  "	E	= 0	(mV)		: reversal potential\n"
  "	mg	= 1    (mM)		: external magnesium concentration\n"
  "	gmax = 2 (umho)		: normally 2\n"
  "	gampafactor = 0.001 (1)\n"
  "	nmdafactor = 0.07 (1)\n"
  "	ltdinvl = 250 (ms)		: longer intervals, no change\n"
  "	ltpinvl = 33.33 (ms)		: shorter interval, LTP\n"
  "	sighalf = 25 (1)\n"
  "	sigslope = 5 (1)\n"
  "        sigexp = 4\n"
  "	ampatau = 3 (ms)\n"
  "	x = 0 (um) : cartesian synapse location\n"
  "	mgid = -1 : associated mitral gid\n"
  "	ggid = -1 : associated granule gid\n"
  "	srcgid = -1 : gid of the mitral detector\n"
  "	training = 1 : is on\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)		: postsynaptic voltage\n"
  "	i 		(nA)		: total current = iampa+inmda\n"
  "	inmda 		(nA)		: current = gnmda*(v - E)\n"
  "	iampa 		(nA)		: current = gampa*(v - E)\n"
  "	gnmda 		(umho)		: \n"
  "	Rinf				: steady state channels open\n"
  "	Rtau		(ms)		: time constant of channel binding\n"
  "	synon\n"
  "}\n"
  "\n"
  "STATE {Ron Roff\n"
  "	gampa 		(umho)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	Rinf = Alpha / (Alpha + Beta)\n"
  "	Rtau = 1 / (Alpha + Beta)\n"
  "	synon = 0\n"
  "	gampa = 0\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE release METHOD cnexp\n"
  "	gnmda = mgblock(v)*(Ron + Roff)*gmax*nmdafactor\n"
  "	inmda = gnmda*(v - E)\n"
  "	iampa = gampa*(v - E)\n"
  "	i = iampa + inmda\n"
  "}\n"
  "\n"
  "DERIVATIVE release {\n"
  "	Ron' = (synon*Rinf - Ron)/Rtau\n"
  "	Roff' = -Beta*Roff\n"
  "	gampa' = -gampa/ampatau\n"
  "}\n"
  "\n"
  ": following supports both saturation from single input and\n"
  ": summation from multiple inputs\n"
  ": if spike occurs during CDur then new off time is t + CDur\n"
  ": ie. transmitter concatenates but does not summate\n"
  ": Note: automatic initialization of all reference args to 0 except first\n"
  "\n"
  "\n"
  "FUNCTION mgblock(v(mV)) {\n"
  "	:TABLE\n"
  "	:DEPEND mg\n"
  "	:FROM -140 TO 80 WITH 1000\n"
  "\n"
  "	: from Jahr & Stevens\n"
  "\n"
  "	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))\n"
  "}\n"
  "\n"
  "FUNCTION plast(step(1))(1) {\n"
  "	plast = (1 - 1/(1 + exp((step - sighalf)/sigslope)))^sigexp\n"
  "\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight, s, w, tlast (ms), r0, t0 (ms)) {\n"
  "	INITIAL {\n"
  "		:s = 0 \n"
  "                if(s == 0) {\n"
  "                  w = 0\n"
  "                } else {\n"
  "		  w = weight*plast(s)\n"
  "                }\n"
  "		tlast = -1e9 (ms)\n"
  "		r0 = 0\n"
  "		t0 = -1e9 (ms)\n"
  "	}\n"
  "	: flag is an implicit argument of NET_RECEIVE and  normally 0\n"
  "        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse\n"
  "		: plasticity affects this spike. If desired to affect\n"
  "		: the next spike then put following group after\n"
  "		: net_send\n"
  "		if (t - tlast < ltpinvl) { : LTP\n"
  "                        if (training) {\n"
  "			  s = s + 1\n"
  "			  :if (s > 2*sighalf) { s = 2*sighalf }\n"
  "                          if(s>75) {s=75}\n"
  "                        }\n"
  "		}else if (t - tlast > ltdinvl) { : no change\n"
  "		}else{ : LTD\n"
  "                      if (training) {\n"
  "			s = s - 1\n"
  "			if (s < 0) { s = 0 }\n"
  "                      }\n"
  "		}\n"
  "                         \n"
  "		tlast = t\n"
  "                         \n"
  "                if(s == 0) {\n"
  "                    w = 0\n"
  "                } else {\n"
  "                    w = weight*plast(s)\n"
  "                }\n"
  "                         \n"
  "		gampa = gampa + w*gmax*gampafactor\n"
  "		r0 = r0*exp(-Beta*(t - t0))\n"
  "		t0 = t\n"
  "		synon = synon + w\n"
  "		Ron = Ron + r0\n"
  "		Roff = Roff - r0\n"
  "		: come again in Cdur with flag = current value of w+1\n"
  "		net_send(Cdur, w + 1)\n"
  "        }else{ : turn off what was added Cdur ago\n"
  "		r0 = (flag-1)*Rinf + (r0 - (flag-1)*Rinf)*exp(-(t - t0)/Rtau)\n"
  "		t0 = t\n"
  "		synon = synon - (flag-1)\n"
  "		Ron = Ron - r0\n"
  "		Roff = Roff + r0\n"
  "	}\n"
  "}\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
