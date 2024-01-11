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
static constexpr auto number_of_floating_point_variables = 27;
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
 
#define nrn_init _nrn_init__FastInhibSTDP
#define _nrn_initial _nrn_initial__FastInhibSTDP
#define nrn_cur _nrn_cur__FastInhibSTDP
#define _nrn_current _nrn_current__FastInhibSTDP
#define nrn_jacob _nrn_jacob__FastInhibSTDP
#define nrn_state _nrn_state__FastInhibSTDP
#define _net_receive _net_receive__FastInhibSTDP 
#define state state__FastInhibSTDP 
 
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
#define tau1 _ml->template fpfield<0>(_iml)
#define tau1_columnindex 0
#define tau2 _ml->template fpfield<1>(_iml)
#define tau2_columnindex 1
#define gmax _ml->template fpfield<2>(_iml)
#define gmax_columnindex 2
#define e _ml->template fpfield<3>(_iml)
#define e_columnindex 3
#define wmax _ml->template fpfield<4>(_iml)
#define wmax_columnindex 4
#define aLTP _ml->template fpfield<5>(_iml)
#define aLTP_columnindex 5
#define aLTD _ml->template fpfield<6>(_iml)
#define aLTD_columnindex 6
#define mgid _ml->template fpfield<7>(_iml)
#define mgid_columnindex 7
#define ggid _ml->template fpfield<8>(_iml)
#define ggid_columnindex 8
#define srcgid _ml->template fpfield<9>(_iml)
#define srcgid_columnindex 9
#define i _ml->template fpfield<10>(_iml)
#define i_columnindex 10
#define interval _ml->template fpfield<11>(_iml)
#define interval_columnindex 11
#define tlast_pre _ml->template fpfield<12>(_iml)
#define tlast_pre_columnindex 12
#define tlast_post _ml->template fpfield<13>(_iml)
#define tlast_post_columnindex 13
#define M _ml->template fpfield<14>(_iml)
#define M_columnindex 14
#define P _ml->template fpfield<15>(_iml)
#define P_columnindex 15
#define deltaw _ml->template fpfield<16>(_iml)
#define deltaw_columnindex 16
#define wsyn _ml->template fpfield<17>(_iml)
#define wsyn_columnindex 17
#define A _ml->template fpfield<18>(_iml)
#define A_columnindex 18
#define B _ml->template fpfield<19>(_iml)
#define B_columnindex 19
#define g _ml->template fpfield<20>(_iml)
#define g_columnindex 20
#define factor _ml->template fpfield<21>(_iml)
#define factor_columnindex 21
#define DA _ml->template fpfield<22>(_iml)
#define DA_columnindex 22
#define DB _ml->template fpfield<23>(_iml)
#define DB_columnindex 23
#define v _ml->template fpfield<24>(_iml)
#define v_columnindex 24
#define _g _ml->template fpfield<25>(_iml)
#define _g_columnindex 25
#define _tsav _ml->template fpfield<26>(_iml)
#define _tsav_columnindex 26
#define _nd_area *_ml->dptr_field<0>(_iml)
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
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
 {0, 0}
};
 /* declare global and static user variables */
#define on on_FastInhibSTDP
 double on = 1;
#define tauLTD tauLTD_FastInhibSTDP
 double tauLTD = 20;
#define tauLTP tauLTP_FastInhibSTDP
 double tauLTP = 20;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {"tau2", 1e-09, 1e+09},
 {"tau1", 1e-09, 1e+09},
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"tauLTP_FastInhibSTDP", "ms"},
 {"tauLTD_FastInhibSTDP", "ms"},
 {"tau1", "ms"},
 {"tau2", "ms"},
 {"gmax", "uS"},
 {"e", "mV"},
 {"i", "nA"},
 {"interval", "ms"},
 {"tlast_pre", "ms"},
 {"tlast_post", "ms"},
 {0, 0}
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"tauLTP_FastInhibSTDP", &tauLTP_FastInhibSTDP},
 {"tauLTD_FastInhibSTDP", &tauLTD_FastInhibSTDP},
 {"on_FastInhibSTDP", &on_FastInhibSTDP},
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
 
#define _cvode_ieq _ppvar[2].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"FastInhibSTDP",
 "tau1",
 "tau2",
 "gmax",
 "e",
 "wmax",
 "aLTP",
 "aLTD",
 "mgid",
 "ggid",
 "srcgid",
 0,
 "i",
 "interval",
 "tlast_pre",
 "tlast_post",
 "M",
 "P",
 "deltaw",
 "wsyn",
 0,
 "A",
 "B",
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
   _ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 27);
 	/*initialize range parameters*/
 	tau1 = 1;
 	tau2 = 200;
 	gmax = 0.003;
 	e = -80;
 	wmax = 1;
 	aLTP = 0.001;
 	aLTD = 0.00106;
 	mgid = -1;
 	ggid = -1;
 	srcgid = -1;
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 27);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _fi_stdp_reg() {
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
                                       _nrn_mechanism_field<double>{"tau1"} /* 0 */,
                                       _nrn_mechanism_field<double>{"tau2"} /* 1 */,
                                       _nrn_mechanism_field<double>{"gmax"} /* 2 */,
                                       _nrn_mechanism_field<double>{"e"} /* 3 */,
                                       _nrn_mechanism_field<double>{"wmax"} /* 4 */,
                                       _nrn_mechanism_field<double>{"aLTP"} /* 5 */,
                                       _nrn_mechanism_field<double>{"aLTD"} /* 6 */,
                                       _nrn_mechanism_field<double>{"mgid"} /* 7 */,
                                       _nrn_mechanism_field<double>{"ggid"} /* 8 */,
                                       _nrn_mechanism_field<double>{"srcgid"} /* 9 */,
                                       _nrn_mechanism_field<double>{"i"} /* 10 */,
                                       _nrn_mechanism_field<double>{"interval"} /* 11 */,
                                       _nrn_mechanism_field<double>{"tlast_pre"} /* 12 */,
                                       _nrn_mechanism_field<double>{"tlast_post"} /* 13 */,
                                       _nrn_mechanism_field<double>{"M"} /* 14 */,
                                       _nrn_mechanism_field<double>{"P"} /* 15 */,
                                       _nrn_mechanism_field<double>{"deltaw"} /* 16 */,
                                       _nrn_mechanism_field<double>{"wsyn"} /* 17 */,
                                       _nrn_mechanism_field<double>{"A"} /* 18 */,
                                       _nrn_mechanism_field<double>{"B"} /* 19 */,
                                       _nrn_mechanism_field<double>{"g"} /* 20 */,
                                       _nrn_mechanism_field<double>{"factor"} /* 21 */,
                                       _nrn_mechanism_field<double>{"DA"} /* 22 */,
                                       _nrn_mechanism_field<double>{"DB"} /* 23 */,
                                       _nrn_mechanism_field<double>{"v"} /* 24 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 25 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 26 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 2 */);
  hoc_register_prop_size(_mechtype, 27, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 FastInhibSTDP /mnt/mydata/cerebellum/mod/fi_stdp.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[2], _dlist1[2];
 static int state(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state (_internalthreadargsproto_) { {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  Prop* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
   _thread = nullptr; _nt = (NrnThread*)_pnt->_vnt;   _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   if ( _args[0] >= 0.0 ) {
     P = P * exp ( ( tlast_pre - t ) / tauLTP ) + aLTP ;
     interval = tlast_post - t ;
     tlast_pre = t ;
     deltaw = wmax * M * exp ( interval / tauLTD ) ;
     }
   else {
     M = M * exp ( ( tlast_post - t ) / tauLTD ) - aLTD ;
     interval = t - tlast_pre ;
     tlast_post = t ;
     deltaw = wmax * P * exp ( - interval / tauLTP ) ;
     }
   if ( on ) {
     wsyn = wsyn + deltaw ;
     if ( wsyn > wmax ) {
       wsyn = wmax ;
       }
     if ( wsyn < 0.0 ) {
       wsyn = 0.0 ;
       }
     }
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + wsyn * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = A + wsyn * factor ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + wsyn * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + wsyn * factor ;
     }
 } }
 
static int _ode_count(int _type){ return 2;}
 
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
  for (int _i=0; _i < 2; ++_i) {
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
  A = A0;
  B = B0;
 {
   double _ltp ;
 if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   interval = 0.0 ;
   tlast_pre = 0.0 ;
   tlast_post = 0.0 ;
   M = 0.0 ;
   P = 0.0 ;
   deltaw = 0.0 ;
   wsyn = 0.0 ;
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
   g = ( B - A ) * gmax ;
   i = g * ( v - e ) ;
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
 {   state(_threadargs_);
  }}}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {A_columnindex, 0};  _dlist1[0] = {DA_columnindex, 0};
 _slist1[1] = {B_columnindex, 0};  _dlist1[1] = {DB_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/fi_stdp.mod";
    const char* nmodl_file_text = 
  ": Weight adjuster portion based on stdwa_songabbott.mod in ModeDB 64261\n"
  ": Conductance portion based on Exp2Syn.\n"
  ": This model is intended for use as the mitral side of the reciprocal\n"
  ": synapse and as such the pre events come from the granule side ThreshDetect\n"
  ": instance of the Mitral Granule Reciprocal Synapse (MGRS) with non-negative\n"
  ": weight (positive delay required),\n"
  ": and the post events come from the mitral side ThreshDetect instance\n"
  ": of the MGRS with negative weight (0 delay allowed).\n"
  "\n"
  "COMMENT\n"
  "Spike Timing Dependent Weight Adjuster\n"
  "based on Song and Abbott, 2001.\n"
  "Andrew Davison, UNIC, CNRS, 2003-2004\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS FastInhibSTDP\n"
  "\n"
  "	: conductance\n"
  "	RANGE tau1, tau2, e, i\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE gmax\n"
  "	RANGE mgid, ggid, srcgid\n"
  "\n"
  "	: weight adjuster\n"
  "	RANGE interval, tlast_pre, tlast_post, M, P\n"
  "	RANGE deltaw, wmax, aLTP, aLTD\n"
  "	RANGE wsyn\n"
  "	GLOBAL tauLTP, tauLTD, on\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	: conductance\n"
  "	tau1=1 (ms) <1e-9,1e9>\n"
  "	tau2 = 200 (ms) <1e-9,1e9>\n"
  "	gmax = .003 (uS) \n"
  "	e = -80	(mV)\n"
  "\n"
  "	: weight adjuster\n"
  "	tauLTP  = 20	(ms)    : decay time for LTP part ( values from           )\n"
  "	tauLTD  = 20	(ms)    : decay time for LTD part ( Song and Abbott, 2001 )\n"
  "	wmax    = 1		: min and max values of synaptic weight\n"
  "	aLTP    = 0.001		: amplitude of LTP steps\n"
  "	aLTD    = 0.00106	: amplitude of LTD steps\n"
  "	on	= 1		: allows learning to be turned on and off globally\n"
  "\n"
  "	: administrative\n"
  "	mgid = -1 : associated mitral gid\n"
  "	ggid = -1 : associated granule gid\n"
  "	srcgid = -1 : the gid of the granule detector\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	: conductance\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (uS)\n"
  "	factor\n"
  "\n"
  "	: weight adjuster\n"
  "	interval	(ms)	: since last spike of the other kind\n"
  "	tlast_pre	(ms)	: time of last presynaptic spike\n"
  "	tlast_post	(ms)	: time of last postsynaptic spike\n"
  "	M			: LTD function\n"
  "	P			: LTP function\n"
  "	deltaw			: change in weight\n"
  "	wsyn			: weight of the synapse\n"
  "\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A\n"
  "	B\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	: conductance\n"
  "	LOCAL tp\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "\n"
  "	: weight adjuster\n"
  "	interval = 0\n"
  "	tlast_pre = 0\n"
  "	tlast_post = 0\n"
  "	M = 0\n"
  "	P = 0\n"
  "	deltaw = 0\n"
  "	wsyn = 0\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	g = (B - A)*gmax\n"
  "	i = g*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  "NET_RECEIVE (w) {\n"
  "	if (w >= 0) {				: this is a pre-synaptic spike\n"
  "		P = P*exp((tlast_pre-t)/tauLTP) + aLTP\n"
  "		interval = tlast_post - t	: interval is negative\n"
  "		tlast_pre = t\n"
  "		deltaw = wmax * M * exp(interval/tauLTD)\n"
  "	} else {				: this is a post-synaptic spike\n"
  "		M = M*exp((tlast_post-t)/tauLTD) - aLTD\n"
  "		interval = t - tlast_pre	: interval is positive\n"
  "		tlast_post = t\n"
  "		deltaw = wmax * P * exp(-interval/tauLTP)\n"
  "	}\n"
  "	if (on) {\n"
  "		wsyn = wsyn + deltaw\n"
  "		if (wsyn > wmax) {\n"
  "			wsyn = wmax\n"
  "		}\n"
  "		if (wsyn < 0) {\n"
  "			wsyn = 0\n"
  "		}\n"
  "	}\n"
  "	A = A + wsyn*factor\n"
  "	B = B + wsyn*factor\n"
  "}\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
