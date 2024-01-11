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
static constexpr auto number_of_datum_variables = 3;
static constexpr auto number_of_floating_point_variables = 18;
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
 
#define nrn_init _nrn_init__Kir2_3
#define _nrn_initial _nrn_initial__Kir2_3
#define nrn_cur _nrn_cur__Kir2_3
#define _nrn_current _nrn_current__Kir2_3
#define nrn_jacob _nrn_jacob__Kir2_3
#define nrn_state _nrn_state__Kir2_3
#define _net_receive _net_receive__Kir2_3 
#define _f_rate _f_rate__Kir2_3 
#define rate rate__Kir2_3 
#define states states__Kir2_3 
 
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
#define Aalpha_d _ml->template fpfield<0>(_iml)
#define Aalpha_d_columnindex 0
#define Kalpha_d _ml->template fpfield<1>(_iml)
#define Kalpha_d_columnindex 1
#define V0alpha_d _ml->template fpfield<2>(_iml)
#define V0alpha_d_columnindex 2
#define Abeta_d _ml->template fpfield<3>(_iml)
#define Abeta_d_columnindex 3
#define Kbeta_d _ml->template fpfield<4>(_iml)
#define Kbeta_d_columnindex 4
#define V0beta_d _ml->template fpfield<5>(_iml)
#define V0beta_d_columnindex 5
#define gkbar _ml->template fpfield<6>(_iml)
#define gkbar_columnindex 6
#define ik _ml->template fpfield<7>(_iml)
#define ik_columnindex 7
#define d_inf _ml->template fpfield<8>(_iml)
#define d_inf_columnindex 8
#define tau_d _ml->template fpfield<9>(_iml)
#define tau_d_columnindex 9
#define g _ml->template fpfield<10>(_iml)
#define g_columnindex 10
#define alpha_d _ml->template fpfield<11>(_iml)
#define alpha_d_columnindex 11
#define beta_d _ml->template fpfield<12>(_iml)
#define beta_d_columnindex 12
#define d _ml->template fpfield<13>(_iml)
#define d_columnindex 13
#define ek _ml->template fpfield<14>(_iml)
#define ek_columnindex 14
#define Dd _ml->template fpfield<15>(_iml)
#define Dd_columnindex 15
#define v _ml->template fpfield<16>(_iml)
#define v_columnindex 16
#define _g _ml->template fpfield<17>(_iml)
#define _g_columnindex 17
#define _ion_ek *(_ml->dptr_field<0>(_iml))
#define _p_ion_ek static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_ik *(_ml->dptr_field<1>(_iml))
#define _p_ion_ik static_cast<neuron::container::data_handle<double>>(_ppvar[1])
#define _ion_dikdv *(_ml->dptr_field<2>(_iml))
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alp_d(void);
 static void _hoc_bet_d(void);
 static void _hoc_rate(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 _prop_id = _nrn_get_prop_id(_prop);
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_Kir2_3", _hoc_setdata},
 {"alp_d_Kir2_3", _hoc_alp_d},
 {"bet_d_Kir2_3", _hoc_bet_d},
 {"rate_Kir2_3", _hoc_rate},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_alp_d(Prop*);
 static double _npy_bet_d(Prop*);
 static double _npy_rate(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"alp_d", _npy_alp_d},
 {"bet_d", _npy_bet_d},
 {"rate", _npy_rate},
 {0, 0}
};
#define alp_d alp_d_Kir2_3
#define bet_d bet_d_Kir2_3
 extern double alp_d( _internalthreadargsprotocomma_ double );
 extern double bet_d( _internalthreadargsprotocomma_ double );
 
static void _check_rate(_internalthreadargsproto_); 
static void _check_table_thread(_threadargsprotocomma_ int _type, _nrn_model_sorted_token const& _sorted_token) {
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml, _type};
  {
    auto* const _ml = &_lmr;
   _check_rate(_threadargs_);
   }
}
 /* declare global and static user variables */
#define usetable usetable_Kir2_3
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {"usetable_Kir2_3", 0, 1},
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"Aalpha_d_Kir2_3", "/ms"},
 {"Kalpha_d_Kir2_3", "mV"},
 {"V0alpha_d_Kir2_3", "mV"},
 {"Abeta_d_Kir2_3", "/ms"},
 {"Kbeta_d_Kir2_3", "mV"},
 {"V0beta_d_Kir2_3", "mV"},
 {"gkbar_Kir2_3", "mho/cm2"},
 {"ik_Kir2_3", "mA/cm2"},
 {"tau_d_Kir2_3", "ms"},
 {"g_Kir2_3", "mho/cm2"},
 {"alpha_d_Kir2_3", "/ms"},
 {"beta_d_Kir2_3", "/ms"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double d0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"usetable_Kir2_3", &usetable_Kir2_3},
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
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[3].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Kir2_3",
 "Aalpha_d_Kir2_3",
 "Kalpha_d_Kir2_3",
 "V0alpha_d_Kir2_3",
 "Abeta_d_Kir2_3",
 "Kbeta_d_Kir2_3",
 "V0beta_d_Kir2_3",
 "gkbar_Kir2_3",
 0,
 "ik_Kir2_3",
 "d_inf_Kir2_3",
 "tau_d_Kir2_3",
 "g_Kir2_3",
 "alpha_d_Kir2_3",
 "beta_d_Kir2_3",
 0,
 "d_Kir2_3",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 18);
 	/*initialize range parameters*/
 	Aalpha_d = 0.13289;
 	Kalpha_d = -24.3902;
 	V0alpha_d = -83.94;
 	Abeta_d = 0.16994;
 	Kbeta_d = 35.714;
 	V0beta_d = -83.94;
 	gkbar = 0.0009;
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 18);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 0); /* ek */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ik */
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _Kir23_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 5);
  _extcall_thread.resize(4);
  _thread_mem_init(_extcall_thread.data());
 _mechtype = nrn_get_mechtype(_mechanism[1]);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"Aalpha_d"} /* 0 */,
                                       _nrn_mechanism_field<double>{"Kalpha_d"} /* 1 */,
                                       _nrn_mechanism_field<double>{"V0alpha_d"} /* 2 */,
                                       _nrn_mechanism_field<double>{"Abeta_d"} /* 3 */,
                                       _nrn_mechanism_field<double>{"Kbeta_d"} /* 4 */,
                                       _nrn_mechanism_field<double>{"V0beta_d"} /* 5 */,
                                       _nrn_mechanism_field<double>{"gkbar"} /* 6 */,
                                       _nrn_mechanism_field<double>{"ik"} /* 7 */,
                                       _nrn_mechanism_field<double>{"d_inf"} /* 8 */,
                                       _nrn_mechanism_field<double>{"tau_d"} /* 9 */,
                                       _nrn_mechanism_field<double>{"g"} /* 10 */,
                                       _nrn_mechanism_field<double>{"alpha_d"} /* 11 */,
                                       _nrn_mechanism_field<double>{"beta_d"} /* 12 */,
                                       _nrn_mechanism_field<double>{"d"} /* 13 */,
                                       _nrn_mechanism_field<double>{"ek"} /* 14 */,
                                       _nrn_mechanism_field<double>{"Dd"} /* 15 */,
                                       _nrn_mechanism_field<double>{"v"} /* 16 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 17 */,
                                       _nrn_mechanism_field<double*>{"_ion_ek", "k_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_ik", "k_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_dikdv", "k_ion"} /* 2 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 3 */);
  hoc_register_prop_size(_mechtype, 18, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Kir2_3 /mnt/mydata/cerebellum/mod/Kir23.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_d_inf;
 static double *_t_tau_d;
static int _reset;
static const char *modelname = "Cerebellum Granule Cell Model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rate(_internalthreadargsprotocomma_ double);
static int rate(_internalthreadargsprotocomma_ double);
 
#define _deriv1_advance _thread[0].literal_value<int>()
#define _dith1 1
#define _recurse _thread[2].literal_value<int>()
#define _newtonspace1 _thread[3].literal_value<NewtonSpace*>()
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static void _n_rate(_internalthreadargsprotocomma_ double _lv);
 static neuron::container::field_index _slist2[1];
  static neuron::container::field_index _slist1[1], _dlist1[1];
 static int states(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   rate ( _threadargscomma_ v ) ;
   Dd = ( d_inf - d ) / tau_d ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 rate ( _threadargscomma_ v ) ;
 Dd = Dd  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_d )) ;
  return 0;
}
 /*END CVODE*/
 
static int states (_internalthreadargsproto_) {
  int _reset=0;
  int error = 0;
 {
  auto* _savstate1 =_thread[_dith1].get<double*>();
  auto* _dlist2 = _thread[_dith1].get<double*>() + 1;
  int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 for(int _id=0; _id < 1; _id++) {
  _savstate1[_id] = _ml->data(_iml, _slist1[_id]);
}
 error = nrn_newton_thread(_newtonspace1, 1, _slist2, neuron::scopmath::row_view{_ml, _iml}, states, _dlist2, _ml, _iml, _ppvar, _thread, _nt);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   rate ( _threadargscomma_ v ) ;
   Dd = ( d_inf - d ) / tau_d ;
   {int _id; for(_id=0; _id < 1; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _ml->data(_iml, _dlist1[_id]) - (_ml->data(_iml, _slist1[_id]) - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _ml->data(_iml, _slist1[_id]) - _savstate1[_id];}}}
 } }
 return _reset;}
 
double alp_d ( _internalthreadargsprotocomma_ double _lv ) {
   double _lalp_d;
 double _lQ10 ;
 _lQ10 = pow( 3.0 , ( ( celsius - 20.0 ) / 10.0 ) ) ;
   _lalp_d = _lQ10 * Aalpha_d * exp ( ( _lv - V0alpha_d ) / Kalpha_d ) ;
   
return _lalp_d;
 }
 
static void _hoc_alp_d(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for alp_d_Kir2_3. Requires prior call to setdata_Kir2_3 and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  alp_d ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_alp_d(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  alp_d ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 
double bet_d ( _internalthreadargsprotocomma_ double _lv ) {
   double _lbet_d;
 double _lQ10 ;
 _lQ10 = pow( 3.0 , ( ( celsius - 20.0 ) / 10.0 ) ) ;
   _lbet_d = _lQ10 * Abeta_d * exp ( ( _lv - V0beta_d ) / Kbeta_d ) ;
   
return _lbet_d;
 }
 
static void _hoc_bet_d(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for bet_d_Kir2_3. Requires prior call to setdata_Kir2_3 and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  bet_d ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_bet_d(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  bet_d ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 static double _mfac_rate, _tmin_rate;
  static void _check_rate(_internalthreadargsproto_) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_Aalpha_d;
  static double _sav_Kalpha_d;
  static double _sav_V0alpha_d;
  static double _sav_Abeta_d;
  static double _sav_Kbeta_d;
  static double _sav_V0beta_d;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_Aalpha_d != Aalpha_d) { _maktable = 1;}
  if (_sav_Kalpha_d != Kalpha_d) { _maktable = 1;}
  if (_sav_V0alpha_d != V0alpha_d) { _maktable = 1;}
  if (_sav_Abeta_d != Abeta_d) { _maktable = 1;}
  if (_sav_Kbeta_d != Kbeta_d) { _maktable = 1;}
  if (_sav_V0beta_d != V0beta_d) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rate =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rate)/200.; _mfac_rate = 1./_dx;
   for (_i=0, _x=_tmin_rate; _i < 201; _x += _dx, _i++) {
    _f_rate(_threadargscomma_ _x);
    _t_d_inf[_i] = d_inf;
    _t_tau_d[_i] = tau_d;
   }
   _sav_Aalpha_d = Aalpha_d;
   _sav_Kalpha_d = Kalpha_d;
   _sav_V0alpha_d = V0alpha_d;
   _sav_Abeta_d = Abeta_d;
   _sav_Kbeta_d = Kbeta_d;
   _sav_V0beta_d = V0beta_d;
   _sav_celsius = celsius;
  }
 }

 static int rate(_internalthreadargsprotocomma_ double _lv) { 
#if 0
_check_rate(_threadargs_);
#endif
 _n_rate(_threadargscomma_ _lv);
 return 0;
 }

 static void _n_rate(_internalthreadargsprotocomma_ double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rate(_threadargscomma_ _lv); return; 
}
 _xi = _mfac_rate * (_lv - _tmin_rate);
 if (std::isnan(_xi)) {
  d_inf = _xi;
  tau_d = _xi;
  return;
 }
 if (_xi <= 0.) {
 d_inf = _t_d_inf[0];
 tau_d = _t_tau_d[0];
 return; }
 if (_xi >= 200.) {
 d_inf = _t_d_inf[200];
 tau_d = _t_tau_d[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 d_inf = _t_d_inf[_i] + _theta*(_t_d_inf[_i+1] - _t_d_inf[_i]);
 tau_d = _t_tau_d[_i] + _theta*(_t_tau_d[_i+1] - _t_tau_d[_i]);
 }

 
static int  _f_rate ( _internalthreadargsprotocomma_ double _lv ) {
   double _la_d , _lb_d ;
 _la_d = alp_d ( _threadargscomma_ _lv ) ;
   _lb_d = bet_d ( _threadargscomma_ _lv ) ;
   tau_d = 1.0 / ( _la_d + _lb_d ) ;
   d_inf = _la_d / ( _la_d + _lb_d ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for rate_Kir2_3. Requires prior call to setdata_Kir2_3 and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 
#if 1
 _check_rate(_threadargs_);
#endif
 _r = 1.;
 rate ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_rate(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 
#if 1
 _check_rate(_threadargs_);
#endif
 _r = 1.;
 rate ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
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
  ek = _ion_ek;
     _ode_spec1 (_threadargs_);
  }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 1; ++_i) {
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
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[_dith1] = new double[2]{};
   _newtonspace1 = nrn_cons_newtonspace(1);
 }
 
static void _thread_cleanup(Datum* _thread) {
   delete[] _thread[_dith1].get<double*>();
   nrn_destroy_newtonspace(_newtonspace1);
 }

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  d = d0;
 {
   rate ( _threadargscomma_ v ) ;
   d = d_inf ;
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

#if 0
 _check_rate(_threadargs_);
#endif
   _v = _vec_v[_ni[_iml]];
 v = _v;
  ek = _ion_ek;
 initmodel(_threadargs_);
 }
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   g = gkbar * d ;
   ik = g * ( v - ek ) ;
   alpha_d = alp_d ( _threadargscomma_ v ) ;
   beta_d = bet_d ( _threadargscomma_ v ) ;
   }
 _current += ik;

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
  ek = _ion_ek;
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ik += ik ;
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
_ni = _ml_arg->_nodeindices;
size_t _cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
for (size_t _iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
  ek = _ion_ek;
 {  _deriv1_advance = 1;
 derivimplicit_thread(1, _slist1, _dlist1, neuron::scopmath::row_view{_ml, _iml}, states, _ml, _iml, _ppvar, _thread, _nt);
_deriv1_advance = 0;
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 1; ++_i) {
      _ml->data(_iml, _slist1[_i]) += dt*_ml->data(_iml, _dlist1[_i]);
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {d_columnindex, 0};  _dlist1[0] = {Dd_columnindex, 0};
 _slist2[0] = {d_columnindex, 0};
   _t_d_inf = makevector(201*sizeof(double));
   _t_tau_d = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/Kir23.mod";
    const char* nmodl_file_text = 
  "TITLE Cerebellum Granule Cell Model\n"
  "\n"
  "COMMENT\n"
  "	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental\n"
  "	Evidence and Modeling of a Slow K+-Dependent Mechanism\n"
  "	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,\n"
  "	Andrea Fontana and Giovanni Naldi\n"
  "\n"
  "Suffix from Ubc_Kir to Kir2_3\n"
  "ENDCOMMENT\n"
  " \n"
  "NEURON { \n"
  "	SUFFIX Kir2_3\n"
  "	USEION k READ ek WRITE ik \n"
  "	RANGE gkbar, ik, g, alpha_d, beta_d, ek\n"
  "	RANGE Aalpha_d, Kalpha_d, V0alpha_d\n"
  "	RANGE Abeta_d, Kbeta_d, V0beta_d\n"
  "	RANGE d_inf, tau_d \n"
  "} \n"
  " \n"
  "UNITS { \n"
  "	(mA) = (milliamp) \n"
  "	(mV) = (millivolt) \n"
  "} \n"
  " \n"
  "PARAMETER { \n"
  "	Aalpha_d = 0.13289 (/ms)\n"
  "\n"
  "	\n"
  "	Kalpha_d = -24.3902 (mV)\n"
  "\n"
  "	V0alpha_d = -83.94 (mV)\n"
  "	Abeta_d = 0.16994 (/ms)\n"
  "\n"
  "	\n"
  "	Kbeta_d = 35.714 (mV)\n"
  "\n"
  "	V0beta_d = -83.94 (mV)\n"
  "	v (mV) \n"
  "	gkbar = 0.0009 (mho/cm2) \n"
  "	ek (mV) \n"
  "	celsius = 30 (degC) \n"
  "} \n"
  "\n"
  "STATE { \n"
  "	d \n"
  "} \n"
  "\n"
  "ASSIGNED { \n"
  "	ik (mA/cm2) \n"
  "	d_inf \n"
  "	tau_d (ms) \n"
  "	g (mho/cm2) \n"
  "	alpha_d (/ms) \n"
  "	beta_d (/ms) \n"
  "} \n"
  " \n"
  "INITIAL { \n"
  "	rate(v) \n"
  "	d = d_inf \n"
  "} \n"
  " \n"
  "BREAKPOINT { \n"
  "	SOLVE states METHOD derivimplicit\n"
  "	g = gkbar*d   \n"
  "	ik = g*(v - ek) \n"
  "	alpha_d = alp_d(v) \n"
  "	beta_d = bet_d(v) \n"
  "} \n"
  " \n"
  "DERIVATIVE states { \n"
  "	rate(v) \n"
  "	d' =(d_inf - d)/tau_d \n"
  "} \n"
  " \n"
  "FUNCTION alp_d(v(mV))(/ms) { LOCAL Q10\n"
  "	Q10 = 3^((celsius-20(degC))/10(degC))\n"
  "	alp_d = Q10*Aalpha_d*exp((v-V0alpha_d)/Kalpha_d) \n"
  "} \n"
  " \n"
  "FUNCTION bet_d(v(mV))(/ms) { LOCAL Q10\n"
  "	Q10 = 3^((celsius-20(degC))/10(degC))\n"
  "	bet_d = Q10*Abeta_d*exp((v-V0beta_d)/Kbeta_d) \n"
  "} \n"
  " \n"
  "PROCEDURE rate(v (mV)) {LOCAL a_d, b_d \n"
  "	TABLE d_inf, tau_d  \n"
  "	DEPEND Aalpha_d, Kalpha_d, V0alpha_d, \n"
  "	       Abeta_d, Kbeta_d, V0beta_d, celsius FROM -100 TO 100 WITH 200 \n"
  "	a_d = alp_d(v)  \n"
  "	b_d = bet_d(v) \n"
  "	tau_d = 1/(a_d + b_d) \n"
  "	d_inf = a_d/(a_d + b_d) \n"
  "} \n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
