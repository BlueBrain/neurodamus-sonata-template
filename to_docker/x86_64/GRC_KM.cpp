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
static constexpr auto number_of_floating_point_variables = 20;
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
 
#define nrn_init _nrn_init__GRC_KM
#define _nrn_initial _nrn_initial__GRC_KM
#define nrn_cur _nrn_cur__GRC_KM
#define _nrn_current _nrn_current__GRC_KM
#define nrn_jacob _nrn_jacob__GRC_KM
#define nrn_state _nrn_state__GRC_KM
#define _net_receive _net_receive__GRC_KM 
#define _f_rate _f_rate__GRC_KM 
#define rate rate__GRC_KM 
#define states states__GRC_KM 
 
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
#define Aalpha_n _ml->template fpfield<0>(_iml)
#define Aalpha_n_columnindex 0
#define Kalpha_n _ml->template fpfield<1>(_iml)
#define Kalpha_n_columnindex 1
#define V0alpha_n _ml->template fpfield<2>(_iml)
#define V0alpha_n_columnindex 2
#define Abeta_n _ml->template fpfield<3>(_iml)
#define Abeta_n_columnindex 3
#define Kbeta_n _ml->template fpfield<4>(_iml)
#define Kbeta_n_columnindex 4
#define V0beta_n _ml->template fpfield<5>(_iml)
#define V0beta_n_columnindex 5
#define V0_ninf _ml->template fpfield<6>(_iml)
#define V0_ninf_columnindex 6
#define B_ninf _ml->template fpfield<7>(_iml)
#define B_ninf_columnindex 7
#define gkbar _ml->template fpfield<8>(_iml)
#define gkbar_columnindex 8
#define ik _ml->template fpfield<9>(_iml)
#define ik_columnindex 9
#define n_inf _ml->template fpfield<10>(_iml)
#define n_inf_columnindex 10
#define tau_n _ml->template fpfield<11>(_iml)
#define tau_n_columnindex 11
#define g _ml->template fpfield<12>(_iml)
#define g_columnindex 12
#define alpha_n _ml->template fpfield<13>(_iml)
#define alpha_n_columnindex 13
#define beta_n _ml->template fpfield<14>(_iml)
#define beta_n_columnindex 14
#define n _ml->template fpfield<15>(_iml)
#define n_columnindex 15
#define ek _ml->template fpfield<16>(_iml)
#define ek_columnindex 16
#define Dn _ml->template fpfield<17>(_iml)
#define Dn_columnindex 17
#define v _ml->template fpfield<18>(_iml)
#define v_columnindex 18
#define _g _ml->template fpfield<19>(_iml)
#define _g_columnindex 19
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
 static void _hoc_alp_n(void);
 static void _hoc_bet_n(void);
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
 {"setdata_GRC_KM", _hoc_setdata},
 {"alp_n_GRC_KM", _hoc_alp_n},
 {"bet_n_GRC_KM", _hoc_bet_n},
 {"rate_GRC_KM", _hoc_rate},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_alp_n(Prop*);
 static double _npy_bet_n(Prop*);
 static double _npy_rate(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"alp_n", _npy_alp_n},
 {"bet_n", _npy_bet_n},
 {"rate", _npy_rate},
 {0, 0}
};
#define alp_n alp_n_GRC_KM
#define bet_n bet_n_GRC_KM
 extern double alp_n( _internalthreadargsprotocomma_ double );
 extern double bet_n( _internalthreadargsprotocomma_ double );
 
static void _check_rate(_internalthreadargsproto_); 
static void _check_table_thread(_threadargsprotocomma_ int _type, _nrn_model_sorted_token const& _sorted_token) {
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml, _type};
  {
    auto* const _ml = &_lmr;
   _check_rate(_threadargs_);
   }
}
 /* declare global and static user variables */
#define usetable usetable_GRC_KM
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {"usetable_GRC_KM", 0, 1},
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"Aalpha_n_GRC_KM", "/ms"},
 {"Kalpha_n_GRC_KM", "mV"},
 {"V0alpha_n_GRC_KM", "mV"},
 {"Abeta_n_GRC_KM", "/ms"},
 {"Kbeta_n_GRC_KM", "mV"},
 {"V0beta_n_GRC_KM", "mV"},
 {"V0_ninf_GRC_KM", "mV"},
 {"B_ninf_GRC_KM", "mV"},
 {"gkbar_GRC_KM", "mho/cm2"},
 {"ik_GRC_KM", "mA/cm2"},
 {"tau_n_GRC_KM", "ms"},
 {"g_GRC_KM", "mho/cm2"},
 {"alpha_n_GRC_KM", "/ms"},
 {"beta_n_GRC_KM", "/ms"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"usetable_GRC_KM", &usetable_GRC_KM},
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
"GRC_KM",
 "Aalpha_n_GRC_KM",
 "Kalpha_n_GRC_KM",
 "V0alpha_n_GRC_KM",
 "Abeta_n_GRC_KM",
 "Kbeta_n_GRC_KM",
 "V0beta_n_GRC_KM",
 "V0_ninf_GRC_KM",
 "B_ninf_GRC_KM",
 "gkbar_GRC_KM",
 0,
 "ik_GRC_KM",
 "n_inf_GRC_KM",
 "tau_n_GRC_KM",
 "g_GRC_KM",
 "alpha_n_GRC_KM",
 "beta_n_GRC_KM",
 0,
 "n_GRC_KM",
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
    assert(_nrn_mechanism_get_num_vars(_prop) == 20);
 	/*initialize range parameters*/
 	Aalpha_n = 0.0033;
 	Kalpha_n = 40;
 	V0alpha_n = -30;
 	Abeta_n = 0.0033;
 	Kbeta_n = -20;
 	V0beta_n = -30;
 	V0_ninf = -35;
 	B_ninf = 6;
 	gkbar = 0.00025;
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 20);
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

 extern "C" void _GRC_KM_reg() {
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
                                       _nrn_mechanism_field<double>{"Aalpha_n"} /* 0 */,
                                       _nrn_mechanism_field<double>{"Kalpha_n"} /* 1 */,
                                       _nrn_mechanism_field<double>{"V0alpha_n"} /* 2 */,
                                       _nrn_mechanism_field<double>{"Abeta_n"} /* 3 */,
                                       _nrn_mechanism_field<double>{"Kbeta_n"} /* 4 */,
                                       _nrn_mechanism_field<double>{"V0beta_n"} /* 5 */,
                                       _nrn_mechanism_field<double>{"V0_ninf"} /* 6 */,
                                       _nrn_mechanism_field<double>{"B_ninf"} /* 7 */,
                                       _nrn_mechanism_field<double>{"gkbar"} /* 8 */,
                                       _nrn_mechanism_field<double>{"ik"} /* 9 */,
                                       _nrn_mechanism_field<double>{"n_inf"} /* 10 */,
                                       _nrn_mechanism_field<double>{"tau_n"} /* 11 */,
                                       _nrn_mechanism_field<double>{"g"} /* 12 */,
                                       _nrn_mechanism_field<double>{"alpha_n"} /* 13 */,
                                       _nrn_mechanism_field<double>{"beta_n"} /* 14 */,
                                       _nrn_mechanism_field<double>{"n"} /* 15 */,
                                       _nrn_mechanism_field<double>{"ek"} /* 16 */,
                                       _nrn_mechanism_field<double>{"Dn"} /* 17 */,
                                       _nrn_mechanism_field<double>{"v"} /* 18 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 19 */,
                                       _nrn_mechanism_field<double*>{"_ion_ek", "k_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_ik", "k_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_dikdv", "k_ion"} /* 2 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 3 */);
  hoc_register_prop_size(_mechtype, 20, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GRC_KM /mnt/mydata/cerebellum/mod/GRC_KM.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_n_inf;
 static double *_t_tau_n;
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
   Dn = ( n_inf - n ) / tau_n ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 rate ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_n )) ;
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
   Dn = ( n_inf - n ) / tau_n ;
   {int _id; for(_id=0; _id < 1; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _ml->data(_iml, _dlist1[_id]) - (_ml->data(_iml, _slist1[_id]) - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _ml->data(_iml, _slist1[_id]) - _savstate1[_id];}}}
 } }
 return _reset;}
 
double alp_n ( _internalthreadargsprotocomma_ double _lv ) {
   double _lalp_n;
 double _lQ10 ;
 _lQ10 = pow( 3.0 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   _lalp_n = _lQ10 * Aalpha_n * exp ( ( _lv - V0alpha_n ) / Kalpha_n ) ;
   
return _lalp_n;
 }
 
static void _hoc_alp_n(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for alp_n_GRC_KM. Requires prior call to setdata_GRC_KM and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  alp_n ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_alp_n(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  alp_n ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 
double bet_n ( _internalthreadargsprotocomma_ double _lv ) {
   double _lbet_n;
 double _lQ10 ;
 _lQ10 = pow( 3.0 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   _lbet_n = _lQ10 * Abeta_n * exp ( ( _lv - V0beta_n ) / Kbeta_n ) ;
   
return _lbet_n;
 }
 
static void _hoc_bet_n(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for bet_n_GRC_KM. Requires prior call to setdata_GRC_KM and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  bet_n ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_bet_n(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  bet_n ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 static double _mfac_rate, _tmin_rate;
  static void _check_rate(_internalthreadargsproto_) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_Aalpha_n;
  static double _sav_Kalpha_n;
  static double _sav_V0alpha_n;
  static double _sav_Abeta_n;
  static double _sav_Kbeta_n;
  static double _sav_V0beta_n;
  static double _sav_V0_ninf;
  static double _sav_B_ninf;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_Aalpha_n != Aalpha_n) { _maktable = 1;}
  if (_sav_Kalpha_n != Kalpha_n) { _maktable = 1;}
  if (_sav_V0alpha_n != V0alpha_n) { _maktable = 1;}
  if (_sav_Abeta_n != Abeta_n) { _maktable = 1;}
  if (_sav_Kbeta_n != Kbeta_n) { _maktable = 1;}
  if (_sav_V0beta_n != V0beta_n) { _maktable = 1;}
  if (_sav_V0_ninf != V0_ninf) { _maktable = 1;}
  if (_sav_B_ninf != B_ninf) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rate =  - 100.0 ;
   _tmax =  30.0 ;
   _dx = (_tmax - _tmin_rate)/13000.; _mfac_rate = 1./_dx;
   for (_i=0, _x=_tmin_rate; _i < 13001; _x += _dx, _i++) {
    _f_rate(_threadargscomma_ _x);
    _t_n_inf[_i] = n_inf;
    _t_tau_n[_i] = tau_n;
   }
   _sav_Aalpha_n = Aalpha_n;
   _sav_Kalpha_n = Kalpha_n;
   _sav_V0alpha_n = V0alpha_n;
   _sav_Abeta_n = Abeta_n;
   _sav_Kbeta_n = Kbeta_n;
   _sav_V0beta_n = V0beta_n;
   _sav_V0_ninf = V0_ninf;
   _sav_B_ninf = B_ninf;
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
  n_inf = _xi;
  tau_n = _xi;
  return;
 }
 if (_xi <= 0.) {
 n_inf = _t_n_inf[0];
 tau_n = _t_tau_n[0];
 return; }
 if (_xi >= 13000.) {
 n_inf = _t_n_inf[13000];
 tau_n = _t_tau_n[13000];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 n_inf = _t_n_inf[_i] + _theta*(_t_n_inf[_i+1] - _t_n_inf[_i]);
 tau_n = _t_tau_n[_i] + _theta*(_t_tau_n[_i+1] - _t_tau_n[_i]);
 }

 
static int  _f_rate ( _internalthreadargsprotocomma_ double _lv ) {
   double _la_n , _lb_n ;
 _la_n = alp_n ( _threadargscomma_ _lv ) ;
   _lb_n = bet_n ( _threadargscomma_ _lv ) ;
   tau_n = 1.0 / ( _la_n + _lb_n ) ;
   n_inf = 1.0 / ( 1.0 + exp ( - ( _lv - V0_ninf ) / B_ninf ) ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for rate_GRC_KM. Requires prior call to setdata_GRC_KM and that the specified mechanism instance still be in existence.", NULL);
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
  n = n0;
 {
   rate ( _threadargscomma_ v ) ;
   n = n_inf ;
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
   g = gkbar * n ;
   ik = g * ( v - ek ) ;
   alpha_n = alp_n ( _threadargscomma_ v ) ;
   beta_n = bet_n ( _threadargscomma_ v ) ;
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
 _slist1[0] = {n_columnindex, 0};  _dlist1[0] = {Dn_columnindex, 0};
 _slist2[0] = {n_columnindex, 0};
   _t_n_inf = makevector(13001*sizeof(double));
   _t_tau_n = makevector(13001*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/GRC_KM.mod";
    const char* nmodl_file_text = 
  "TITLE Cerebellum Granule Cell Model\n"
  "\n"
  "COMMENT\n"
  "        KM channel\n"
  "   \n"
  "	Author: A. Fontana\n"
  "	CoAuthor: T.Nieus Last revised: 20.11.99\n"
  "	\n"
  "ENDCOMMENT\n"
  " \n"
  "NEURON { \n"
  "	SUFFIX GRC_KM \n"
  "	USEION k READ ek WRITE ik \n"
  "	RANGE gkbar, ik, g, alpha_n, beta_n \n"
  "	RANGE Aalpha_n, Kalpha_n, V0alpha_n\n"
  "	RANGE Abeta_n, Kbeta_n, V0beta_n\n"
  "	RANGE V0_ninf, B_ninf\n"
  "	RANGE n_inf, tau_n \n"
  "} \n"
  " \n"
  "UNITS { \n"
  "	(mA) = (milliamp) \n"
  "	(mV) = (millivolt) \n"
  "} \n"
  " \n"
  "PARAMETER { \n"
  "	Aalpha_n = 0.0033 (/ms)\n"
  "	Kalpha_n = 40 (mV)\n"
  "\n"
  "	V0alpha_n = -30 (mV)\n"
  "	Abeta_n = 0.0033 (/ms)\n"
  "	Kbeta_n = -20 (mV)\n"
  "\n"
  "	V0beta_n = -30 (mV)\n"
  "	V0_ninf = -35 (mV)	:-30\n"
  "	B_ninf = 6 (mV)		:6:4 rimesso a 6 dopo calibrazione febbraio 2003	\n"
  "	v (mV) \n"
  "	gkbar= 0.00025 (mho/cm2) :0.0001\n"
  "	ek = -84.69 (mV) \n"
  "	celsius = 30 (degC) \n"
  "} \n"
  "\n"
  "STATE { \n"
  "	n \n"
  "} \n"
  "\n"
  "ASSIGNED { \n"
  "	ik (mA/cm2) \n"
  "	n_inf \n"
  "	tau_n (ms) \n"
  "	g (mho/cm2) \n"
  "	alpha_n (/ms) \n"
  "	beta_n (/ms) \n"
  "} \n"
  " \n"
  "INITIAL { \n"
  "	rate(v) \n"
  "	n = n_inf \n"
  "} \n"
  " \n"
  "BREAKPOINT { \n"
  "	SOLVE states METHOD derivimplicit \n"
  "	g = gkbar*n \n"
  "	ik = g*(v - ek) \n"
  "	alpha_n = alp_n(v) \n"
  "	beta_n = bet_n(v) \n"
  "} \n"
  " \n"
  "DERIVATIVE states { \n"
  "	rate(v) \n"
  "	n' =(n_inf - n)/tau_n \n"
  "} \n"
  " \n"
  "FUNCTION alp_n(v(mV))(/ms) { LOCAL Q10\n"
  "	Q10 = 3^((celsius-22(degC))/10(degC)) \n"
  "	alp_n = Q10*Aalpha_n*exp((v-V0alpha_n)/Kalpha_n) \n"
  "} \n"
  " \n"
  "FUNCTION bet_n(v(mV))(/ms) { LOCAL Q10\n"
  "	Q10 = 3^((celsius-22(degC))/10(degC)) \n"
  "	bet_n = Q10*Abeta_n*exp((v-V0beta_n)/Kbeta_n) \n"
  "} \n"
  " \n"
  "PROCEDURE rate(v (mV)) {LOCAL a_n, b_n \n"
  "	TABLE n_inf, tau_n \n"
  "	DEPEND Aalpha_n, Kalpha_n, V0alpha_n, \n"
  "	       Abeta_n, Kbeta_n, V0beta_n, V0_ninf, B_ninf, celsius FROM -100 TO 30 WITH 13000 \n"
  "	a_n = alp_n(v)  \n"
  "	b_n = bet_n(v) \n"
  "	tau_n = 1/(a_n + b_n) \n"
  ":	n_inf = a_n/(a_n + b_n) \n"
  "	n_inf = 1/(1+exp(-(v-V0_ninf)/B_ninf))\n"
  "} \n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
