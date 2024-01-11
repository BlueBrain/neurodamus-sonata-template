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
static constexpr auto number_of_datum_variables = 9;
static constexpr auto number_of_floating_point_variables = 26;
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
 
#define nrn_init _nrn_init__Kv1_5
#define _nrn_initial _nrn_initial__Kv1_5
#define nrn_cur _nrn_cur__Kv1_5
#define _nrn_current _nrn_current__Kv1_5
#define nrn_jacob _nrn_jacob__Kv1_5
#define nrn_state _nrn_state__Kv1_5
#define _net_receive _net_receive__Kv1_5 
#define rates rates__Kv1_5 
#define states states__Kv1_5 
 
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
#define gKur _ml->template fpfield<0>(_iml)
#define gKur_columnindex 0
#define Tauact _ml->template fpfield<1>(_iml)
#define Tauact_columnindex 1
#define Tauinactf _ml->template fpfield<2>(_iml)
#define Tauinactf_columnindex 2
#define Tauinacts _ml->template fpfield<3>(_iml)
#define Tauinacts_columnindex 3
#define gnonspec _ml->template fpfield<4>(_iml)
#define gnonspec_columnindex 4
#define ik _ml->template fpfield<5>(_iml)
#define ik_columnindex 5
#define minf _ml->template fpfield<6>(_iml)
#define minf_columnindex 6
#define ninf _ml->template fpfield<7>(_iml)
#define ninf_columnindex 7
#define uinf _ml->template fpfield<8>(_iml)
#define uinf_columnindex 8
#define mtau _ml->template fpfield<9>(_iml)
#define mtau_columnindex 9
#define ntau _ml->template fpfield<10>(_iml)
#define ntau_columnindex 10
#define utau _ml->template fpfield<11>(_iml)
#define utau_columnindex 11
#define ino _ml->template fpfield<12>(_iml)
#define ino_columnindex 12
#define m _ml->template fpfield<13>(_iml)
#define m_columnindex 13
#define n _ml->template fpfield<14>(_iml)
#define n_columnindex 14
#define u _ml->template fpfield<15>(_iml)
#define u_columnindex 15
#define Dm _ml->template fpfield<16>(_iml)
#define Dm_columnindex 16
#define Dn _ml->template fpfield<17>(_iml)
#define Dn_columnindex 17
#define Du _ml->template fpfield<18>(_iml)
#define Du_columnindex 18
#define ek _ml->template fpfield<19>(_iml)
#define ek_columnindex 19
#define ki _ml->template fpfield<20>(_iml)
#define ki_columnindex 20
#define ko _ml->template fpfield<21>(_iml)
#define ko_columnindex 21
#define nai _ml->template fpfield<22>(_iml)
#define nai_columnindex 22
#define nao _ml->template fpfield<23>(_iml)
#define nao_columnindex 23
#define v _ml->template fpfield<24>(_iml)
#define v_columnindex 24
#define _g _ml->template fpfield<25>(_iml)
#define _g_columnindex 25
#define _ion_ek *(_ml->dptr_field<0>(_iml))
#define _p_ion_ek static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_ki *(_ml->dptr_field<1>(_iml))
#define _p_ion_ki static_cast<neuron::container::data_handle<double>>(_ppvar[1])
#define _ion_ko *(_ml->dptr_field<2>(_iml))
#define _p_ion_ko static_cast<neuron::container::data_handle<double>>(_ppvar[2])
#define _ion_ik *(_ml->dptr_field<3>(_iml))
#define _p_ion_ik static_cast<neuron::container::data_handle<double>>(_ppvar[3])
#define _ion_dikdv *(_ml->dptr_field<4>(_iml))
#define _ion_nai *(_ml->dptr_field<5>(_iml))
#define _p_ion_nai static_cast<neuron::container::data_handle<double>>(_ppvar[5])
#define _ion_nao *(_ml->dptr_field<6>(_iml))
#define _p_ion_nao static_cast<neuron::container::data_handle<double>>(_ppvar[6])
#define _ion_ino *(_ml->dptr_field<7>(_iml))
#define _p_ion_ino static_cast<neuron::container::data_handle<double>>(_ppvar[7])
#define _ion_dinodv *(_ml->dptr_field<8>(_iml))
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alp(void);
 static void _hoc_bet(void);
 static void _hoc_ce(void);
 static void _hoc_rates(void);
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
 {"setdata_Kv1_5", _hoc_setdata},
 {"alp_Kv1_5", _hoc_alp},
 {"bet_Kv1_5", _hoc_bet},
 {"ce_Kv1_5", _hoc_ce},
 {"rates_Kv1_5", _hoc_rates},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_alp(Prop*);
 static double _npy_bet(Prop*);
 static double _npy_ce(Prop*);
 static double _npy_rates(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"alp", _npy_alp},
 {"bet", _npy_bet},
 {"ce", _npy_ce},
 {"rates", _npy_rates},
 {0, 0}
};
#define alp alp_Kv1_5
#define bet bet_Kv1_5
#define ce ce_Kv1_5
 extern double alp( _internalthreadargsprotocomma_ double , double );
 extern double bet( _internalthreadargsprotocomma_ double , double );
 extern double ce( _internalthreadargsprotocomma_ double , double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {"gnonspec_Kv1_5", 0, 1e+09},
 {"gKur_Kv1_5", 0, 1e+09},
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"gKur_Kv1_5", "S/cm2"},
 {"Tauact_Kv1_5", "ms"},
 {"Tauinactf_Kv1_5", "ms"},
 {"Tauinacts_Kv1_5", "ms"},
 {"gnonspec_Kv1_5", "S/cm2"},
 {"ik_Kv1_5", "mA/cm2"},
 {"mtau_Kv1_5", "ms"},
 {"ntau_Kv1_5", "ms"},
 {"utau_Kv1_5", "ms"},
 {"ino_Kv1_5", "mA/cm2"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double m0 = 0;
 static double n0 = 0;
 static double u0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 
#define _cvode_ieq _ppvar[9].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Kv1_5",
 "gKur_Kv1_5",
 "Tauact_Kv1_5",
 "Tauinactf_Kv1_5",
 "Tauinacts_Kv1_5",
 "gnonspec_Kv1_5",
 0,
 "ik_Kv1_5",
 "minf_Kv1_5",
 "ninf_Kv1_5",
 "uinf_Kv1_5",
 "mtau_Kv1_5",
 "ntau_Kv1_5",
 "utau_Kv1_5",
 "ino_Kv1_5",
 0,
 "m_Kv1_5",
 "n_Kv1_5",
 "u_Kv1_5",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _na_sym;
 static Symbol* _no_sym;
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 10, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 26);
 	/*initialize range parameters*/
 	gKur = 0.00013195;
 	Tauact = 1;
 	Tauinactf = 1;
 	Tauinacts = 1;
 	gnonspec = 0;
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 26);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 1, 1);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 0); /* ek */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 1); /* ki */
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 2); /* ko */
 	_ppvar[3] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ik */
 	_ppvar[4] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dikdv */
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[5] = _nrn_mechanism_get_param_handle(prop_ion, 1); /* nai */
 	_ppvar[6] = _nrn_mechanism_get_param_handle(prop_ion, 2); /* nao */
 prop_ion = need_memb(_no_sym);
 	_ppvar[7] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ino */
 	_ppvar[8] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dinodv */
 
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

 extern "C" void _Kv15_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("na", -10000.);
 	ion_reg("no", 1.0);
 	_k_sym = hoc_lookup("k_ion");
 	_na_sym = hoc_lookup("na_ion");
 	_no_sym = hoc_lookup("no_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 5);
  _extcall_thread.resize(4);
  _thread_mem_init(_extcall_thread.data());
 _mechtype = nrn_get_mechtype(_mechanism[1]);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"gKur"} /* 0 */,
                                       _nrn_mechanism_field<double>{"Tauact"} /* 1 */,
                                       _nrn_mechanism_field<double>{"Tauinactf"} /* 2 */,
                                       _nrn_mechanism_field<double>{"Tauinacts"} /* 3 */,
                                       _nrn_mechanism_field<double>{"gnonspec"} /* 4 */,
                                       _nrn_mechanism_field<double>{"ik"} /* 5 */,
                                       _nrn_mechanism_field<double>{"minf"} /* 6 */,
                                       _nrn_mechanism_field<double>{"ninf"} /* 7 */,
                                       _nrn_mechanism_field<double>{"uinf"} /* 8 */,
                                       _nrn_mechanism_field<double>{"mtau"} /* 9 */,
                                       _nrn_mechanism_field<double>{"ntau"} /* 10 */,
                                       _nrn_mechanism_field<double>{"utau"} /* 11 */,
                                       _nrn_mechanism_field<double>{"ino"} /* 12 */,
                                       _nrn_mechanism_field<double>{"m"} /* 13 */,
                                       _nrn_mechanism_field<double>{"n"} /* 14 */,
                                       _nrn_mechanism_field<double>{"u"} /* 15 */,
                                       _nrn_mechanism_field<double>{"Dm"} /* 16 */,
                                       _nrn_mechanism_field<double>{"Dn"} /* 17 */,
                                       _nrn_mechanism_field<double>{"Du"} /* 18 */,
                                       _nrn_mechanism_field<double>{"ek"} /* 19 */,
                                       _nrn_mechanism_field<double>{"ki"} /* 20 */,
                                       _nrn_mechanism_field<double>{"ko"} /* 21 */,
                                       _nrn_mechanism_field<double>{"nai"} /* 22 */,
                                       _nrn_mechanism_field<double>{"nao"} /* 23 */,
                                       _nrn_mechanism_field<double>{"v"} /* 24 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 25 */,
                                       _nrn_mechanism_field<double*>{"_ion_ek", "k_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_ki", "k_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_ko", "k_ion"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"_ion_ik", "k_ion"} /* 3 */,
                                       _nrn_mechanism_field<double*>{"_ion_dikdv", "k_ion"} /* 4 */,
                                       _nrn_mechanism_field<double*>{"_ion_nai", "na_ion"} /* 5 */,
                                       _nrn_mechanism_field<double*>{"_ion_nao", "na_ion"} /* 6 */,
                                       _nrn_mechanism_field<double*>{"_ion_ino", "no_ion"} /* 7 */,
                                       _nrn_mechanism_field<double*>{"_ion_dinodv", "no_ion"} /* 8 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 9 */);
  hoc_register_prop_size(_mechtype, 26, 10);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "no_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "no_ion");
  hoc_register_dparam_semantics(_mechtype, 9, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Kv1_5 /mnt/mydata/cerebellum/mod/Kv15.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double F = 0x1.78e555060882cp+16;
 static double R = 0x1.0a1013e8990bep+3;
static int _reset;
static const char *modelname = "Cardiac IKur  current & nonspec cation current with identical kinetics";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_internalthreadargsprotocomma_ double);
 
#define _deriv1_advance _thread[0].literal_value<int>()
#define _dith1 1
#define _recurse _thread[2].literal_value<int>()
#define _newtonspace1 _thread[3].literal_value<NewtonSpace*>()
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist2[3];
  static neuron::container::field_index _slist1[3], _dlist1[3];
 static int states(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dn = ( ninf - n ) / ntau ;
   Du = ( uinf - u ) / utau ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
 Du = Du  / (1. - dt*( ( ( ( - 1.0 ) ) ) / utau )) ;
  return 0;
}
 /*END CVODE*/
 
static int states (_internalthreadargsproto_) {
  int _reset=0;
  int error = 0;
 {
  auto* _savstate1 =_thread[_dith1].get<double*>();
  auto* _dlist2 = _thread[_dith1].get<double*>() + 3;
  int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 for(int _id=0; _id < 3; _id++) {
  _savstate1[_id] = _ml->data(_iml, _slist1[_id]);
}
 error = nrn_newton_thread(_newtonspace1, 3, _slist2, neuron::scopmath::row_view{_ml, _iml}, states, _dlist2, _ml, _iml, _ppvar, _thread, _nt);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dn = ( ninf - n ) / ntau ;
   Du = ( uinf - u ) / utau ;
   {int _id; for(_id=0; _id < 3; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _ml->data(_iml, _dlist1[_id]) - (_ml->data(_iml, _slist1[_id]) - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _ml->data(_iml, _slist1[_id]) - _savstate1[_id];}}}
 } }
 return _reset;}
 
double alp ( _internalthreadargsprotocomma_ double _lv , double _li ) {
   double _lalp;
 double _lq10 ;
 _lv = _lv ;
   _lq10 = pow( 2.2 , ( ( celsius - 37.0 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lalp = _lq10 * 0.65 / ( exp ( - ( _lv + 10.0 ) / 8.5 ) + exp ( - ( _lv - 30.0 ) / 59.0 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = 0.001 * _lq10 / ( 2.4 + 10.9 * exp ( - ( _lv + 90.0 ) / 78.0 ) ) ;
     }
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  Prop* _local_prop = _prop_id ? _extcall_prop : nullptr;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  alp ( _threadargscomma_ *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static double _npy_alp(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  alp ( _threadargscomma_ *getarg(1) , *getarg(2) );
 return(_r);
}
 
double bet ( _internalthreadargsprotocomma_ double _lv , double _li ) {
   double _lbet;
 double _lq10 ;
 _lv = _lv ;
   _lq10 = pow( 2.2 , ( ( celsius - 37.0 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lbet = _lq10 * 0.65 / ( 2.5 + exp ( ( _lv + 82.0 ) / 17.0 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = _lq10 * 0.001 * exp ( ( _lv - 168.0 ) / 16.0 ) ;
     }
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  Prop* _local_prop = _prop_id ? _extcall_prop : nullptr;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  bet ( _threadargscomma_ *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static double _npy_bet(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  bet ( _threadargscomma_ *getarg(1) , *getarg(2) );
 return(_r);
}
 
double ce ( _internalthreadargsprotocomma_ double _lv , double _li ) {
   double _lce;
 _lv = _lv ;
   if ( _li  == 0.0 ) {
     _lce = 1.0 / ( 1.0 + exp ( - ( _lv + 30.3 ) / 9.6 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lce = 1.0 * ( 0.25 + 1.0 / ( 1.35 + exp ( ( _lv + 7.0 ) / 14.0 ) ) ) ;
     }
   else if ( _li  == 2.0 ) {
     _lce = 1.0 * ( 0.1 + 1.0 / ( 1.1 + exp ( ( _lv + 7.0 ) / 14.0 ) ) ) ;
     }
   
return _lce;
 }
 
static void _hoc_ce(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  Prop* _local_prop = _prop_id ? _extcall_prop : nullptr;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  ce ( _threadargscomma_ *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static double _npy_ce(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r =  ce ( _threadargscomma_ *getarg(1) , *getarg(2) );
 return(_r);
}
 
static int  rates ( _internalthreadargsprotocomma_ double _lv ) {
   double _la , _lb , _lc ;
 _la = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 0.0 ) ;
   _lc = ce ( _threadargscomma_ _lv , 0.0 ) ;
   mtau = 1.0 / ( _la + _lb ) / 3.0 * Tauact ;
   minf = _lc ;
   _la = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 1.0 ) ;
   _lc = ce ( _threadargscomma_ _lv , 1.0 ) ;
   ntau = 1.0 / ( _la + _lb ) / 3.0 * Tauinactf ;
   ninf = _lc ;
   _lc = ce ( _threadargscomma_ _lv , 2.0 ) ;
   uinf = _lc ;
   utau = 6800.0 * Tauinacts ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for rates_Kv1_5. Requires prior call to setdata_Kv1_5 and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r = 1.;
 rates ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_rates(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
_nt = nrn_threads;
 _r = 1.;
 rates ( _threadargscomma_ *getarg(1) );
 return(_r);
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
  ek = _ion_ek;
  ki = _ion_ki;
  ko = _ion_ko;
  nai = _ion_nai;
  nao = _ion_nao;
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
  ek = _ion_ek;
  ki = _ion_ki;
  ko = _ion_ko;
  nai = _ion_nai;
  nao = _ion_nao;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[_dith1] = new double[6]{};
   _newtonspace1 = nrn_cons_newtonspace(3);
 }
 
static void _thread_cleanup(Datum* _thread) {
   delete[] _thread[_dith1].get<double*>();
   nrn_destroy_newtonspace(_newtonspace1);
 }

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  m = m0;
  n = n0;
  u = u0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   n = ninf ;
   u = uinf ;
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
   _v = _vec_v[_ni[_iml]];
 v = _v;
  ek = _ion_ek;
  ki = _ion_ki;
  ko = _ion_ko;
  nai = _ion_nai;
  nao = _ion_nao;
 initmodel(_threadargs_);
  }
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   double _lz ;
 _lz = ( R * ( celsius + 273.15 ) ) / F ;
   ik = gKur * ( 0.1 + 1.0 / ( 1.0 + exp ( - ( v - 15.0 ) / 13.0 ) ) ) * m * m * m * n * u * ( v - ek ) ;
   ino = gnonspec * ( 0.1 + 1.0 / ( 1.0 + exp ( - ( v - 15.0 ) / 13.0 ) ) ) * m * m * m * n * u * ( v - _lz * log ( ( nao + ko ) / ( nai + ki ) ) ) ;
   }
 _current += ik;
 _current += ino;

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
  ki = _ion_ki;
  ko = _ion_ko;
  nai = _ion_nai;
  nao = _ion_nao;
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ double _dino;
 double _dik;
  _dik = ik;
  _dino = ino;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dikdv += (_dik - ik)/.001 ;
  _ion_dinodv += (_dino - ino)/.001 ;
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ik += ik ;
  _ion_ino += ino ;
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
  ki = _ion_ki;
  ko = _ion_ko;
  nai = _ion_nai;
  nao = _ion_nao;
 {  _deriv1_advance = 1;
 derivimplicit_thread(3, _slist1, _dlist1, neuron::scopmath::row_view{_ml, _iml}, states, _ml, _iml, _ppvar, _thread, _nt);
_deriv1_advance = 0;
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _ml->data(_iml, _slist1[_i]) += dt*_ml->data(_iml, _dlist1[_i]);
    }}
 }  }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {m_columnindex, 0};  _dlist1[0] = {Dm_columnindex, 0};
 _slist1[1] = {n_columnindex, 0};  _dlist1[1] = {Dn_columnindex, 0};
 _slist1[2] = {u_columnindex, 0};  _dlist1[2] = {Du_columnindex, 0};
 _slist2[0] = {m_columnindex, 0};
 _slist2[1] = {n_columnindex, 0};
 _slist2[2] = {u_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/Kv15.mod";
    const char* nmodl_file_text = 
  "TITLE Cardiac IKur  current & nonspec cation current with identical kinetics\n"
  ": Hodgkin - Huxley type channels, modified to fit IKur data from Feng et al Am J Physiol 1998 275:H1717 - H 1725\n"
  ": Suffix from Kv15 to Kv1_5\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX Kv1_5\n"
  "	USEION k READ ek,ki,ko WRITE ik\n"
  "	USEION na READ nai,nao\n"
  "	USEION no WRITE ino VALENCE 1: nonspecific cation current\n"
  "	RANGE gKur, ik, ino, Tauact, Tauinactf,Tauinacts, gnonspec, nao, nai, ko,ki\n"
  "	RANGE minf, ninf, uinf, mtau , ntau, utau\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "        (mM) = (milli/liter)\n"
  "	F = (faraday) (coulombs)\n"
  "	R 	= (k-mole)	(joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	 gKur=0.13195e-3 (S/cm2) <0,1e9>\n"
  "	Tauact=1 (ms)\n"
  "	Tauinactf=1 (ms)\n"
  "	Tauinacts=1 (ms)\n"
  "	gnonspec=0   (S/cm2) <0,1e9>\n"
  "}\n"
  "STATE {\n"
  "	 m n u\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	celsius (degC) : 37\n"
  "       	ik (mA/cm2)\n"
  "	minf ninf uinf\n"
  "	mtau (ms)\n"
  "        ntau (ms)\n"
  "	utau (ms)\n"
  "	ek (mV)\n"
  "	ino (mA/cm2)\n"
  "	ki (mM)\n"
  "	ko (mM)\n"
  "	nai (mM)\n"
  "	nao (mM)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = minf\n"
  "        n = ninf\n"
  "	u = uinf\n"
  "}\n"
  "\n"
  "BREAKPOINT { LOCAL z\n"
  "	z = (R*(celsius+273.15))/F\n"
  "	SOLVE states METHOD derivimplicit\n"
  "		ik = gKur*(0.1 + 1/(1 + exp(-(v - 15)/13)))*m*m*m*n*u*(v - ek)\n"
  "	ino=gnonspec*(0.1 + 1/(1 + exp(-(v - 15)/13)))*m*m*m*n*u*(v - z*log((nao+ko)/(nai+ki)))\n"
  "}\n"
  "\n"
  "DERIVATIVE states {	: exact when v held constant\n"
  "	rates(v)\n"
  "	m' = (minf - m)/mtau\n"
  "        n' = (ninf - n)/ntau\n"
  "	u' = (uinf - u)/utau\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "FUNCTION alp(v(mV),i) { LOCAL q10 : order m n\n"
  "	v = v\n"
  "	q10 = 2.2^((celsius - 37)/10)\n"
  "       if (i==0) {\n"
  "	          alp = q10*0.65/(exp(-(v + 10)/8.5) + exp(-(v - 30)/59))\n"
  "          } else if (i==1) {\n"
  "                   alp = 0.001*q10/(2.4 +10.9* exp(-(v + 90)/78))\n"
  "          }\n"
  "	\n"
  "}\n"
  "\n"
  "FUNCTION bet(v(mV),i) (/ms) { LOCAL q10 : order m n u\n"
  "	v = v \n"
  "	q10 = 2.2^((celsius - 37)/10)\n"
  "        if (i==0){\n"
  "	         bet = q10*0.65/(2.5 + exp((v + 82)/17))\n"
  "        }else if (i==1){\n"
  "                  bet = q10*0.001*exp((v - 168)/16)\n"
  "        }\n"
  "}\n"
  "                \n"
  "FUNCTION ce(v(mV),i)(/ms) {   :  order m n u \n"
  "        v = v\n"
  "       \n"
  "        if (i==0) {\n"
  "                ce = 1/(1 + exp(-(v + 30.3)/9.6))\n"
  "        }else if (i==1){\n"
  "                ce = 1*(0.25+1/(1.35 + exp((v + 7)/14)))\n"
  "       \n"
  "	}else if (i==2){\n"
  "                ce = 1*(0.1+1/(1.1 + exp((v + 7)/14)))\n"
  "        }\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(v) {LOCAL a,b,c :\n"
  "	\n"
  "		a = alp(v,0)  b=bet(v,0) c = ce(v,0)\n"
  "		mtau = 1/(a + b)/3*Tauact\n"
  "		minf = c\n"
  "               a = alp(v,1)  b=bet(v,1) c = ce(v,1)\n"
  "		ntau = 1/(a + b)/3*Tauinactf\n"
  "		ninf = c\n"
  "		c = ce(v,2)\n"
  "		uinf = c\n"
  "		utau = 6800*Tauinacts\n"
  "}\n"
  "UNITSON\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
