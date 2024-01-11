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
 
#define nrn_init _nrn_init__orn
#define _nrn_initial _nrn_initial__orn
#define nrn_cur _nrn_cur__orn
#define _nrn_current _nrn_current__orn
#define nrn_jacob _nrn_jacob__orn
#define nrn_state _nrn_state__orn
#define _net_receive _net_receive__orn 
#define initstream initstream__orn 
#define noiseFromRandom123 noiseFromRandom123__orn 
#define noiseFromRandom noiseFromRandom__orn 
#define oup oup__orn 
#define states states__orn 
 
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
#define g_e_max _ml->template fpfield<0>(_iml)
#define g_e_max_columnindex 0
#define cc_peak _ml->template fpfield<1>(_iml)
#define cc_peak_columnindex 1
#define g_e_baseline _ml->template fpfield<2>(_iml)
#define g_e_baseline_columnindex 2
#define std_e _ml->template fpfield<3>(_iml)
#define std_e_columnindex 3
#define tau_e _ml->template fpfield<4>(_iml)
#define tau_e_columnindex 4
#define i _ml->template fpfield<5>(_iml)
#define i_columnindex 5
#define g_e _ml->template fpfield<6>(_iml)
#define g_e_columnindex 6
#define D_e _ml->template fpfield<7>(_iml)
#define D_e_columnindex 7
#define O _ml->template fpfield<8>(_iml)
#define O_columnindex 8
#define C _ml->template fpfield<9>(_iml)
#define C_columnindex 9
#define D _ml->template fpfield<10>(_iml)
#define D_columnindex 10
#define g_e1 _ml->template fpfield<11>(_iml)
#define g_e1_columnindex 11
#define exp_e _ml->template fpfield<12>(_iml)
#define exp_e_columnindex 12
#define amp_e _ml->template fpfield<13>(_iml)
#define amp_e_columnindex 13
#define DO _ml->template fpfield<14>(_iml)
#define DO_columnindex 14
#define DC _ml->template fpfield<15>(_iml)
#define DC_columnindex 15
#define DD _ml->template fpfield<16>(_iml)
#define DD_columnindex 16
#define v _ml->template fpfield<17>(_iml)
#define v_columnindex 17
#define _g _ml->template fpfield<18>(_iml)
#define _g_columnindex 18
#define _tsav _ml->template fpfield<19>(_iml)
#define _tsav_columnindex 19
#define _nd_area *_ml->dptr_field<0>(_iml)
#define donotuse	*_ppvar[2].get<double*>()
#define _p_donotuse _ppvar[2].literal_value<void*>()
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  2;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_initstream(void*);
 static double _hoc_noiseFromRandom123(void*);
 static double _hoc_noiseFromRandom(void*);
 static double _hoc_normrand123(void*);
 static double _hoc_oup(void*);
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
 {"initstream", _hoc_initstream},
 {"noiseFromRandom123", _hoc_noiseFromRandom123},
 {"noiseFromRandom", _hoc_noiseFromRandom},
 {"normrand123", _hoc_normrand123},
 {"oup", _hoc_oup},
 {0, 0}
};
#define normrand123 normrand123_orn
 extern double normrand123( _internalthreadargsproto_ );
 /* declare global and static user variables */
#define E_e E_e_orn
 double E_e = 0;
#define net_receive_on net_receive_on_orn
 double net_receive_on = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"E_e_orn", "mV"},
 {"g_e_max", "umho"},
 {"g_e_baseline", "umho"},
 {"std_e", "umho"},
 {"tau_e", "ms"},
 {"i", "nA"},
 {"g_e", "umho"},
 {"D_e", "umho umho /ms"},
 {0, 0}
};
 static double C0 = 0;
 static double D0 = 0;
 static double O0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"E_e_orn", &E_e_orn},
 {"net_receive_on_orn", &net_receive_on_orn},
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
"orn",
 "g_e_max",
 "cc_peak",
 "g_e_baseline",
 "std_e",
 "tau_e",
 0,
 "i",
 "g_e",
 "D_e",
 0,
 "O",
 "C",
 "D",
 0,
 "donotuse",
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
    assert(_nrn_mechanism_get_num_vars(_prop) == 20);
 	/*initialize range parameters*/
 	g_e_max = 0.075;
 	cc_peak = 0;
 	g_e_baseline = 0;
 	std_e = 0.001;
 	tau_e = 400;
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 20);
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
 static void _net_init(Point_process*, double*, double);
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void bbcore_write(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_write(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 static void bbcore_read(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _orn_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 5,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
  _extcall_thread.resize(4);
  _thread_mem_init(_extcall_thread.data());
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"g_e_max"} /* 0 */,
                                       _nrn_mechanism_field<double>{"cc_peak"} /* 1 */,
                                       _nrn_mechanism_field<double>{"g_e_baseline"} /* 2 */,
                                       _nrn_mechanism_field<double>{"std_e"} /* 3 */,
                                       _nrn_mechanism_field<double>{"tau_e"} /* 4 */,
                                       _nrn_mechanism_field<double>{"i"} /* 5 */,
                                       _nrn_mechanism_field<double>{"g_e"} /* 6 */,
                                       _nrn_mechanism_field<double>{"D_e"} /* 7 */,
                                       _nrn_mechanism_field<double>{"O"} /* 8 */,
                                       _nrn_mechanism_field<double>{"C"} /* 9 */,
                                       _nrn_mechanism_field<double>{"D"} /* 10 */,
                                       _nrn_mechanism_field<double>{"g_e1"} /* 11 */,
                                       _nrn_mechanism_field<double>{"exp_e"} /* 12 */,
                                       _nrn_mechanism_field<double>{"amp_e"} /* 13 */,
                                       _nrn_mechanism_field<double>{"DO"} /* 14 */,
                                       _nrn_mechanism_field<double>{"DC"} /* 15 */,
                                       _nrn_mechanism_field<double>{"DD"} /* 16 */,
                                       _nrn_mechanism_field<double>{"v"} /* 17 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 18 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 19 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"donotuse", "bbcorepointer"} /* 2 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 3 */);
  hoc_register_prop_size(_mechtype, 20, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 4;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 orn /mnt/mydata/cerebellum/mod/orn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "Fluctuating    conductances";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int initstream(_internalthreadargsproto_);
static int noiseFromRandom123(_internalthreadargsproto_);
static int noiseFromRandom(_internalthreadargsproto_);
static int oup(_internalthreadargsproto_);
 
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
   double _lKO , _lKC1 , _lKC2 , _lKD1 , _lKD2 ;
 _lKO = 1.0 / 100.0 ;
   _lKC1 = 1.0 / 100.0 ;
   _lKC2 = 1e-4 ;
   _lKD1 = 1.0 / 6000.0 ;
   _lKD2 = 1.0 / 100.0 ;
   DO = _lKO * ( 1.0 - C - O ) ;
   DC = _lKC1 * ( 1.0 - C ) * C + _lKC2 * ( 1.0 - C ) ;
   DD = _lKD1 * O * ( 1.0 - D ) - _lKD2 * D * ( 1.0 - O ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 double _lKO , _lKC1 , _lKC2 , _lKD1 , _lKD2 ;
 _lKO = 1.0 / 100.0 ;
 _lKC1 = 1.0 / 100.0 ;
 _lKC2 = 1e-4 ;
 _lKD1 = 1.0 / 6000.0 ;
 _lKD2 = 1.0 / 100.0 ;
 DO = DO  / (1. - dt*( ( _lKO )*( ( ( - 1.0 ) ) ) )) ;
 DC = DC  / (1. - dt*( (( ( _lKC1 )*( ( ( - 1.0 ) ) ) )*( C ) + ( _lKC1 * ( 1.0 - C ) )*( 1.0 )) + ( _lKC2 )*( ( ( - 1.0 ) ) ) )) ;
 DD = DD  / (1. - dt*( ( _lKD1 * O )*( ( ( - 1.0 ) ) ) - ( ( _lKD2 )*( 1.0 ) )*( ( 1.0 - O ) ) )) ;
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
   double _lKO , _lKC1 , _lKC2 , _lKD1 , _lKD2 ;
 _lKO = 1.0 / 100.0 ;
   _lKC1 = 1.0 / 100.0 ;
   _lKC2 = 1e-4 ;
   _lKD1 = 1.0 / 6000.0 ;
   _lKD2 = 1.0 / 100.0 ;
   DO = _lKO * ( 1.0 - C - O ) ;
   DC = _lKC1 * ( 1.0 - C ) * C + _lKC2 * ( 1.0 - C ) ;
   DD = _lKD1 * O * ( 1.0 - D ) - _lKD2 * D * ( 1.0 - O ) ;
   {int _id; for(_id=0; _id < 3; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _ml->data(_iml, _dlist1[_id]) - (_ml->data(_iml, _slist1[_id]) - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _ml->data(_iml, _slist1[_id]) - _savstate1[_id];}}}
 } }
 return _reset;}
 
static int  oup ( _internalthreadargsproto_ ) {
   if ( tau_e  != 0.0 ) {
     g_e1 = exp_e * g_e1 + amp_e * normrand123 ( _threadargs_ ) ;
     }
    return 0; }
 
static double _hoc_oup(void* _vptr) {
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
 _r = 1.;
 oup ( _threadargs_ );
 return(_r);
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  Prop* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
   _thread = nullptr; _nt = (NrnThread*)_pnt->_vnt;   _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 3;
    double __state = C;
    double __primary_delta = (0.0) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _ml->data(_iml, _dlist1[__i]) = 0.0;
    }
    _ml->data(_iml, _dlist1[1]) = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _ml->data(_iml, _slist1[__i]) += _ml->data(_iml, _dlist1[__i]);
    }
  } else {
 C = 0.0 ;
     }
 if ( net_receive_on ) {
     cc_peak = _args[0] ;
     g_e_baseline = _args[1] ;
     g_e_max = _args[2] ;
     std_e = _args[3] ;
     }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
     _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  Datum* _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  Datum* _thread = (Datum*)0;
  NrnThread* _nt = (NrnThread*)_pnt->_vnt;
 }
 
/*VERBATIM*/
#if !NRNBBCORE /* running in NEURON */
/*
   1 means noiseFromRandom was called when _ran_compat was previously 0 .
   2 means noiseFromRandom123 was called when _ran_compat was previously 0.
*/
static int _ran_compat; /* specifies the noise style for all instances */
#endif /* running in NEURON */
 
static int  initstream ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
  if (_p_donotuse) {
    uint32_t id1, id2, id3;
    #if NRNBBCORE
      nrnran123_setseq((nrnran123_State*)_p_donotuse, 0, 0);
    #else
      if (_ran_compat == 1) {
        nrn_random_reset((Rand*)_p_donotuse);
      }else{
        nrnran123_setseq((nrnran123_State*)_p_donotuse, 0, 0);
      }
    #endif
  }	
  return 0; }
 
static double _hoc_initstream(void* _vptr) {
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
 _r = 1.;
 initstream ( _threadargs_ );
 return(_r);
}
 
double normrand123 ( _internalthreadargsproto_ ) {
   double _lnormrand123;
 
/*VERBATIM*/
  if (_p_donotuse) {
    /*
      :Supports separate independent but reproducible streams for
      : each instance. However, the corresponding hoc Random
      : distribution MUST be set to Random.negexp(1)
    */
    #if !NRNBBCORE
      if(_ran_compat == 1) {
        _lnormrand123 = nrn_random_pick((Rand*)_p_donotuse);
      }else{
        _lnormrand123 = nrnran123_normal((nrnran123_State*)_p_donotuse);
      }
    #else
      #pragma acc routine(nrnran123_normal) seq
      _lnormrand123 = nrnran123_normal((nrnran123_State*)_p_donotuse);
    #endif
  }else{
    /* only use Random123 */
    assert(0);
  }
 
return _lnormrand123;
 }
 
static double _hoc_normrand123(void* _vptr) {
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
 _r =  normrand123 ( _threadargs_ );
 return(_r);
}
 
static int  noiseFromRandom ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
#if !NRNBBCORE
 {
	void** pv = (void**)(&_p_donotuse);
	if (_ran_compat == 2) {
		fprintf(stderr, "orn.noiseFromRandom123 was previously called\n");
		assert(0);
	}
	_ran_compat = 1;
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
#endif
  return 0; }
 
static double _hoc_noiseFromRandom(void* _vptr) {
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
 _r = 1.;
 noiseFromRandom ( _threadargs_ );
 return(_r);
}
 
static int  noiseFromRandom123 ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
#if !NRNBBCORE
 {
        nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);
        if (_ran_compat == 1) {
          fprintf(stderr, "orn.noiseFromRandom was previously called\n");
          assert(0);
        }
        _ran_compat = 2;
        if (*pv) {
          nrnran123_deletestream(*pv);
          *pv = (nrnran123_State*)0;
        }
        if (ifarg(3)) {
	      *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));
        }else if (ifarg(2)) {
	      *pv = nrnran123_newstream((uint32_t)*getarg(1), (uint32_t)*getarg(2));
        }
 }
#endif
  return 0; }
 
static double _hoc_noiseFromRandom123(void* _vptr) {
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
 _r = 1.;
 noiseFromRandom123 ( _threadargs_ );
 return(_r);
}
 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* xx, int *offset, _threadargsproto_) {
#if !NRNBBCORE
	/* error if using the legacy normrand */
	if (!_p_donotuse) {
		fprintf(stderr, "orn: cannot use the legacy normrand generator for the random stream.\n");
		assert(0);
	}
	if (d) {
		uint32_t* di = ((uint32_t*)d) + *offset;
		if (_ran_compat == 1) {
			Rand** pv = (Rand**)(&_p_donotuse);
			/* error if not using Random123 generator */
			if (!nrn_random_isran123(*pv, di, di+1, di+2)) {
				fprintf(stderr, "orn: Random123 generator is required\n");
				assert(0);
			}
		}else{
			nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);
			nrnran123_getids3(*pv, di, di+1, di+2);
		}
		/*printf("orn bbcore_write %d %d %d\n", di[0], di[1], di[3]);*/
	}
#endif
	*offset += 3;
}

static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
	uint32_t* di = ((uint32_t*)d) + *offset;
	nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);
#if !NRNBBCORE
    if(*pv) {
        nrnran123_deletestream(*pv);
    }
#endif
	*pv = nrnran123_newstream3(di[0], di[1], di[2]);
	*offset += 3;
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
  C = C0;
  D = D0;
  O = O0;
 {
   initstream ( _threadargs_ ) ;
   g_e1 = 0.0 ;
   if ( tau_e  != 0.0 ) {
     D_e = 2.0 * std_e * std_e / tau_e ;
     exp_e = exp ( - dt / tau_e ) ;
     amp_e = std_e * sqrt ( ( 1.0 - exp ( - 2.0 * dt / tau_e ) ) ) ;
     }
   O = 0.0 ;
   C = 1.0 ;
   D = 0.0 ;
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
   double _lSORN ;
 _lSORN = O * ( 1.0 - D ) ;
   g_e = g_e1 + _lSORN * cc_peak * g_e_max + g_e_baseline ;
   if ( g_e < 0.0 ) {
     g_e = 0.0 ;
     }
   i = g_e * ( v - E_e ) ;
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
 {  { oup(_threadargs_); }
  } {  _deriv1_advance = 1;
 derivimplicit_thread(3, _slist1, _dlist1, neuron::scopmath::row_view{_ml, _iml}, states, _ml, _iml, _ppvar, _thread, _nt);
_deriv1_advance = 0;
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _ml->data(_iml, _slist1[_i]) += dt*_ml->data(_iml, _dlist1[_i]);
    }}
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {O_columnindex, 0};  _dlist1[0] = {DO_columnindex, 0};
 _slist1[1] = {C_columnindex, 0};  _dlist1[1] = {DC_columnindex, 0};
 _slist1[2] = {D_columnindex, 0};  _dlist1[2] = {DD_columnindex, 0};
 _slist2[0] = {C_columnindex, 0};
 _slist2[1] = {D_columnindex, 0};
 _slist2[2] = {O_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/orn.mod";
    const char* nmodl_file_text = 
  "TITLE Fluctuating    conductances\n"
  "\n"
  "COMMENT\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "	Fluctuating conductance model for synaptic bombardment\n"
  "	======================================================\n"
  "\n"
  "THEORY\n"
  "\n"
  "  Synaptic bombardment is represented by a stochastic model containing\n"
  "  two fluctuating conductances g_e(t) and g_i(t) descibed by:\n"
  "\n"
  "     Isyn = g_e(t) * [V - E_e] + g_i(t) * [V - E_i]\n"
  "     d g_e / dt = -(g_e - g_e0) / tau_e + sqrt(D_e) * Ft\n"
  "     d g_i / dt = -(g_i - g_i0) / tau_i + sqrt(D_i) * Ft\n"
  "\n"
  "  where E_e, E_i are the reversal potentials, g_e0, g_i0 are the average\n"
  "  conductances, tau_e, tau_i are time constants, D_e, D_i are noise diffusion\n"
  "  coefficients and Ft is a gaussian white noise of unit standard deviation.\n"
  "\n"
  "  g_e and g_i are described by an Ornstein-Uhlenbeck (OU) stochastic process\n"
  "  where tau_e and tau_i represent the \"correlation\" (if tau_e and tau_i are \n"
  "  zero, g_e and g_i are white noise).  The estimation of OU parameters can\n"
  "  be made from the power spectrum:\n"
  "\n"
  "     S(w) =  2 * D * tau^2 / (1 + w^2 * tau^2)\n"
  "\n"
  "  and the diffusion coeffient D is estimated from the variance:\n"
  "\n"
  "     D = 2 * sigma^2 / tau\n"
  "\n"
  "\n"
  "NUMERICAL RESOLUTION\n"
  "\n"
  "  The numerical scheme for integration of OU processes takes advantage \n"
  "  of the fact that these processes are gaussian, which led to an exact\n"
  "  update rule independent of the time step dt (see Gillespie DT, Am J Phys \n"
  "  64: 225, 1996):\n"
  "\n"
  "     x(t+dt) = x(t) * exp(-dt/tau) + A * N(0,1)\n"
  "\n"
  "  where A = sqrt( D*tau/2 * (1-exp(-2*dt/tau)) ) and N(0,1) is a normal\n"
  "  random number (avg=0, sigma=1)\n"
  "\n"
  "\n"
  "IMPLEMENTATION\n"
  "\n"
  "  This mechanism is implemented as a nonspecific current defined as a\n"
  "  point process.\n"
  "\n"
  "\n"
  "PARAMETERS\n"
  "\n"
  "  The mechanism takes the following parameters:\n"
  "\n"
  "     E_e = 0  (mV)		: reversal potential of excitatory conductance\n"
  "\n"
  "     g_e0 = 0.0121 (umho)	: average excitatory conductance\n"
  "\n"
  "     std_e = 0.0030 (umho)	: standard dev of excitatory conductance\n"
  "\n"
  "     tau_e = 2.728 (ms)		: time constant of excitatory conductance\n"
  "\n"
  "\n"
  "Gfluct3: conductance cannot be negative\n"
  "\n"
  "\n"
  "REFERENCE\n"
  "\n"
  "  Destexhe, A., Rudolph, M., Fellous, J-M. and Sejnowski, T.J.  \n"
  "  Fluctuating synaptic conductances recreate in-vivo--like activity in\n"
  "  neocortical neurons. Neuroscience 107: 13-24 (2001).\n"
  "\n"
  "  (electronic copy available at http://cns.iaf.cnrs-gif.fr)\n"
  "\n"
  "\n"
  "  A. Destexhe, 1999\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "        THREADSAFE : only true if every instance has its own distinct Random\n"
  "	POINT_PROCESS orn\n"
  "	RANGE g_e, g_e_max, cc_peak, g_e_baseline\n"
  "	RANGE std_e, tau_e, D_e\n"
  "	GLOBAL net_receive_on\n"
  "	NONSPECIFIC_CURRENT i\n"
  "        BBCOREPOINTER donotuse\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp) \n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	dt		  (ms)\n"
  "        \n"
  "	E_e	= 0 	  (mV)            : reversal potential of excitatory conductance\n"
  "	g_e_max	= 75e-3 (umho)          : average excitatory conductance\n"
  "        cc_peak = 0     : (affinity*odor cc)\n"
  "        \n"
  "        g_e_baseline    = 0       (umho)          : background noise\n"
  "	std_e	= 1e-3    (umho)	  : standard dev of excitatory conductance\n"
  "	tau_e	= 400     (ms)            : time constant of excitatory conductance\n"
  "	net_receive_on = 0\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV)		: membrane voltage\n"
  "	i 	(nA)		: fluctuating current\n"
  "	g_e	(umho)		: total excitatory conductance\n"
  "	g_e1	(umho)		: fluctuating excitatory conductance\n"
  "	D_e	(umho umho /ms) : excitatory diffusion coefficient\n"
  "	exp_e\n"
  "	amp_e	(umho)\n"
  "        \n"
  "        donotuse\n"
  "}\n"
  "\n"
  "STATE { O C D }\n"
  "\n"
  "INITIAL {\n"
  "	initstream()\n"
  "	g_e1 = 0\n"
  "	if(tau_e != 0) {\n"
  "		D_e = 2 * std_e * std_e / tau_e\n"
  "		exp_e = exp(-dt/tau_e)\n"
  "		amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )\n"
  "	}\n"
  "        \n"
  "        O = 0\n"
  "        C = 1\n"
  "        D = 0\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "        LOCAL SORN\n"
  "	SOLVE oup\n"
  "        SOLVE states METHOD derivimplicit\n"
  "        \n"
  "        SORN = O * (1-D)\n"
  "        \n"
  "        g_e = g_e1 + SORN * cc_peak * g_e_max + g_e_baseline\n"
  "          \n"
  "\n"
  "        if(g_e < 0) {\n"
  "            g_e = 0\n"
  "        }\n"
  "\n"
  "        i = g_e * (v - E_e)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "  LOCAL KO, KC1, KC2, KD1, KD2\n"
  "  KO = 1/100\n"
  "  KC1 = 1/100\n"
  "  KC2 = 1e-4\n"
  "  KD1 = 1/6000\n"
  "  KD2 = 1/100\n"
  "  O' = KO*(1-C-O)\n"
  "  C' = KC1*(1-C)*C + KC2*(1-C)\n"
  "  D' = KD1*O*(1-D) - KD2*D*(1-O)\n"
  "}\n"
  "\n"
  "PROCEDURE oup() {		: use Scop function normrand(mean, std_dev)\n"
  "   if(tau_e!=0) {\n"
  "	g_e1 =  exp_e * g_e1 + amp_e * normrand123()\n"
  "   }\n"
  "}\n"
  "\n"
  "NET_RECEIVE(peak, base, gemax, stde) {\n"
  "  INITIAL {} : do not want to initialize base, etc. to 0.0\n"
  "  C = 0\n"
  "  : off for original implementation where OdorStim fills in following\n"
  "  : at the interpreter level and then does a direct NetCon.event\n"
  "  : This is here to allow debugging with prcellstate.\n"
  "  if (net_receive_on) {\n"
  "    cc_peak = peak\n"
  "    g_e_baseline = base\n"
  "    g_e_max = gemax\n"
  "    std_e = stde\n"
  "  }\n"
  "}\n"
  "\n"
  "\n"
  "VERBATIM\n"
  "#if !NRNBBCORE /* running in NEURON */\n"
  "/*\n"
  "   1 means noiseFromRandom was called when _ran_compat was previously 0 .\n"
  "   2 means noiseFromRandom123 was called when _ran_compat was previously 0.\n"
  "*/\n"
  "static int _ran_compat; /* specifies the noise style for all instances */\n"
  "#endif /* running in NEURON */\n"
  "ENDVERBATIM\n"
  "\n"
  "PROCEDURE initstream() {\n"
  "VERBATIM\n"
  "  if (_p_donotuse) {\n"
  "    uint32_t id1, id2, id3;\n"
  "    #if NRNBBCORE\n"
  "      nrnran123_setseq((nrnran123_State*)_p_donotuse, 0, 0);\n"
  "    #else\n"
  "      if (_ran_compat == 1) {\n"
  "        nrn_random_reset((Rand*)_p_donotuse);\n"
  "      }else{\n"
  "        nrnran123_setseq((nrnran123_State*)_p_donotuse, 0, 0);\n"
  "      }\n"
  "    #endif\n"
  "  }	\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION normrand123() {\n"
  "VERBATIM\n"
  "  if (_p_donotuse) {\n"
  "    /*\n"
  "      :Supports separate independent but reproducible streams for\n"
  "      : each instance. However, the corresponding hoc Random\n"
  "      : distribution MUST be set to Random.negexp(1)\n"
  "    */\n"
  "    #if !NRNBBCORE\n"
  "      if(_ran_compat == 1) {\n"
  "        _lnormrand123 = nrn_random_pick((Rand*)_p_donotuse);\n"
  "      }else{\n"
  "        _lnormrand123 = nrnran123_normal((nrnran123_State*)_p_donotuse);\n"
  "      }\n"
  "    #else\n"
  "      #pragma acc routine(nrnran123_normal) seq\n"
  "      _lnormrand123 = nrnran123_normal((nrnran123_State*)_p_donotuse);\n"
  "    #endif\n"
  "  }else{\n"
  "    /* only use Random123 */\n"
  "    assert(0);\n"
  "  }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE noiseFromRandom() {\n"
  "VERBATIM\n"
  "#if !NRNBBCORE\n"
  " {\n"
  "	void** pv = (void**)(&_p_donotuse);\n"
  "	if (_ran_compat == 2) {\n"
  "		fprintf(stderr, \"orn.noiseFromRandom123 was previously called\\n\");\n"
  "		assert(0);\n"
  "	}\n"
  "	_ran_compat = 1;\n"
  "	if (ifarg(1)) {\n"
  "		*pv = nrn_random_arg(1);\n"
  "	}else{\n"
  "		*pv = (void*)0;\n"
  "	}\n"
  " }\n"
  "#endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE noiseFromRandom123() {\n"
  "VERBATIM\n"
  "#if !NRNBBCORE\n"
  " {\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);\n"
  "        if (_ran_compat == 1) {\n"
  "          fprintf(stderr, \"orn.noiseFromRandom was previously called\\n\");\n"
  "          assert(0);\n"
  "        }\n"
  "        _ran_compat = 2;\n"
  "        if (*pv) {\n"
  "          nrnran123_deletestream(*pv);\n"
  "          *pv = (nrnran123_State*)0;\n"
  "        }\n"
  "        if (ifarg(3)) {\n"
  "	      *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));\n"
  "        }else if (ifarg(2)) {\n"
  "	      *pv = nrnran123_newstream((uint32_t)*getarg(1), (uint32_t)*getarg(2));\n"
  "        }\n"
  " }\n"
  "#endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "static void bbcore_write(double* x, int* d, int* xx, int *offset, _threadargsproto_) {\n"
  "#if !NRNBBCORE\n"
  "	/* error if using the legacy normrand */\n"
  "	if (!_p_donotuse) {\n"
  "		fprintf(stderr, \"orn: cannot use the legacy normrand generator for the random stream.\\n\");\n"
  "		assert(0);\n"
  "	}\n"
  "	if (d) {\n"
  "		uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "		if (_ran_compat == 1) {\n"
  "			Rand** pv = (Rand**)(&_p_donotuse);\n"
  "			/* error if not using Random123 generator */\n"
  "			if (!nrn_random_isran123(*pv, di, di+1, di+2)) {\n"
  "				fprintf(stderr, \"orn: Random123 generator is required\\n\");\n"
  "				assert(0);\n"
  "			}\n"
  "		}else{\n"
  "			nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);\n"
  "			nrnran123_getids3(*pv, di, di+1, di+2);\n"
  "		}\n"
  "		/*printf(\"orn bbcore_write %d %d %d\\n\", di[0], di[1], di[3]);*/\n"
  "	}\n"
  "#endif\n"
  "	*offset += 3;\n"
  "}\n"
  "\n"
  "static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {\n"
  "	uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "	nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);\n"
  "#if !NRNBBCORE\n"
  "    if(*pv) {\n"
  "        nrnran123_deletestream(*pv);\n"
  "    }\n"
  "#endif\n"
  "	*pv = nrnran123_newstream3(di[0], di[1], di[2]);\n"
  "	*offset += 3;\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
