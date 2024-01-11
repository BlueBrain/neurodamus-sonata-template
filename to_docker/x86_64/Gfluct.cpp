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
static constexpr auto number_of_floating_point_variables = 21;
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
 
#define nrn_init _nrn_init__Gfluct
#define _nrn_initial _nrn_initial__Gfluct
#define nrn_cur _nrn_cur__Gfluct
#define _nrn_current _nrn_current__Gfluct
#define nrn_jacob _nrn_jacob__Gfluct
#define nrn_state _nrn_state__Gfluct
#define _net_receive _net_receive__Gfluct 
#define noiseFromRandom noiseFromRandom__Gfluct 
#define oup oup__Gfluct 
 
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
#define E_e _ml->template fpfield<0>(_iml)
#define E_e_columnindex 0
#define E_i _ml->template fpfield<1>(_iml)
#define E_i_columnindex 1
#define g_e0 _ml->template fpfield<2>(_iml)
#define g_e0_columnindex 2
#define g_i0 _ml->template fpfield<3>(_iml)
#define g_i0_columnindex 3
#define std_e _ml->template fpfield<4>(_iml)
#define std_e_columnindex 4
#define std_i _ml->template fpfield<5>(_iml)
#define std_i_columnindex 5
#define tau_e _ml->template fpfield<6>(_iml)
#define tau_e_columnindex 6
#define tau_i _ml->template fpfield<7>(_iml)
#define tau_i_columnindex 7
#define i _ml->template fpfield<8>(_iml)
#define i_columnindex 8
#define g_e _ml->template fpfield<9>(_iml)
#define g_e_columnindex 9
#define g_i _ml->template fpfield<10>(_iml)
#define g_i_columnindex 10
#define g_e1 _ml->template fpfield<11>(_iml)
#define g_e1_columnindex 11
#define g_i1 _ml->template fpfield<12>(_iml)
#define g_i1_columnindex 12
#define D_e _ml->template fpfield<13>(_iml)
#define D_e_columnindex 13
#define D_i _ml->template fpfield<14>(_iml)
#define D_i_columnindex 14
#define exp_e _ml->template fpfield<15>(_iml)
#define exp_e_columnindex 15
#define exp_i _ml->template fpfield<16>(_iml)
#define exp_i_columnindex 16
#define amp_e _ml->template fpfield<17>(_iml)
#define amp_e_columnindex 17
#define amp_i _ml->template fpfield<18>(_iml)
#define amp_i_columnindex 18
#define v _ml->template fpfield<19>(_iml)
#define v_columnindex 19
#define _g _ml->template fpfield<20>(_iml)
#define _g_columnindex 20
#define _nd_area *_ml->dptr_field<0>(_iml)
#define donotuse	*_ppvar[2].get<double*>()
#define _p_donotuse _ppvar[2].literal_value<void*>()
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  2;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
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
 {"noiseFromRandom", _hoc_noiseFromRandom},
 {"normrand123", _hoc_normrand123},
 {"oup", _hoc_oup},
 {0, 0}
};
#define normrand123 normrand123_Gfluct
 extern double normrand123( _internalthreadargsproto_ );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"E_e", "mV"},
 {"E_i", "mV"},
 {"g_e0", "umho"},
 {"g_i0", "umho"},
 {"std_e", "umho"},
 {"std_i", "umho"},
 {"tau_e", "ms"},
 {"tau_i", "ms"},
 {"i", "nA"},
 {"g_e", "umho"},
 {"g_i", "umho"},
 {"g_e1", "umho"},
 {"g_i1", "umho"},
 {"D_e", "umho umho /ms"},
 {"D_i", "umho umho /ms"},
 {0, 0}
};
 static double delta_t = 0.01;
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
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Gfluct",
 "E_e",
 "E_i",
 "g_e0",
 "g_i0",
 "std_e",
 "std_i",
 "tau_e",
 "tau_i",
 0,
 "i",
 "g_e",
 "g_i",
 "g_e1",
 "g_i1",
 "D_e",
 "D_i",
 0,
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
   _ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 21);
 	/*initialize range parameters*/
 	E_e = 0;
 	E_i = -75;
 	g_e0 = 0.0121;
 	g_i0 = 0.0573;
 	std_e = 0.003;
 	std_i = 0.0066;
 	tau_e = 3;
 	tau_i = 10.49;
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 21);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 static void bbcore_write(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_write(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 static void bbcore_read(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _Gfluct_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"E_e"} /* 0 */,
                                       _nrn_mechanism_field<double>{"E_i"} /* 1 */,
                                       _nrn_mechanism_field<double>{"g_e0"} /* 2 */,
                                       _nrn_mechanism_field<double>{"g_i0"} /* 3 */,
                                       _nrn_mechanism_field<double>{"std_e"} /* 4 */,
                                       _nrn_mechanism_field<double>{"std_i"} /* 5 */,
                                       _nrn_mechanism_field<double>{"tau_e"} /* 6 */,
                                       _nrn_mechanism_field<double>{"tau_i"} /* 7 */,
                                       _nrn_mechanism_field<double>{"i"} /* 8 */,
                                       _nrn_mechanism_field<double>{"g_e"} /* 9 */,
                                       _nrn_mechanism_field<double>{"g_i"} /* 10 */,
                                       _nrn_mechanism_field<double>{"g_e1"} /* 11 */,
                                       _nrn_mechanism_field<double>{"g_i1"} /* 12 */,
                                       _nrn_mechanism_field<double>{"D_e"} /* 13 */,
                                       _nrn_mechanism_field<double>{"D_i"} /* 14 */,
                                       _nrn_mechanism_field<double>{"exp_e"} /* 15 */,
                                       _nrn_mechanism_field<double>{"exp_i"} /* 16 */,
                                       _nrn_mechanism_field<double>{"amp_e"} /* 17 */,
                                       _nrn_mechanism_field<double>{"amp_i"} /* 18 */,
                                       _nrn_mechanism_field<double>{"v"} /* 19 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 20 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"donotuse", "bbcorepointer"} /* 2 */);
  hoc_register_prop_size(_mechtype, 21, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Gfluct /mnt/mydata/cerebellum/mod/Gfluct.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "Fluctuating conductances";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int noiseFromRandom(_internalthreadargsproto_);
static int oup(_internalthreadargsproto_);
 
static int  oup ( _internalthreadargsproto_ ) {
   if ( tau_e  != 0.0 ) {
     g_e1 = exp_e * g_e1 + amp_e * normrand123 ( _threadargs_ ) ;
     }
   if ( tau_i  != 0.0 ) {
     g_i1 = exp_i * g_i1 + amp_i * normrand123 ( _threadargs_ ) ;
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
		_lnormrand123 = nrn_random_pick((Rand*)_p_donotuse);
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
 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* xx, int *offset, _threadargsproto_) {
#if !NRNBBCORE
	/* error if using the legacy normrand */
	if (!_p_donotuse) {
		fprintf(stderr, "Gfluct: cannot use the legacy normrand generator for the random stream.\n");
		assert(0);
	}
	if (d) {
		uint32_t* di = ((uint32_t*)d) + *offset;
		Rand** pv = (Rand**)(&_p_donotuse);
		/* error if not using Random123 generator */
		if (!nrn_random_isran123(*pv, di, di+1, di+2)) {
			fprintf(stderr, "Gfluct: Random123 generator is required\n");
			assert(0);
		}
		/*printf("Gfluct bbcore_write %d %d %d\n", di[0], di[1], di[3]);*/
	}
	*offset += 3;
#endif
}

static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
	assert(!_p_donotuse);
	uint32_t* di = ((uint32_t*)d) + *offset;
	nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);
	*pv = nrnran123_newstream3(di[0], di[1], di[2]);
	*offset += 3;
}
 
static int _ode_count(int _type){ hoc_execerror("Gfluct", "cannot be used with CVODE"); return 0;}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
 {
   g_e1 = 0.0 ;
   g_i1 = 0.0 ;
   if ( tau_e  != 0.0 ) {
     D_e = 2.0 * std_e * std_e / tau_e ;
     exp_e = exp ( - dt / tau_e ) ;
     amp_e = std_e * sqrt ( ( 1.0 - exp ( - 2.0 * dt / tau_e ) ) ) ;
     }
   if ( tau_i  != 0.0 ) {
     D_i = 2.0 * std_i * std_i / tau_i ;
     exp_i = exp ( - dt / tau_i ) ;
     amp_i = std_i * sqrt ( ( 1.0 - exp ( - 2.0 * dt / tau_i ) ) ) ;
     }
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
 initmodel(_threadargs_);
}
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   if ( tau_e  == 0.0 ) {
     g_e = std_e * normrand123 ( _threadargs_ ) ;
     }
   if ( tau_i  == 0.0 ) {
     g_i = std_i * normrand123 ( _threadargs_ ) ;
     }
   g_e = g_e0 + g_e1 ;
   if ( g_e < 0.0 ) {
     g_e = 0.0 ;
     }
   g_i = g_i0 + g_i1 ;
   if ( g_i < 0.0 ) {
     g_i = 0.0 ;
     }
   i = g_e * ( v - E_e ) + g_i * ( v - E_i ) ;
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
 {  { oup(_threadargs_); }
  }}}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/Gfluct.mod";
    const char* nmodl_file_text = 
  "TITLE Fluctuating conductances\n"
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
  "     E_i = -75 (mV)		: reversal potential of inhibitory conductance\n"
  "\n"
  "     g_e0 = 0.0121 (umho)	: average excitatory conductance\n"
  "     g_i0 = 0.0573 (umho)	: average inhibitory conductance\n"
  "\n"
  "     std_e = 0.0030 (umho)	: standard dev of excitatory conductance\n"
  "     std_i = 0.0066 (umho)	: standard dev of inhibitory conductance\n"
  "\n"
  "     tau_e = 2.728 (ms)		: time constant of excitatory conductance\n"
  "     tau_i = 10.49 (ms)		: time constant of inhibitory conductance\n"
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
  "	POINT_PROCESS Gfluct\n"
  "	RANGE g_e, g_i, E_e, E_i, g_e0, g_i0, g_e1, g_i1\n"
  "	RANGE std_e, std_i, tau_e, tau_i, D_e, D_i\n"
  "	NONSPECIFIC_CURRENT i\n"
  "        \n"
  "        THREADSAFE : only true if every instance has its own distinct Random\n"
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
  "	dt		(ms)\n"
  "\n"
  "	E_e	= 0 	(mV)	: reversal potential of excitatory conductance\n"
  "	E_i	= -75 	(mV)	: reversal potential of inhibitory conductance\n"
  "\n"
  "	g_e0	= 0.0121 (umho)	: average excitatory conductance\n"
  "	g_i0	= 0.0573 (umho)	: average inhibitory conductance\n"
  "\n"
  "	std_e	= 0.0030 (umho)	: standard dev of excitatory conductance\n"
  "	std_i	= 0.0066 (umho)	: standard dev of inhibitory conductance\n"
  "\n"
  "	tau_e	= 3	(ms)	: time constant of excitatory conductance\n"
  "	tau_i	= 10.49	(ms)	: time constant of inhibitory conductance\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV)		: membrane voltage\n"
  "	i 	(nA)		: fluctuating current\n"
  "	g_e	(umho)		: total excitatory conductance\n"
  "	g_i	(umho)		: total inhibitory conductance\n"
  "	g_e1	(umho)		: fluctuating excitatory conductance\n"
  "	g_i1	(umho)		: fluctuating inhibitory conductance\n"
  "	D_e	(umho umho /ms) : excitatory diffusion coefficient\n"
  "	D_i	(umho umho /ms) : inhibitory diffusion coefficient\n"
  "	exp_e\n"
  "	exp_i\n"
  "	amp_e	(umho)\n"
  "	amp_i	(umho)\n"
  "\n"
  "        donotuse\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	g_e1 = 0\n"
  "	g_i1 = 0\n"
  "	if(tau_e != 0) {\n"
  "		D_e = 2 * std_e * std_e / tau_e\n"
  "		exp_e = exp(-dt/tau_e)\n"
  "		amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )\n"
  "	}\n"
  "	if(tau_i != 0) {\n"
  "		D_i = 2 * std_i * std_i / tau_i\n"
  "		exp_i = exp(-dt/tau_i)\n"
  "		amp_i = std_i * sqrt( (1-exp(-2*dt/tau_i)) )\n"
  "	}\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE oup\n"
  "	if(tau_e==0) {\n"
  "	   g_e = std_e * normrand123()\n"
  "	}\n"
  "	if(tau_i==0) {\n"
  "	   g_i = std_i * normrand123()\n"
  "	}\n"
  "	g_e = g_e0 + g_e1\n"
  "	if(g_e < 0) { g_e = 0 }\n"
  "	g_i = g_i0 + g_i1\n"
  "	if(g_i < 0) { g_i = 0 }\n"
  "	i = g_e * (v - E_e) + g_i * (v - E_i)\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE oup() {		: use Scop function normrand(mean, std_dev)\n"
  "   if(tau_e!=0) {\n"
  "	g_e1 =  exp_e * g_e1 + amp_e * normrand123()\n"
  "   }\n"
  "   if(tau_i!=0) {\n"
  "	g_i1 =  exp_i * g_i1 + amp_i * normrand123()\n"
  "   }\n"
  "}\n"
  "\n"
  "FUNCTION normrand123() {\n"
  "VERBATIM\n"
  "	if (_p_donotuse) {\n"
  "		/*\n"
  "		:Supports separate independent but reproducible streams for\n"
  "		: each instance. However, the corresponding hoc Random\n"
  "		: distribution MUST be set to Random.negexp(1)\n"
  "		*/\n"
  "        #if !NRNBBCORE\n"
  "		_lnormrand123 = nrn_random_pick((Rand*)_p_donotuse);\n"
  "        #else\n"
  "        #pragma acc routine(nrnran123_normal) seq\n"
  "        _lnormrand123 = nrnran123_normal((nrnran123_State*)_p_donotuse);\n"
  "        #endif\n"
  "	}else{\n"
  "		/* only use Random123 */\n"
  "        assert(0);\n"
  "	}\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE noiseFromRandom() {\n"
  "VERBATIM\n"
  "#if !NRNBBCORE\n"
  " {\n"
  "	void** pv = (void**)(&_p_donotuse);\n"
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
  "VERBATIM\n"
  "static void bbcore_write(double* x, int* d, int* xx, int *offset, _threadargsproto_) {\n"
  "#if !NRNBBCORE\n"
  "	/* error if using the legacy normrand */\n"
  "	if (!_p_donotuse) {\n"
  "		fprintf(stderr, \"Gfluct: cannot use the legacy normrand generator for the random stream.\\n\");\n"
  "		assert(0);\n"
  "	}\n"
  "	if (d) {\n"
  "		uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "		Rand** pv = (Rand**)(&_p_donotuse);\n"
  "		/* error if not using Random123 generator */\n"
  "		if (!nrn_random_isran123(*pv, di, di+1, di+2)) {\n"
  "			fprintf(stderr, \"Gfluct: Random123 generator is required\\n\");\n"
  "			assert(0);\n"
  "		}\n"
  "		/*printf(\"Gfluct bbcore_write %d %d %d\\n\", di[0], di[1], di[3]);*/\n"
  "	}\n"
  "	*offset += 3;\n"
  "#endif\n"
  "}\n"
  "\n"
  "static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {\n"
  "	assert(!_p_donotuse);\n"
  "	uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "	nrnran123_State** pv = (nrnran123_State**)(&_p_donotuse);\n"
  "	*pv = nrnran123_newstream3(di[0], di[1], di[2]);\n"
  "	*offset += 3;\n"
  "}\n"
  "ENDVERBATIM\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
