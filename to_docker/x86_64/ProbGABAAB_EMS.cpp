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
static constexpr auto number_of_datum_variables = 5;
static constexpr auto number_of_floating_point_variables = 43;
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
 
#define nrn_init _nrn_init__ProbGABAAB_EMS
#define _nrn_initial _nrn_initial__ProbGABAAB_EMS
#define nrn_cur _nrn_cur__ProbGABAAB_EMS
#define _nrn_current _nrn_current__ProbGABAAB_EMS
#define nrn_jacob _nrn_jacob__ProbGABAAB_EMS
#define nrn_state _nrn_state__ProbGABAAB_EMS
#define _net_receive _net_receive__ProbGABAAB_EMS 
#define clearRNG clearRNG__ProbGABAAB_EMS 
#define setRNG setRNG__ProbGABAAB_EMS 
#define state state__ProbGABAAB_EMS 
#define setup_delay_vecs setup_delay_vecs__ProbGABAAB_EMS 
 
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
#define tau_r_GABAA _ml->template fpfield<0>(_iml)
#define tau_r_GABAA_columnindex 0
#define tau_d_GABAA _ml->template fpfield<1>(_iml)
#define tau_d_GABAA_columnindex 1
#define Use _ml->template fpfield<2>(_iml)
#define Use_columnindex 2
#define Dep _ml->template fpfield<3>(_iml)
#define Dep_columnindex 3
#define Fac _ml->template fpfield<4>(_iml)
#define Fac_columnindex 4
#define e_GABAA _ml->template fpfield<5>(_iml)
#define e_GABAA_columnindex 5
#define e_GABAB _ml->template fpfield<6>(_iml)
#define e_GABAB_columnindex 6
#define u0 _ml->template fpfield<7>(_iml)
#define u0_columnindex 7
#define Nrrp _ml->template fpfield<8>(_iml)
#define Nrrp_columnindex 8
#define synapseID _ml->template fpfield<9>(_iml)
#define synapseID_columnindex 9
#define verboseLevel _ml->template fpfield<10>(_iml)
#define verboseLevel_columnindex 10
#define selected_for_report _ml->template fpfield<11>(_iml)
#define selected_for_report_columnindex 11
#define GABAB_ratio _ml->template fpfield<12>(_iml)
#define GABAB_ratio_columnindex 12
#define conductance _ml->template fpfield<13>(_iml)
#define conductance_columnindex 13
#define i _ml->template fpfield<14>(_iml)
#define i_columnindex 14
#define i_GABAA _ml->template fpfield<15>(_iml)
#define i_GABAA_columnindex 15
#define i_GABAB _ml->template fpfield<16>(_iml)
#define i_GABAB_columnindex 16
#define g_GABAA _ml->template fpfield<17>(_iml)
#define g_GABAA_columnindex 17
#define g_GABAB _ml->template fpfield<18>(_iml)
#define g_GABAB_columnindex 18
#define g _ml->template fpfield<19>(_iml)
#define g_columnindex 19
#define A_GABAA_step _ml->template fpfield<20>(_iml)
#define A_GABAA_step_columnindex 20
#define B_GABAA_step _ml->template fpfield<21>(_iml)
#define B_GABAA_step_columnindex 21
#define A_GABAB_step _ml->template fpfield<22>(_iml)
#define A_GABAB_step_columnindex 22
#define B_GABAB_step _ml->template fpfield<23>(_iml)
#define B_GABAB_step_columnindex 23
#define unoccupied _ml->template fpfield<24>(_iml)
#define unoccupied_columnindex 24
#define occupied _ml->template fpfield<25>(_iml)
#define occupied_columnindex 25
#define tsyn _ml->template fpfield<26>(_iml)
#define tsyn_columnindex 26
#define u _ml->template fpfield<27>(_iml)
#define u_columnindex 27
#define next_delay _ml->template fpfield<28>(_iml)
#define next_delay_columnindex 28
#define A_GABAA _ml->template fpfield<29>(_iml)
#define A_GABAA_columnindex 29
#define B_GABAA _ml->template fpfield<30>(_iml)
#define B_GABAA_columnindex 30
#define A_GABAB _ml->template fpfield<31>(_iml)
#define A_GABAB_columnindex 31
#define B_GABAB _ml->template fpfield<32>(_iml)
#define B_GABAB_columnindex 32
#define factor_GABAA _ml->template fpfield<33>(_iml)
#define factor_GABAA_columnindex 33
#define factor_GABAB _ml->template fpfield<34>(_iml)
#define factor_GABAB_columnindex 34
#define usingR123 _ml->template fpfield<35>(_iml)
#define usingR123_columnindex 35
#define DA_GABAA _ml->template fpfield<36>(_iml)
#define DA_GABAA_columnindex 36
#define DB_GABAA _ml->template fpfield<37>(_iml)
#define DB_GABAA_columnindex 37
#define DA_GABAB _ml->template fpfield<38>(_iml)
#define DA_GABAB_columnindex 38
#define DB_GABAB _ml->template fpfield<39>(_iml)
#define DB_GABAB_columnindex 39
#define v _ml->template fpfield<40>(_iml)
#define v_columnindex 40
#define _g _ml->template fpfield<41>(_iml)
#define _g_columnindex 41
#define _tsav _ml->template fpfield<42>(_iml)
#define _tsav_columnindex 42
#define _nd_area *_ml->dptr_field<0>(_iml)
#define rng	*_ppvar[2].get<double*>()
#define _p_rng _ppvar[2].literal_value<void*>()
#define delay_times	*_ppvar[3].get<double*>()
#define _p_delay_times _ppvar[3].literal_value<void*>()
#define delay_weights	*_ppvar[4].get<double*>()
#define _p_delay_weights _ppvar[4].literal_value<void*>()
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  2;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_bbsavestate(void*);
 static double _hoc_clearRNG(void*);
 static double _hoc_setRNG(void*);
 static double _hoc_state(void*);
 static double _hoc_setup_delay_vecs(void*);
 static double _hoc_toggleVerbose(void*);
 static double _hoc_urand(void*);
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
 {"bbsavestate", _hoc_bbsavestate},
 {"clearRNG", _hoc_clearRNG},
 {"setRNG", _hoc_setRNG},
 {"state", _hoc_state},
 {"setup_delay_vecs", _hoc_setup_delay_vecs},
 {"toggleVerbose", _hoc_toggleVerbose},
 {"urand", _hoc_urand},
 {0, 0}
};
#define bbsavestate bbsavestate_ProbGABAAB_EMS
#define toggleVerbose toggleVerbose_ProbGABAAB_EMS
#define urand urand_ProbGABAAB_EMS
 extern double bbsavestate( _internalthreadargsproto_ );
 extern double toggleVerbose( _internalthreadargsproto_ );
 extern double urand( _internalthreadargsproto_ );
 /* declare global and static user variables */
#define gmax gmax_ProbGABAAB_EMS
 double gmax = 0.001;
#define init_depleted init_depleted_ProbGABAAB_EMS
 double init_depleted = 0;
#define minis_single_vesicle minis_single_vesicle_ProbGABAAB_EMS
 double minis_single_vesicle = 0;
#define nc_type_param nc_type_param_ProbGABAAB_EMS
 double nc_type_param = 4;
#define tau_d_GABAB tau_d_GABAB_ProbGABAAB_EMS
 double tau_d_GABAB = 260.9;
#define tau_r_GABAB tau_r_GABAB_ProbGABAAB_EMS
 double tau_r_GABAB = 3.5;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"tau_r_GABAB_ProbGABAAB_EMS", "ms"},
 {"tau_d_GABAB_ProbGABAAB_EMS", "ms"},
 {"gmax_ProbGABAAB_EMS", "uS"},
 {"tau_r_GABAA", "ms"},
 {"tau_d_GABAA", "ms"},
 {"Use", "1"},
 {"Dep", "ms"},
 {"Fac", "ms"},
 {"e_GABAA", "mV"},
 {"e_GABAB", "mV"},
 {"Nrrp", "1"},
 {"GABAB_ratio", "1"},
 {"i", "nA"},
 {"i_GABAA", "nA"},
 {"i_GABAB", "nA"},
 {"g_GABAA", "uS"},
 {"g_GABAB", "uS"},
 {"g", "uS"},
 {"unoccupied", "1"},
 {"occupied", "1"},
 {"tsyn", "ms"},
 {"u", "1"},
 {"next_delay", "ms"},
 {0, 0}
};
 static double A_GABAB0 = 0;
 static double A_GABAA0 = 0;
 static double B_GABAB0 = 0;
 static double B_GABAA0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"tau_r_GABAB_ProbGABAAB_EMS", &tau_r_GABAB_ProbGABAAB_EMS},
 {"tau_d_GABAB_ProbGABAAB_EMS", &tau_d_GABAB_ProbGABAAB_EMS},
 {"gmax_ProbGABAAB_EMS", &gmax_ProbGABAAB_EMS},
 {"nc_type_param_ProbGABAAB_EMS", &nc_type_param_ProbGABAAB_EMS},
 {"minis_single_vesicle_ProbGABAAB_EMS", &minis_single_vesicle_ProbGABAAB_EMS},
 {"init_depleted_ProbGABAAB_EMS", &init_depleted_ProbGABAAB_EMS},
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
"ProbGABAAB_EMS",
 "tau_r_GABAA",
 "tau_d_GABAA",
 "Use",
 "Dep",
 "Fac",
 "e_GABAA",
 "e_GABAB",
 "u0",
 "Nrrp",
 "synapseID",
 "verboseLevel",
 "selected_for_report",
 "GABAB_ratio",
 "conductance",
 0,
 "i",
 "i_GABAA",
 "i_GABAB",
 "g_GABAA",
 "g_GABAB",
 "g",
 "A_GABAA_step",
 "B_GABAA_step",
 "A_GABAB_step",
 "B_GABAB_step",
 "unoccupied",
 "occupied",
 "tsyn",
 "u",
 "next_delay",
 0,
 "A_GABAA",
 "B_GABAA",
 "A_GABAB",
 "B_GABAB",
 0,
 "rng",
 "delay_times",
 "delay_weights",
 0};
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
  if (nrn_point_prop_) {
    _nrn_mechanism_access_alloc_seq(_prop) = _nrn_mechanism_access_alloc_seq(nrn_point_prop_);
    _ppvar = _nrn_mechanism_access_dparam(nrn_point_prop_);
  } else {
   _ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 43);
 	/*initialize range parameters*/
 	tau_r_GABAA = 0.2;
 	tau_d_GABAA = 8;
 	Use = 1;
 	Dep = 100;
 	Fac = 10;
 	e_GABAA = -80;
 	e_GABAB = -97;
 	u0 = 0;
 	Nrrp = 1;
 	synapseID = 0;
 	verboseLevel = 0;
 	selected_for_report = 0;
 	GABAB_ratio = 0;
 	conductance = 0;
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 43);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[5])
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 static void bbcore_write(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_write(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 static void bbcore_read(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _ProbGABAAB_EMS_reg() {
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
                                       _nrn_mechanism_field<double>{"tau_r_GABAA"} /* 0 */,
                                       _nrn_mechanism_field<double>{"tau_d_GABAA"} /* 1 */,
                                       _nrn_mechanism_field<double>{"Use"} /* 2 */,
                                       _nrn_mechanism_field<double>{"Dep"} /* 3 */,
                                       _nrn_mechanism_field<double>{"Fac"} /* 4 */,
                                       _nrn_mechanism_field<double>{"e_GABAA"} /* 5 */,
                                       _nrn_mechanism_field<double>{"e_GABAB"} /* 6 */,
                                       _nrn_mechanism_field<double>{"u0"} /* 7 */,
                                       _nrn_mechanism_field<double>{"Nrrp"} /* 8 */,
                                       _nrn_mechanism_field<double>{"synapseID"} /* 9 */,
                                       _nrn_mechanism_field<double>{"verboseLevel"} /* 10 */,
                                       _nrn_mechanism_field<double>{"selected_for_report"} /* 11 */,
                                       _nrn_mechanism_field<double>{"GABAB_ratio"} /* 12 */,
                                       _nrn_mechanism_field<double>{"conductance"} /* 13 */,
                                       _nrn_mechanism_field<double>{"i"} /* 14 */,
                                       _nrn_mechanism_field<double>{"i_GABAA"} /* 15 */,
                                       _nrn_mechanism_field<double>{"i_GABAB"} /* 16 */,
                                       _nrn_mechanism_field<double>{"g_GABAA"} /* 17 */,
                                       _nrn_mechanism_field<double>{"g_GABAB"} /* 18 */,
                                       _nrn_mechanism_field<double>{"g"} /* 19 */,
                                       _nrn_mechanism_field<double>{"A_GABAA_step"} /* 20 */,
                                       _nrn_mechanism_field<double>{"B_GABAA_step"} /* 21 */,
                                       _nrn_mechanism_field<double>{"A_GABAB_step"} /* 22 */,
                                       _nrn_mechanism_field<double>{"B_GABAB_step"} /* 23 */,
                                       _nrn_mechanism_field<double>{"unoccupied"} /* 24 */,
                                       _nrn_mechanism_field<double>{"occupied"} /* 25 */,
                                       _nrn_mechanism_field<double>{"tsyn"} /* 26 */,
                                       _nrn_mechanism_field<double>{"u"} /* 27 */,
                                       _nrn_mechanism_field<double>{"next_delay"} /* 28 */,
                                       _nrn_mechanism_field<double>{"A_GABAA"} /* 29 */,
                                       _nrn_mechanism_field<double>{"B_GABAA"} /* 30 */,
                                       _nrn_mechanism_field<double>{"A_GABAB"} /* 31 */,
                                       _nrn_mechanism_field<double>{"B_GABAB"} /* 32 */,
                                       _nrn_mechanism_field<double>{"factor_GABAA"} /* 33 */,
                                       _nrn_mechanism_field<double>{"factor_GABAB"} /* 34 */,
                                       _nrn_mechanism_field<double>{"usingR123"} /* 35 */,
                                       _nrn_mechanism_field<double>{"DA_GABAA"} /* 36 */,
                                       _nrn_mechanism_field<double>{"DB_GABAA"} /* 37 */,
                                       _nrn_mechanism_field<double>{"DA_GABAB"} /* 38 */,
                                       _nrn_mechanism_field<double>{"DB_GABAB"} /* 39 */,
                                       _nrn_mechanism_field<double>{"v"} /* 40 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 41 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 42 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"rng", "bbcorepointer"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"delay_times", "bbcorepointer"} /* 3 */,
                                       _nrn_mechanism_field<double*>{"delay_weights", "bbcorepointer"} /* 4 */,
                                       _nrn_mechanism_field<void*>{"_tqitem", "netsend"} /* 5 */);
  hoc_register_prop_size(_mechtype, 43, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 4, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 5, "netsend");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ProbGABAAB_EMS /mnt/mydata/cerebellum/mod/ProbGABAAB_EMS.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "GABAAB receptor with presynaptic short-term plasticity";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int clearRNG(_internalthreadargsproto_);
static int setRNG(_internalthreadargsproto_);
static int state(_internalthreadargsproto_);
static int setup_delay_vecs(_internalthreadargsproto_);
 
/*VERBATIM*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#ifndef NRN_VERSION_GTEQ_8_2_0
#include "nrnran123.h"

#ifndef CORENEURON_BUILD
extern int ifarg(int iarg);

extern void* vector_arg(int iarg);
extern double* vector_vec(void* vv);
extern int vector_capacity(void* vv);
#endif

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif

 
static int  setup_delay_vecs ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
    void** vv_delay_times = (void**)(&_p_delay_times);
    void** vv_delay_weights = (void**)(&_p_delay_weights);
    *vv_delay_times = (void*)NULL;
    *vv_delay_weights = (void*)NULL;
    if (ifarg(1)) {
        *vv_delay_times = vector_arg(1);
    }
    if (ifarg(2)) {
        *vv_delay_weights = vector_arg(2);
    }
#endif
  return 0; }
 
static double _hoc_setup_delay_vecs(void* _vptr) {
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
 setup_delay_vecs ( _threadargs_ );
 return(_r);
}
 
static int  state ( _internalthreadargsproto_ ) {
   A_GABAA = A_GABAA * A_GABAA_step ;
   B_GABAA = B_GABAA * B_GABAA_step ;
   A_GABAB = A_GABAB * A_GABAB_step ;
   B_GABAB = B_GABAB * B_GABAB_step ;
    return 0; }
 
static double _hoc_state(void* _vptr) {
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
 state ( _threadargs_ );
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
   double _lresult , _lves , _loccu ;
 _args[1] = _args[0] ;
   _args[2] = _args[0] * GABAB_ratio ;
   if ( _lflag  == 1.0 ) {
     
/*VERBATIM*/
        // setup self events for delayed connections to change weights
        IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
        if (vv_delay_weights && vector_capacity(vv_delay_weights)>=next_delay) {
            double* weights_v = vector_vec(vv_delay_weights);
            double next_delay_weight = weights_v[(int)next_delay];
 _args[0] = conductance * next_delay_weight ;
     next_delay = next_delay + 1.0 ;
     
/*VERBATIM*/
        }
        return;
 }
   if ( _args[0] <= 0.0  || t < 0.0 ) {
     
/*VERBATIM*/
        return;
 }
   if ( Fac > 0.0 ) {
     u = u * exp ( - ( t - tsyn ) / Fac ) ;
     }
   else {
     u = Use ;
     }
   if ( Fac > 0.0 ) {
     u = u + Use * ( 1.0 - u ) ;
     }
   {int  _lcounter ;for ( _lcounter = 0 ; _lcounter <= ( ((int) unoccupied ) - 1 ) ; _lcounter ++ ) {
     _args[3] = exp ( - ( t - tsyn ) / Dep ) ;
     _lresult = urand ( _threadargs_ ) ;
     if ( _lresult > _args[3] ) {
       occupied = occupied + 1.0 ;
       if ( verboseLevel > 0.0 ) {
          printf ( "[Syn %.0f] Recovered! t = %g, Psurv = %g, urand = %g\n" , synapseID , t , _args[3] , _lresult ) ;
          }
       }
     } }
   _lves = 0.0 ;
   _loccu = occupied ;
   if ( _loccu > 1.0  && minis_single_vesicle  && _args[4]  == 1.0 ) {
     _loccu = 1.0 ;
     }
   {int  _lcounter ;for ( _lcounter = 0 ; _lcounter <= ( ((int) _loccu ) - 1 ) ; _lcounter ++ ) {
     _lresult = urand ( _threadargs_ ) ;
     if ( _lresult < u ) {
       occupied = occupied - 1.0 ;
       _lves = _lves + 1.0 ;
       }
     } }
   unoccupied = Nrrp - occupied ;
   tsyn = t ;
   if ( _lves > 0.0 ) {
     A_GABAA = A_GABAA + _lves / Nrrp * _args[1] * factor_GABAA ;
     B_GABAA = B_GABAA + _lves / Nrrp * _args[1] * factor_GABAA ;
     A_GABAB = A_GABAB + _lves / Nrrp * _args[2] * factor_GABAB ;
     B_GABAB = B_GABAB + _lves / Nrrp * _args[2] * factor_GABAB ;
     if ( verboseLevel > 0.0 ) {
        printf ( "[Syn %.0f] Release! t = %g, vals: %g %g %g %g\n" , synapseID , t , A_GABAA , _args[1] , factor_GABAA , _args[0] ) ;
        }
     }
   else {
     if ( verboseLevel > 0.0 ) {
        printf ( "[Syn %.0f] Failure! t = %g: urand = %g\n" , synapseID , t , _lresult ) ;
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
 if ( _args[4]  == 0.0 ) {
     
/*VERBATIM*/
            // setup self events for delayed connections to change weights
            IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
            IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

            if (vv_delay_times && vector_capacity(vv_delay_times)>=1) {
                double* deltm_el = vector_vec(vv_delay_times);
                int delay_times_idx;
                next_delay = 0;
                for (delay_times_idx = 0; delay_times_idx < vector_capacity(vv_delay_times); ++delay_times_idx) {
                    double next_delay_t = deltm_el[delay_times_idx];
 net_send ( _tqitem, _args, _pnt, t +  next_delay_t , 1.0 ) ;
     
/*VERBATIM*/
                }
            }
 }
   }
 
static int  setRNG ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
    #ifndef CORENEURON_BUILD
    // For compatibility, allow for either MCellRan4 or Random123
    // Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        uint32_t a3 = 0;

        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {   // not a double, so assume hoc object type
        void** pv = (void**)(&_p_rng);
        *pv = nrn_random_arg(1);
    } else {  // no arg, so clear pointer
        void** pv = (void**)(&_p_rng);
        *pv = (void*)0;
    }
    #endif
  return 0; }
 
static double _hoc_setRNG(void* _vptr) {
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
 setRNG ( _threadargs_ );
 return(_r);
}
 
static int  clearRNG ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
    #ifndef CORENEURON_BUILD
    if (usingR123) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
    } else {
        void** pv = (void**)(&_p_rng);
        if (*pv) {
            *pv = (void*)0;
        }
    }
    #endif
  return 0; }
 
static double _hoc_clearRNG(void* _vptr) {
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
 clearRNG ( _threadargs_ );
 return(_r);
}
 
double urand ( _internalthreadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
    double value = 0.0;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(RANDCAST _p_rng);
        #endif
    } else {
        // Note: prior versions used scop_random(1), but since we never use this model without configuring the rng.  Maybe should throw error?
        value = 0.0;
    }
    _lurand = value;
 
return _lurand;
 }
 
static double _hoc_urand(void* _vptr) {
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
 _r =  urand ( _threadargs_ );
 return(_r);
}
 
double bbsavestate ( _internalthreadargsproto_ ) {
   double _lbbsavestate;
 _lbbsavestate = 0.0 ;
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval;
#ifndef NRN_VERSION_GTEQ_8_2_0
        double *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
#endif
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng) {
            // tell how many items need saving
            if (*xdir == -1) {  // count items
                if( usingR123 ) {
                    *xdir = 2.0;
                } else {
                    *xdir = 1.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  // save
                if( usingR123 ) {
                    uint32_t seq;
                    char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(RANDCAST _p_rng);
                }
            } else {  // restore
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(RANDCAST _p_rng, (long)(xval[0]));
                }
            }
        }
#endif
 
return _lbbsavestate;
 }
 
static double _hoc_bbsavestate(void* _vptr) {
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
 _r =  bbsavestate ( _threadargs_ );
 return(_r);
}
 
double toggleVerbose ( _internalthreadargsproto_ ) {
   double _ltoggleVerbose;
 verboseLevel = 1.0 - verboseLevel ;
   
return _ltoggleVerbose;
 }
 
static double _hoc_toggleVerbose(void* _vptr) {
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
 _r =  toggleVerbose ( _threadargs_ );
 return(_r);
}
 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {
  IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
  IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

  if (d) {
    uint32_t* di = ((uint32_t*)d) + *d_offset;
    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
    nrnran123_getids3(*pv, di, di+1, di+2);

    // write strem sequence
    char which;
    nrnran123_getseq(*pv, di+3, &which);
    di[4] = (int)which;
    //printf("ProbGABAAB_EMS bbcore_write %d %d %d\n", di[0], di[1], di[2]);

  }
  // reserve random123 parameters on serialization buffer
  *d_offset += 5;

  // serialize connection delay vectors
  if (vv_delay_times && vv_delay_weights &&
     (vector_capacity(vv_delay_times) >= 1) && (vector_capacity(vv_delay_weights) >= 1)) {
    if (d) {
      uint32_t* di = ((uint32_t*)d) + *d_offset;
      // store vector sizes for deserialization
      di[0] = vector_capacity(vv_delay_times);
      di[1] = vector_capacity(vv_delay_weights);
    }
    if (x) {
      double* delay_times_el = vector_vec(vv_delay_times);
      double* delay_weights_el = vector_vec(vv_delay_weights);
      double* x_i = x + *x_offset;
      int delay_vecs_idx;
      int x_idx = 0;
      for(delay_vecs_idx = 0; delay_vecs_idx < vector_capacity(vv_delay_times); ++delay_vecs_idx) {
         x_i[x_idx++] = delay_times_el[delay_vecs_idx];
         x_i[x_idx++] = delay_weights_el[delay_vecs_idx];
      }
    }
    // reserve space for connection delay data on serialization buffer
    *x_offset += vector_capacity(vv_delay_times) + vector_capacity(vv_delay_weights);
  } else {
    if (d) {
      uint32_t* di = ((uint32_t*)d) + *d_offset;
      di[0] = 0;
      di[1] = 0;
    }

  }
  // reserve space for delay vectors (may be 0)
  *d_offset += 2;

}

static void bbcore_read(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {
  // deserialize random123 data
  uint32_t* di = ((uint32_t*)d) + *d_offset;
  if (di[0] != 0 || di[1] != 0 || di[2] != 0) {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
#if !NRNBBCORE
      if(*pv) {
          nrnran123_deletestream(*pv);
      }
#endif
      *pv = nrnran123_newstream3(di[0], di[1], di[2]);
      // restore stream sequence
      char which = (char)di[4];
      nrnran123_setseq(*pv, di[3], which);
  }
  //printf("ProbAMPANMDA_EMS bbcore_read %d %d %d\n", di[0], di[1], di[2]);

  int delay_times_sz = di[5];
  int delay_weights_sz = di[6];
  *d_offset += 7;

  if ((delay_times_sz > 0) && (delay_weights_sz > 0)) {
    double* x_i = x + *x_offset;

    // allocate vectors
    if (!_p_delay_times) {
      _p_delay_times = (double*)vector_new1(delay_times_sz);
    }
    assert(delay_times_sz == vector_capacity((IvocVect*)_p_delay_times));
    if (!_p_delay_weights) {
      _p_delay_weights = (double*)vector_new1(delay_weights_sz);
    }
    assert(delay_weights_sz == vector_capacity((IvocVect*)_p_delay_weights));

    double* delay_times_el = vector_vec((IvocVect*)_p_delay_times);
    double* delay_weights_el = vector_vec((IvocVect*)_p_delay_weights);

    // copy data
    int x_idx;
    int vec_idx = 0;
    for(x_idx = 0; x_idx < delay_times_sz + delay_weights_sz; x_idx += 2) {
      delay_times_el[vec_idx] = x_i[x_idx];
      delay_weights_el[vec_idx++] = x_i[x_idx+1];
    }
    *x_offset += delay_times_sz + delay_weights_sz;

  }
}
 
static int _ode_count(int _type){ hoc_execerror("ProbGABAAB_EMS", "cannot be used with CVODE"); return 0;}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  A_GABAB = A_GABAB0;
  A_GABAA = A_GABAA0;
  B_GABAB = B_GABAB0;
  B_GABAA = B_GABAA0;
 {
   double _ltp_GABAA , _ltp_GABAB ;
 tsyn = 0.0 ;
   u = u0 ;
   if ( init_depleted ) {
     unoccupied = Nrrp ;
     occupied = 0.0 ;
     }
   else {
     unoccupied = 0.0 ;
     occupied = Nrrp ;
     }
   A_GABAA = 0.0 ;
   B_GABAA = 0.0 ;
   A_GABAB = 0.0 ;
   B_GABAB = 0.0 ;
   _ltp_GABAA = ( tau_r_GABAA * tau_d_GABAA ) / ( tau_d_GABAA - tau_r_GABAA ) * log ( tau_d_GABAA / tau_r_GABAA ) ;
   _ltp_GABAB = ( tau_r_GABAB * tau_d_GABAB ) / ( tau_d_GABAB - tau_r_GABAB ) * log ( tau_d_GABAB / tau_r_GABAB ) ;
   factor_GABAA = - exp ( - _ltp_GABAA / tau_r_GABAA ) + exp ( - _ltp_GABAA / tau_d_GABAA ) ;
   factor_GABAA = 1.0 / factor_GABAA ;
   factor_GABAB = - exp ( - _ltp_GABAB / tau_r_GABAB ) + exp ( - _ltp_GABAB / tau_d_GABAB ) ;
   factor_GABAB = 1.0 / factor_GABAB ;
   A_GABAA_step = exp ( dt * ( ( - 1.0 ) / tau_r_GABAA ) ) ;
   B_GABAA_step = exp ( dt * ( ( - 1.0 ) / tau_d_GABAA ) ) ;
   A_GABAB_step = exp ( dt * ( ( - 1.0 ) / tau_r_GABAB ) ) ;
   B_GABAB_step = exp ( dt * ( ( - 1.0 ) / tau_d_GABAB ) ) ;
   
/*VERBATIM*/
        if( usingR123 ) {
            nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
        }
 next_delay = - 1.0 ;
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
   g_GABAA = gmax * ( B_GABAA - A_GABAA ) ;
   g_GABAB = gmax * ( B_GABAB - A_GABAB ) ;
   g = g_GABAA + g_GABAB ;
   i_GABAA = g_GABAA * ( v - e_GABAA ) ;
   i_GABAB = g_GABAB * ( v - e_GABAB ) ;
   i = i_GABAA + i_GABAB ;
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
 {  { state(_threadargs_); }
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
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/ProbGABAAB_EMS.mod";
    const char* nmodl_file_text = 
  "COMMENT\n"
  "/**\n"
  " * @file ProbGABAAB.mod\n"
  " * @brief\n"
  " * @author king, muller\n"
  " * @date 2011-08-17\n"
  " * @remark Copyright \n"
  "\n"
  " BBP/EPFL 2005-2011; All rights reserved. Do not distribute without further notice.\n"
  " */\n"
  "ENDCOMMENT\n"
  "\n"
  "TITLE GABAAB receptor with presynaptic short-term plasticity\n"
  "\n"
  "\n"
  "COMMENT\n"
  "GABAA receptor conductance using a dual-exponential profile\n"
  "presynaptic short-term plasticity based on Fuhrmann et al, 2002\n"
  "Implemented by Srikanth Ramaswamy, Blue Brain Project, March 2009\n"
  "\n"
  "_EMS (Eilif Michael Srikanth)\n"
  "Modification of ProbGABAA: 2-State model by Eilif Muller, Michael Reimann, Srikanth Ramaswamy, Blue Brain Project, August 2011\n"
  "This new model was motivated by the following constraints:\n"
  "\n"
  "1) No consumption on failure.\n"
  "2) No release just after release until recovery.\n"
  "3) Same ensemble averaged trace as deterministic/canonical Tsodyks-Markram\n"
  "   using same parameters determined from experiment.\n"
  "4) Same quantal size as present production probabilistic model.\n"
  "\n"
  "To satisfy these constaints, the synapse is implemented as a\n"
  "uni-vesicular (generalization to multi-vesicular should be\n"
  "straight-forward) 2-state Markov process.  The states are\n"
  "{1=recovered, 0=unrecovered}.\n"
  "\n"
  "For a pre-synaptic spike or external spontaneous release trigger\n"
  "event, the synapse will only release if it is in the recovered state,\n"
  "and with probability u (which follows facilitation dynamics).  If it\n"
  "releases, it will transition to the unrecovered state.  Recovery is as\n"
  "a Poisson process with rate 1/Dep.\n"
  "\n"
  "This model satisys all of (1)-(4).\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "NEURON {\n"
  "    THREADSAFE\n"
  "    POINT_PROCESS ProbGABAAB_EMS\n"
  "\n"
  "    RANGE tau_r_GABAA, tau_d_GABAA\n"
  "    GLOBAL tau_r_GABAB, tau_d_GABAB\n"
  "    RANGE Use, u, Dep, Fac, u0, tsyn\n"
  "    RANGE unoccupied, occupied, Nrrp\n"
  "\n"
  "    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB, GABAB_ratio\n"
  "    RANGE A_GABAA_step, B_GABAA_step, A_GABAB_step, B_GABAB_step\n"
  "\n"
  "    NONSPECIFIC_CURRENT i\n"
  "    BBCOREPOINTER rng\n"
  "    RANGE synapseID, selected_for_report, verboseLevel, conductance\n"
  "    RANGE next_delay\n"
  "    BBCOREPOINTER delay_times, delay_weights\n"
  "    GLOBAL nc_type_param\n"
  "    GLOBAL minis_single_vesicle\n"
  "    GLOBAL init_depleted\n"
  "\n"
  "    :RANGE sgid, tgid  : For debugging\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    tau_r_GABAA = 0.2   (ms)  : dual-exponential conductance profile\n"
  "    tau_d_GABAA = 8     (ms)  : IMPORTANT: tau_r < tau_d\n"
  "    tau_r_GABAB = 3.5   (ms)  : dual-exponential conductance profile :Placeholder value from hippocampal recordings SR\n"
  "    tau_d_GABAB = 260.9 (ms)  : IMPORTANT: tau_r < tau_d  :Placeholder value from hippocampal recordings\n"
  "    Use = 1.0           (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values)\n"
  "    Dep = 100           (ms)  : relaxation time constant from depression\n"
  "    Fac = 10            (ms)  : relaxation time constant from facilitation\n"
  "    e_GABAA = -80       (mV)  : GABAA reversal potential\n"
  "    e_GABAB = -97       (mV)  : GABAB reversal potential\n"
  "    gmax = .001         (uS)  : weight conversion factor (from nS to uS)\n"
  "    u0   = 0                  : initial value of u, which is the running value of release probability\n"
  "    Nrrp = 1            (1)   : Number of total release sites for given contact\n"
  "\n"
  "    synapseID = 0\n"
  "    verboseLevel = 0\n"
  "    selected_for_report = 0\n"
  "    GABAB_ratio = 0     (1)   : The ratio of GABAB to GABAA\n"
  "    conductance = 0.0\n"
  "    nc_type_param = 4\n"
  "    minis_single_vesicle = 0   :// 0 - no limit (old behavior)\n"
  "    init_depleted = 0          :// 0 - init full (old behavior)\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1\n"
  "for comparison with Pr to decide whether to activate the synapse or not\n"
  "ENDCOMMENT\n"
  "\n"
  "VERBATIM\n"
  "\n"
  "#include<stdlib.h>\n"
  "#include<stdio.h>\n"
  "#include<math.h>\n"
  "#ifndef NRN_VERSION_GTEQ_8_2_0\n"
  "#include \"nrnran123.h\"\n"
  "\n"
  "#ifndef CORENEURON_BUILD\n"
  "extern int ifarg(int iarg);\n"
  "\n"
  "extern void* vector_arg(int iarg);\n"
  "extern double* vector_vec(void* vv);\n"
  "extern int vector_capacity(void* vv);\n"
  "#endif\n"
  "\n"
  "double nrn_random_pick(void* r);\n"
  "void* nrn_random_arg(int argpos);\n"
  "#define RANDCAST\n"
  "#else\n"
  "#define RANDCAST (Rand*)\n"
  "#endif\n"
  "\n"
  "ENDVERBATIM\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        i (nA)\n"
  "        i_GABAA (nA)\n"
  "        i_GABAB (nA)\n"
  "        g_GABAA (uS)\n"
  "        g_GABAB (uS)\n"
  "        g (uS)\n"
  "        factor_GABAA\n"
  "        factor_GABAB\n"
  "        A_GABAA_step\n"
  "        B_GABAA_step\n"
  "        A_GABAB_step\n"
  "        B_GABAB_step\n"
  "\n"
  "        rng\n"
  "        usingR123            : TEMPORARY until mcellran4 completely deprecated\n"
  "\n"
  "        : MVR\n"
  "        unoccupied (1) : no. of unoccupied sites following release event\n"
  "        occupied   (1) : no. of occupied sites following one epoch of recovery\n"
  "        tsyn (ms) : the time of the last spike\n"
  "        u (1) : running release probability\n"
  "\n"
  "        : stuff for delayed connections\n"
  "        delay_times\n"
  "        delay_weights\n"
  "        next_delay (ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        A_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAA\n"
  "        B_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAA\n"
  "        A_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAB\n"
  "        B_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAB\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL {\n"
  "        LOCAL tp_GABAA, tp_GABAB\n"
  "\n"
  "        tsyn = 0\n"
  "        u=u0\n"
  "\n"
  "        : MVR\n"
  "        if ( init_depleted ) {\n"
  "            unoccupied = Nrrp\n"
  "            occupied = 0\n"
  "        } else {\n"
  "            unoccupied = 0\n"
  "            occupied = Nrrp\n"
  "        }\n"
  "\n"
  "        A_GABAA = 0\n"
  "        B_GABAA = 0\n"
  "\n"
  "        A_GABAB = 0\n"
  "        B_GABAB = 0\n"
  "\n"
  "        tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA) :time to peak of the conductance\n"
  "        tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) :time to peak of the conductance\n"
  "\n"
  "        factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) :GABAA Normalization factor - so that when t = tp_GABAA, gsyn = gpeak\n"
  "        factor_GABAA = 1/factor_GABAA\n"
  "\n"
  "        factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) :GABAB Normalization factor - so that when t = tp_GABAB, gsyn = gpeak\n"
  "        factor_GABAB = 1/factor_GABAB\n"
  "\n"
  "        A_GABAA_step = exp(dt*(( - 1.0 ) / tau_r_GABAA))\n"
  "        B_GABAA_step = exp(dt*(( - 1.0 ) / tau_d_GABAA))\n"
  "        A_GABAB_step = exp(dt*(( - 1.0 ) / tau_r_GABAB))\n"
  "        B_GABAB_step = exp(dt*(( - 1.0 ) / tau_d_GABAB))\n"
  "\n"
  "    VERBATIM\n"
  "        if( usingR123 ) {\n"
  "            nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);\n"
  "        }\n"
  "    ENDVERBATIM\n"
  "\n"
  "        next_delay = -1\n"
  "\n"
  "}\n"
  "\n"
  "PROCEDURE setup_delay_vecs() {\n"
  "VERBATIM\n"
  "#ifndef CORENEURON_BUILD\n"
  "    void** vv_delay_times = (void**)(&_p_delay_times);\n"
  "    void** vv_delay_weights = (void**)(&_p_delay_weights);\n"
  "    *vv_delay_times = (void*)NULL;\n"
  "    *vv_delay_weights = (void*)NULL;\n"
  "    if (ifarg(1)) {\n"
  "        *vv_delay_times = vector_arg(1);\n"
  "    }\n"
  "    if (ifarg(2)) {\n"
  "        *vv_delay_weights = vector_arg(2);\n"
  "    }\n"
  "#endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE state\n"
  "\n"
  "        g_GABAA = gmax*(B_GABAA-A_GABAA) :compute time varying conductance as the difference of state variables B_GABAA and A_GABAA\n"
  "        g_GABAB = gmax*(B_GABAB-A_GABAB) :compute time varying conductance as the difference of state variables B_GABAB and A_GABAB\n"
  "        g = g_GABAA + g_GABAB\n"
  "        i_GABAA = g_GABAA*(v-e_GABAA) :compute the GABAA driving force based on the time varying conductance, membrane potential, and GABAA reversal\n"
  "        i_GABAB = g_GABAB*(v-e_GABAB) :compute the GABAB driving force based on the time varying conductance, membrane potential, and GABAB reversal\n"
  "        i = i_GABAA + i_GABAB\n"
  "}\n"
  "\n"
  "PROCEDURE state() {\n"
  "        A_GABAA = A_GABAA*A_GABAA_step\n"
  "        B_GABAA = B_GABAA*B_GABAA_step\n"
  "        A_GABAB = A_GABAB*A_GABAB_step\n"
  "        B_GABAB = B_GABAB*B_GABAB_step\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "NET_RECEIVE (weight, weight_GABAA, weight_GABAB, Psurv, nc_type) {\n"
  "    : Psurv - survival probability of unrecovered state\n"
  "    : nc_type:\n"
  "    :   0 = presynaptic netcon\n"
  "    :   1 = spontmini netcon\n"
  "    :   2 = replay netcon\n"
  "\n"
  "    LOCAL result, ves, occu\n"
  "    weight_GABAA = weight\n"
  "    weight_GABAB = weight * GABAB_ratio\n"
  "\n"
  "    INITIAL {\n"
  "        if (nc_type == 0) {  :// presynaptic netcon\n"
  "    VERBATIM\n"
  "            // setup self events for delayed connections to change weights\n"
  "            IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));\n"
  "            IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));\n"
  "\n"
  "            if (vv_delay_times && vector_capacity(vv_delay_times)>=1) {\n"
  "                double* deltm_el = vector_vec(vv_delay_times);\n"
  "                int delay_times_idx;\n"
  "                next_delay = 0;\n"
  "                for (delay_times_idx = 0; delay_times_idx < vector_capacity(vv_delay_times); ++delay_times_idx) {\n"
  "                    double next_delay_t = deltm_el[delay_times_idx];\n"
  "    ENDVERBATIM\n"
  "                    net_send(next_delay_t, 1)\n"
  "    VERBATIM\n"
  "                }\n"
  "            }\n"
  "    ENDVERBATIM\n"
  "        }\n"
  "    }\n"
  "\n"
  "    if (flag == 1) {  :// self event to set next weight at\n"
  "        :UNITSOFF\n"
  "        :    printf( \"gaba: self event at synapse %f time %g\\n\", synapseID, t)\n"
  "        :UNITSON\n"
  "    VERBATIM\n"
  "        // setup self events for delayed connections to change weights\n"
  "        IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));\n"
  "        if (vv_delay_weights && vector_capacity(vv_delay_weights)>=next_delay) {\n"
  "            double* weights_v = vector_vec(vv_delay_weights);\n"
  "            double next_delay_weight = weights_v[(int)next_delay];\n"
  "    ENDVERBATIM\n"
  "            weight = conductance * next_delay_weight\n"
  "            next_delay = next_delay + 1\n"
  "    VERBATIM\n"
  "        }\n"
  "        return;\n"
  "    ENDVERBATIM\n"
  "    }\n"
  "\n"
  "    : [flag == 0] Handle a spike which arrived\n"
  "    :UNITSOFF\n"
  "    :printf(\"[Syn %.0f] Received! (%f -> %f) with weight %g at time %g\\n\", synapseID, sgid, tgid, weight, t)\n"
  "    :UNITSON\n"
  "\n"
  "    : Do not perform any calculations if the synapse (netcon) is deactivated. This avoids drawing from\n"
  "    : random number stream. Also, disable in case of t < 0 (in case of ForwardSkip) which causes numerical\n"
  "    : instability if synapses are activated.\n"
  "    if ( weight <= 0 || t < 0 ) {\n"
  "    VERBATIM\n"
  "        return;\n"
  "    ENDVERBATIM\n"
  "    }\n"
  "\n"
  "    : calc u at event-\n"
  "    if (Fac > 0) {\n"
  "        u = u * exp(-(t - tsyn)/Fac)  :// update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.\n"
  "    } else {\n"
  "        u = Use\n"
  "    }\n"
  "    if(Fac > 0){\n"
  "        u = u + Use*(1-u)  :// update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.\n"
  "    }\n"
  "\n"
  "    : recovery\n"
  "    FROM counter = 0 TO (unoccupied - 1) {\n"
  "        : Iterate over all unoccupied sites and compute how many recover\n"
  "        Psurv = exp(-(t-tsyn)/Dep)\n"
  "        result = urand()\n"
  "        if (result>Psurv) {\n"
  "            occupied = occupied + 1     :// recover a previously unoccupied site\n"
  "            if ( verboseLevel > 0 ) {\n"
  "                UNITSOFF\n"
  "                printf(\"[Syn %.0f] Recovered! t = %g, Psurv = %g, urand = %g\\n\", synapseID, t, Psurv, result)\n"
  "                UNITSON\n"
  "            }\n"
  "        }\n"
  "    }\n"
  "\n"
  "    ves = 0                  :// Initialize the number of released vesicles to 0\n"
  "    occu = occupied          :// Make a copy, so we can update occupied in the loop\n"
  "    if (occu > 1 && minis_single_vesicle && nc_type == 1) {    : // if nc_type is spont_mini consider single vesicle\n"
  "        occu = 1\n"
  "    }\n"
  "    FROM counter = 0 TO (occu - 1) {\n"
  "        : iterate over all occupied sites and compute how many release\n"
  "        result = urand()\n"
  "        if (result < u) {\n"
  "            : release a single site!\n"
  "            occupied = occupied - 1  :// decrease the number of occupied sites by 1\n"
  "            ves = ves + 1            :// increase number of relesed vesicles by 1\n"
  "        }\n"
  "    }\n"
  "\n"
  "    : Update number of unoccupied sites\n"
  "    unoccupied = Nrrp - occupied\n"
  "\n"
  "    : Update tsyn\n"
  "    : tsyn knows about all spikes, not only those that released\n"
  "    : i.e. each spike can increase the u, regardless of recovered state.\n"
  "    :      and each spike trigger an evaluation of recovery\n"
  "    tsyn = t\n"
  "\n"
  "    if (ves > 0) { :no need to evaluate unless we have vesicle release\n"
  "        A_GABAA = A_GABAA + ves/Nrrp*weight_GABAA*factor_GABAA\n"
  "        B_GABAA = B_GABAA + ves/Nrrp*weight_GABAA*factor_GABAA\n"
  "        A_GABAB = A_GABAB + ves/Nrrp*weight_GABAB*factor_GABAB\n"
  "        B_GABAB = B_GABAB + ves/Nrrp*weight_GABAB*factor_GABAB\n"
  "\n"
  "        if ( verboseLevel > 0 ) {\n"
  "            UNITSOFF\n"
  "            printf(\"[Syn %.0f] Release! t = %g, vals: %g %g %g %g\\n\",\n"
  "                   synapseID, t, A_GABAA, weight_GABAA, factor_GABAA, weight)\n"
  "            UNITSON\n"
  "        }\n"
  "\n"
  "    } else {\n"
  "        : total release failure\n"
  "        if ( verboseLevel > 0 ) {\n"
  "            UNITSOFF\n"
  "            printf(\"[Syn %.0f] Failure! t = %g: urand = %g\\n\", synapseID, t, result)\n"
  "            UNITSON\n"
  "        }\n"
  "    }\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE setRNG() {\n"
  "VERBATIM\n"
  "    #ifndef CORENEURON_BUILD\n"
  "    // For compatibility, allow for either MCellRan4 or Random123\n"
  "    // Distinguish by the arg types\n"
  "    // Object => MCellRan4, seeds (double) => Random123\n"
  "    usingR123 = 0;\n"
  "    if( ifarg(1) && hoc_is_double_arg(1) ) {\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "        uint32_t a2 = 0;\n"
  "        uint32_t a3 = 0;\n"
  "\n"
  "        if (*pv) {\n"
  "            nrnran123_deletestream(*pv);\n"
  "            *pv = (nrnran123_State*)0;\n"
  "        }\n"
  "        if (ifarg(2)) {\n"
  "            a2 = (uint32_t)*getarg(2);\n"
  "        }\n"
  "        if (ifarg(3)) {\n"
  "            a3 = (uint32_t)*getarg(3);\n"
  "        }\n"
  "        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);\n"
  "        usingR123 = 1;\n"
  "    } else if( ifarg(1) ) {   // not a double, so assume hoc object type\n"
  "        void** pv = (void**)(&_p_rng);\n"
  "        *pv = nrn_random_arg(1);\n"
  "    } else {  // no arg, so clear pointer\n"
  "        void** pv = (void**)(&_p_rng);\n"
  "        *pv = (void*)0;\n"
  "    }\n"
  "    #endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE clearRNG() {\n"
  "VERBATIM\n"
  "    #ifndef CORENEURON_BUILD\n"
  "    if (usingR123) {\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "        if (*pv) {\n"
  "            nrnran123_deletestream(*pv);\n"
  "            *pv = (nrnran123_State*)0;\n"
  "        }\n"
  "    } else {\n"
  "        void** pv = (void**)(&_p_rng);\n"
  "        if (*pv) {\n"
  "            *pv = (void*)0;\n"
  "        }\n"
  "    }\n"
  "    #endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION urand() {\n"
  "VERBATIM\n"
  "    double value = 0.0;\n"
  "    if ( usingR123 ) {\n"
  "        value = nrnran123_dblpick((nrnran123_State*)_p_rng);\n"
  "    } else if (_p_rng) {\n"
  "        #ifndef CORENEURON_BUILD\n"
  "        value = nrn_random_pick(RANDCAST _p_rng);\n"
  "        #endif\n"
  "    } else {\n"
  "        // Note: prior versions used scop_random(1), but since we never use this model without configuring the rng.  Maybe should throw error?\n"
  "        value = 0.0;\n"
  "    }\n"
  "    _lurand = value;\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION bbsavestate() {\n"
  "        bbsavestate = 0\n"
  "VERBATIM\n"
  "#ifndef CORENEURON_BUILD\n"
  "        /* first arg is direction (0 save, 1 restore), second is array*/\n"
  "        /* if first arg is -1, fill xdir with the size of the array */\n"
  "        double *xdir, *xval;\n"
  "#ifndef NRN_VERSION_GTEQ_8_2_0\n"
  "        double *hoc_pgetarg();\n"
  "        long nrn_get_random_sequence(void* r);\n"
  "        void nrn_set_random_sequence(void* r, int val);\n"
  "#endif\n"
  "        xdir = hoc_pgetarg(1);\n"
  "        xval = hoc_pgetarg(2);\n"
  "        if (_p_rng) {\n"
  "            // tell how many items need saving\n"
  "            if (*xdir == -1) {  // count items\n"
  "                if( usingR123 ) {\n"
  "                    *xdir = 2.0;\n"
  "                } else {\n"
  "                    *xdir = 1.0;\n"
  "                }\n"
  "                return 0.0;\n"
  "            } else if(*xdir ==0 ) {  // save\n"
  "                if( usingR123 ) {\n"
  "                    uint32_t seq;\n"
  "                    char which;\n"
  "                    nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );\n"
  "                    xval[0] = (double) seq;\n"
  "                    xval[1] = (double) which;\n"
  "                } else {\n"
  "                    xval[0] = (double)nrn_get_random_sequence(RANDCAST _p_rng);\n"
  "                }\n"
  "            } else {  // restore\n"
  "                if( usingR123 ) {\n"
  "                    nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );\n"
  "                } else {\n"
  "                    nrn_set_random_sequence(RANDCAST _p_rng, (long)(xval[0]));\n"
  "                }\n"
  "            }\n"
  "        }\n"
  "#endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION toggleVerbose() {\n"
  "    verboseLevel = 1 - verboseLevel\n"
  "}\n"
  "\n"
  "\n"
  "VERBATIM\n"
  "static void bbcore_write(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {\n"
  "  IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));\n"
  "  IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));\n"
  "\n"
  "  if (d) {\n"
  "    uint32_t* di = ((uint32_t*)d) + *d_offset;\n"
  "    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "    nrnran123_getids3(*pv, di, di+1, di+2);\n"
  "\n"
  "    // write strem sequence\n"
  "    char which;\n"
  "    nrnran123_getseq(*pv, di+3, &which);\n"
  "    di[4] = (int)which;\n"
  "    //printf(\"ProbGABAAB_EMS bbcore_write %d %d %d\\n\", di[0], di[1], di[2]);\n"
  "\n"
  "  }\n"
  "  // reserve random123 parameters on serialization buffer\n"
  "  *d_offset += 5;\n"
  "\n"
  "  // serialize connection delay vectors\n"
  "  if (vv_delay_times && vv_delay_weights &&\n"
  "     (vector_capacity(vv_delay_times) >= 1) && (vector_capacity(vv_delay_weights) >= 1)) {\n"
  "    if (d) {\n"
  "      uint32_t* di = ((uint32_t*)d) + *d_offset;\n"
  "      // store vector sizes for deserialization\n"
  "      di[0] = vector_capacity(vv_delay_times);\n"
  "      di[1] = vector_capacity(vv_delay_weights);\n"
  "    }\n"
  "    if (x) {\n"
  "      double* delay_times_el = vector_vec(vv_delay_times);\n"
  "      double* delay_weights_el = vector_vec(vv_delay_weights);\n"
  "      double* x_i = x + *x_offset;\n"
  "      int delay_vecs_idx;\n"
  "      int x_idx = 0;\n"
  "      for(delay_vecs_idx = 0; delay_vecs_idx < vector_capacity(vv_delay_times); ++delay_vecs_idx) {\n"
  "         x_i[x_idx++] = delay_times_el[delay_vecs_idx];\n"
  "         x_i[x_idx++] = delay_weights_el[delay_vecs_idx];\n"
  "      }\n"
  "    }\n"
  "    // reserve space for connection delay data on serialization buffer\n"
  "    *x_offset += vector_capacity(vv_delay_times) + vector_capacity(vv_delay_weights);\n"
  "  } else {\n"
  "    if (d) {\n"
  "      uint32_t* di = ((uint32_t*)d) + *d_offset;\n"
  "      di[0] = 0;\n"
  "      di[1] = 0;\n"
  "    }\n"
  "\n"
  "  }\n"
  "  // reserve space for delay vectors (may be 0)\n"
  "  *d_offset += 2;\n"
  "\n"
  "}\n"
  "\n"
  "static void bbcore_read(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {\n"
  "  // deserialize random123 data\n"
  "  uint32_t* di = ((uint32_t*)d) + *d_offset;\n"
  "  if (di[0] != 0 || di[1] != 0 || di[2] != 0) {\n"
  "      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "#if !NRNBBCORE\n"
  "      if(*pv) {\n"
  "          nrnran123_deletestream(*pv);\n"
  "      }\n"
  "#endif\n"
  "      *pv = nrnran123_newstream3(di[0], di[1], di[2]);\n"
  "      // restore stream sequence\n"
  "      char which = (char)di[4];\n"
  "      nrnran123_setseq(*pv, di[3], which);\n"
  "  }\n"
  "  //printf(\"ProbAMPANMDA_EMS bbcore_read %d %d %d\\n\", di[0], di[1], di[2]);\n"
  "\n"
  "  int delay_times_sz = di[5];\n"
  "  int delay_weights_sz = di[6];\n"
  "  *d_offset += 7;\n"
  "\n"
  "  if ((delay_times_sz > 0) && (delay_weights_sz > 0)) {\n"
  "    double* x_i = x + *x_offset;\n"
  "\n"
  "    // allocate vectors\n"
  "    if (!_p_delay_times) {\n"
  "      _p_delay_times = (double*)vector_new1(delay_times_sz);\n"
  "    }\n"
  "    assert(delay_times_sz == vector_capacity((IvocVect*)_p_delay_times));\n"
  "    if (!_p_delay_weights) {\n"
  "      _p_delay_weights = (double*)vector_new1(delay_weights_sz);\n"
  "    }\n"
  "    assert(delay_weights_sz == vector_capacity((IvocVect*)_p_delay_weights));\n"
  "\n"
  "    double* delay_times_el = vector_vec((IvocVect*)_p_delay_times);\n"
  "    double* delay_weights_el = vector_vec((IvocVect*)_p_delay_weights);\n"
  "\n"
  "    // copy data\n"
  "    int x_idx;\n"
  "    int vec_idx = 0;\n"
  "    for(x_idx = 0; x_idx < delay_times_sz + delay_weights_sz; x_idx += 2) {\n"
  "      delay_times_el[vec_idx] = x_i[x_idx];\n"
  "      delay_weights_el[vec_idx++] = x_i[x_idx+1];\n"
  "    }\n"
  "    *x_offset += delay_times_sz + delay_weights_sz;\n"
  "\n"
  "  }\n"
  "}\n"
  "ENDVERBATIM\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
