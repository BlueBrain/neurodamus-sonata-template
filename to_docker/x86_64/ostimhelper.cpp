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
static constexpr auto number_of_floating_point_variables = 6;
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
 
#define nrn_init _nrn_init__OdorStimHelper
#define _nrn_initial _nrn_initial__OdorStimHelper
#define nrn_cur _nrn_cur__OdorStimHelper
#define _nrn_current _nrn_current__OdorStimHelper
#define nrn_jacob _nrn_jacob__OdorStimHelper
#define nrn_state _nrn_state__OdorStimHelper
#define _net_receive _net_receive__OdorStimHelper 
#define initstream initstream__OdorStimHelper 
#define noiseFromRandom123 noiseFromRandom123__OdorStimHelper 
 
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
#define start _ml->template fpfield<0>(_iml)
#define start_columnindex 0
#define dur _ml->template fpfield<1>(_iml)
#define dur_columnindex 1
#define invl_min _ml->template fpfield<2>(_iml)
#define invl_min_columnindex 2
#define invl_max _ml->template fpfield<3>(_iml)
#define invl_max_columnindex 3
#define v _ml->template fpfield<4>(_iml)
#define v_columnindex 4
#define _tsav _ml->template fpfield<5>(_iml)
#define _tsav_columnindex 5
#define _nd_area *_ml->dptr_field<0>(_iml)
#define space	*_ppvar[2].get<double*>()
#define _p_space _ppvar[2].literal_value<void*>()
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  2;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_invl(void*);
 static double _hoc_initstream(void*);
 static double _hoc_noiseFromRandom123(void*);
 static double _hoc_uniform_pick(void*);
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
 {"invl", _hoc_invl},
 {"initstream", _hoc_initstream},
 {"noiseFromRandom123", _hoc_noiseFromRandom123},
 {"uniform_pick", _hoc_uniform_pick},
 {0, 0}
};
#define invl invl_OdorStimHelper
#define uniform_pick uniform_pick_OdorStimHelper
 extern double invl( _internalthreadargsproto_ );
 extern double uniform_pick( _internalthreadargsproto_ );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"start", "ms"},
 {"dur", "ms"},
 {"invl_min", "ms"},
 {"invl_max", "ms"},
 {0, 0}
};
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
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"OdorStimHelper",
 "start",
 "dur",
 "invl_min",
 "invl_max",
 0,
 0,
 0,
 "space",
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
    assert(_nrn_mechanism_get_num_vars(_prop) == 6);
 	/*initialize range parameters*/
 	start = 0;
 	dur = 0;
 	invl_min = 1;
 	invl_max = 2;
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 6);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[3])
 static void _net_receive(Point_process*, double*, double);
 static void bbcore_write(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_write(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 static void bbcore_read(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _ostimhelper_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nullptr, nullptr, nullptr, nrn_init,
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
                                       _nrn_mechanism_field<double>{"start"} /* 0 */,
                                       _nrn_mechanism_field<double>{"dur"} /* 1 */,
                                       _nrn_mechanism_field<double>{"invl_min"} /* 2 */,
                                       _nrn_mechanism_field<double>{"invl_max"} /* 3 */,
                                       _nrn_mechanism_field<double>{"v"} /* 4 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 5 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"space", "bbcorepointer"} /* 2 */,
                                       _nrn_mechanism_field<void*>{"_tqitem", "netsend"} /* 3 */);
  hoc_register_prop_size(_mechtype, 6, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "netsend");
 add_nrn_artcell(_mechtype, 3);
 add_nrn_has_net_event(_mechtype);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 OdorStimHelper /mnt/mydata/cerebellum/mod/ostimhelper.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int initstream(_internalthreadargsproto_);
static int noiseFromRandom123(_internalthreadargsproto_);
 
double invl ( _internalthreadargsproto_ ) {
   double _linvl;
 _linvl = ( invl_max - invl_min ) * uniform_pick ( _threadargs_ ) + invl_min ;
   
return _linvl;
 }
 
static double _hoc_invl(void* _vptr) {
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
 _r =  invl ( _threadargs_ );
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
   double _ltinv ;
 if ( _lflag  == 1.0 ) {
     net_event ( _pnt, t ) ;
     _ltinv = invl ( _threadargs_ ) ;
     if ( t + _ltinv <= dur ) {
       artcell_net_send ( _tqitem, _args, _pnt, t +  _ltinv , 1.0 ) ;
       }
     }
   } }
 
static int  initstream ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
  if (_p_space) {
    nrnran123_setseq((nrnran123_State*)_p_space, 0, 0);
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
 
double uniform_pick ( _internalthreadargsproto_ ) {
   double _luniform_pick;
 
/*VERBATIM*/
  if (_p_space) {
    _luniform_pick = nrnran123_dblpick((nrnran123_State*)_p_space);
  }else{
    _luniform_pick = 0.5;
  }
 
return _luniform_pick;
 }
 
static double _hoc_uniform_pick(void* _vptr) {
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
 _r =  uniform_pick ( _threadargs_ );
 return(_r);
}
 
static int  noiseFromRandom123 ( _internalthreadargsproto_ ) {
   
/*VERBATIM*/
#if !NRNBBCORE
 {
  nrnran123_State** pv = (nrnran123_State**)(&_p_space);
  if (*pv) {
    nrnran123_deletestream(*pv);
    *pv = (nrnran123_State*)0;
  }
  *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));
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
  if (!_p_space) {
    fprintf(stderr, "OStimHelper: noiseFromRandom123(1,2,3) not called.\n");
    assert(0);
  }
  if (d) {
    uint32_t* di = ((uint32_t*)d) + *offset;
    nrnran123_State** pv = (nrnran123_State**)(&_p_space);
    nrnran123_getids3(*pv, di, di+1, di+2);
  }
#endif
  *offset += 3;
}

static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
  uint32_t* di = ((uint32_t*)d) + *offset;
  nrnran123_State** pv = (nrnran123_State**)(&_p_space);
#if !NRNBBCORE
  if(*pv) {
    nrnran123_deletestream(*pv);
  }
#endif
  *pv = nrnran123_newstream3(di[0], di[1], di[2]);
  *offset += 3;
}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
 {
   initstream ( _threadargs_ ) ;
   if ( t <= start ) {
     artcell_net_send ( _tqitem, nullptr, _ppvar[1].get<Point_process*>(), t +  start , 1.0 ) ;
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
 _tsav = -1e20;
 initmodel(_threadargs_);
}
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{
} return _current;
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
 v=_v;
{
}}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/ostimhelper.mod";
    const char* nmodl_file_text = 
  "NEURON {\n"
  "  THREADSAFE\n"
  "  ARTIFICIAL_CELL OdorStimHelper\n"
  "  RANGE start, dur, invl_min, invl_max\n"
  "  BBCOREPOINTER space\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  start = 0 (ms)\n"
  "  dur = 0 (ms)\n"
  "  invl_min = 1 (ms)\n"
  "  invl_max = 2 (ms)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  space\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "  initstream()\n"
  "  if (t <= start) {\n"
  "    net_send(start, 1)\n"
  "  }\n"
  "}\n"
  "\n"
  "FUNCTION invl()(ms) {\n"
  "  invl = (invl_max - invl_min)*uniform_pick() + invl_min\n"
  "}\n"
  "\n"
  "NET_RECEIVE(dummy) {\n"
  "  LOCAL tinv\n"
  "  if (flag == 1) {\n"
  "    net_event(t)\n"
  "    tinv = invl()\n"
  "    if (t + tinv <= dur) {\n"
  "      net_send(tinv, 1)\n"
  "    }\n"
  "  }\n"
  "}\n"
  "\n"
  ":::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
  "\n"
  "PROCEDURE initstream() {\n"
  "VERBATIM\n"
  "  if (_p_space) {\n"
  "    nrnran123_setseq((nrnran123_State*)_p_space, 0, 0);\n"
  "  }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION uniform_pick() {\n"
  "VERBATIM\n"
  "  if (_p_space) {\n"
  "    _luniform_pick = nrnran123_dblpick((nrnran123_State*)_p_space);\n"
  "  }else{\n"
  "    _luniform_pick = 0.5;\n"
  "  }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE noiseFromRandom123() {\n"
  "VERBATIM\n"
  "#if !NRNBBCORE\n"
  " {\n"
  "  nrnran123_State** pv = (nrnran123_State**)(&_p_space);\n"
  "  if (*pv) {\n"
  "    nrnran123_deletestream(*pv);\n"
  "    *pv = (nrnran123_State*)0;\n"
  "  }\n"
  "  *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));\n"
  " }\n"
  "#endif\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "static void bbcore_write(double* x, int* d, int* xx, int *offset, _threadargsproto_) {\n"
  "#if !NRNBBCORE\n"
  "  /* error if using the legacy normrand */\n"
  "  if (!_p_space) {\n"
  "    fprintf(stderr, \"OStimHelper: noiseFromRandom123(1,2,3) not called.\\n\");\n"
  "    assert(0);\n"
  "  }\n"
  "  if (d) {\n"
  "    uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "    nrnran123_State** pv = (nrnran123_State**)(&_p_space);\n"
  "    nrnran123_getids3(*pv, di, di+1, di+2);\n"
  "  }\n"
  "#endif\n"
  "  *offset += 3;\n"
  "}\n"
  "\n"
  "static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {\n"
  "  uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "  nrnran123_State** pv = (nrnran123_State**)(&_p_space);\n"
  "#if !NRNBBCORE\n"
  "  if(*pv) {\n"
  "    nrnran123_deletestream(*pv);\n"
  "  }\n"
  "#endif\n"
  "  *pv = nrnran123_newstream3(di[0], di[1], di[2]);\n"
  "  *offset += 3;\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
