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
static constexpr auto number_of_datum_variables = 0;
static constexpr auto number_of_floating_point_variables = 1;
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
 
#define nrn_init _nrn_init__SynapseReader
#define _nrn_initial _nrn_initial__SynapseReader
#define nrn_cur _nrn_cur__SynapseReader
#define _nrn_current _nrn_current__SynapseReader
#define nrn_jacob _nrn_jacob__SynapseReader
#define nrn_state _nrn_state__SynapseReader
#define _net_receive _net_receive__SynapseReader 
 
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
#define v _ml->template fpfield<0>(_iml)
#define v_columnindex 0
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
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
 {"setdata_SynapseReader", _hoc_setdata},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {0, 0}
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"SynapseReader",
 0,
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 1);
 	/*initialize range parameters*/
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 1);
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _SynapseReader_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nullptr, nullptr, nullptr, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"v"} /* 0 */);
  hoc_register_prop_size(_mechtype, 1, 0);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 SynapseReader /mnt/mydata/cerebellum/mod/SynapseReader.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{

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
   _v = _vec_v[_ni[_iml]];
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
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/SynapseReader.mod";
    const char* nmodl_file_text = 
  "COMMENT\n"
  "/**\n"
  " * @file syn2reader.mod\n"
  " * @brief A Reader for Synapse files, inc Nrn and SYN2\n"
  " * @author F. Pereira\n"
  " * @date 2018-06\n"
  " * @remark Copyright (C) 2018 EPFL - Blue Brain Project\n"
  " * All rights reserved. Do not distribute without further notice.\n"
  " */\n"
  "\n"
  "// This file is a reminiscent of the Syntool reader mod which used to get synapses via\n"
  "// synapsetool, which supported a number of formats, inc NRN, SYN2 and SONATA edges.\n"
  "//\n"
  "// This section is particularly useful to given insight on the Neurodamus internal\n"
  "// synapse structure. Notice that we historically (from nrn) had 19 fields, out of which\n"
  "// we stored 11, later bumped to 14.\n"
  "\n"
  "// The NRN Synapse fields according to SYN2 transitional spec\n"
  "// ---------------------------------------------------------------------\n"
  "// | Synapse field name     [ND index] |  Gap-J fields name\n"
  "// ---------------------------------------------------------------------\n"
  "//   - connected_neurons_post          | connected_neurons_post (not loaded)\n"
  "//   0 connected_neurons_pre      [0]  | connected_neurons_pre\n"
  "//   1 delay                      [1]  | (N/A)\n"
  "//   2 morpho_section_id_post     [2]  | morpho_section_id_post\n"
  "//   3 morpho_segment_id_post     [3]  | morpho_section_id_post\n"
  "//   4 morpho_offset_segment_post [4]  | morpho_section_id_post\n"
  "//   5 morpho_section_id_pre           | morpho_section_id_pre (unused)\n"
  "//   6 morpho_segment_id_pre           | morpho_section_id_pre (unused)\n"
  "//   7 morpho_offset_segment_pre       | morpho_section_id_pre (unused)\n"
  "//   8 conductance                [5]  | conductance (required)\n"
  "//   9 u_syn                      [6]  | (N/A)\n"
  "//  10 depression_time            [7]  | junction_id_pre\n"
  "//  11 facilitation_time          [8]  | junction_id_post\n"
  "//  12 decay_time                 [9]  | N/A\n"
  "//  13 syn_type_id                [10] | N/A\n"
  "//  14 morpho_type_id_pre              | N/A\n"
  "//  15 morpho_branch_order_dend        | N/A\n"
  "//  16 morpho_branch_order_axon        | N/A\n"
  "//  17 n_rrp_vesicles (optional)  [11] | N/A\n"
  "//  18 morpho_section_type_pos         | N/A\n"
  "//  NA u_hill_coefficient         [12] | N/A\n"
  "//  NA conductance_scale_factor   [13] | N/A\n"
  "\n"
  "\n"
  "// The SYN2 v2 spec fields\n"
  "// : This spec deprecates compat with NRN\n"
  "// ---------------------------------------------------------------------\n"
  "// | Synapse field name     [ND index] |  Gap-J fields name\n"
  "// ---------------------------------------------------------------------\n"
  "//  connected_neurons_post            | connected_neurons_post (not loaded)\n"
  "//  connected_neurons_pre        [0]  | connected_neurons_pre\n"
  "//  delay                        [1]  | (N/A)\n"
  "//  morpho_section_id_post       [2]  | morpho_section_id_post\n"
  "//  (N/A)                        [3]  | (N/A)\n"
  "//  morpho_section_fraction_post [4]  | morpho_section_fraction_post\n"
  "//  conductance                  [5]  | conductance (required)\n"
  "//  u_syn                        [6]  | (N/A)\n"
  "//  depression_time              [7]  | junction_id_pre\n"
  "//  facilitation_time            [8]  | junction_id_post\n"
  "//  decay_time                   [9]  | (N/A)\n"
  "//  syn_type_id                  [10] | (N/A)\n"
  "//  n_rrp_vesicles  (required)   [11] | (N/A)\n"
  "//  u_hill_coefficient           [12] | (N/A)\n"
  "//  conductance_scale_factor     [13] | (N/A)\n"
  "\n"
  "\n"
  "// NEUROGLIAL Field Spec  (for reference only)\n"
  "// -----------------------------------\n"
  "// | Synapse field name     [ND index]\n"
  "// -----------------------------------\n"
  "//  source_node_id     [query field]\n"
  "//  target_node_id               [0]  !! NOTE: target\n"
  "//  synapse_id                   [1]\n"
  "//  morpho_section_id_pre        [2]\n"
  "//  morpho_segment_id_pre        [3]\n"
  "//  morpho_offset_segment_pre    [4]\n"
  "\n"
  "ENDCOMMENT\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
