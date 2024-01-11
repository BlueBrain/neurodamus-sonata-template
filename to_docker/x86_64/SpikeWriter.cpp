/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format on
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
// clang-format off
#include "neuron/cache/mechanism_range.hpp"
#include <vector>
using std::size_t;
static auto& std_cerr_stream = std::cerr;
static constexpr auto number_of_datum_variables = 2;
static constexpr auto number_of_floating_point_variables = 0;
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
 
#define nrn_init _nrn_init__SpikeWriter
#define _nrn_initial _nrn_initial__SpikeWriter
#define nrn_cur _nrn_cur__SpikeWriter
#define _nrn_current _nrn_current__SpikeWriter
#define nrn_jacob _nrn_jacob__SpikeWriter
#define nrn_state _nrn_state__SpikeWriter
#define _net_receive _net_receive__SpikeWriter 
#define write write__SpikeWriter 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _internalthreadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
#define _internalthreadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define _nd_area *_ml->dptr_field<0>(_iml)
 static _nrn_mechanism_cache_instance _ml_real{nullptr};
static _nrn_mechanism_cache_range *_ml{&_ml_real};
static size_t _iml{0};
static Datum *_ppvar;
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_write(void*);
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
 neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
_ppvar = _nrn_mechanism_access_dparam(_prop);
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
 {"write", _hoc_write},
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
 static double v = 0;
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
"SpikeWriter",
 0,
 0,
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
   _ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 0);
 	/*initialize range parameters*/
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 0);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _SpikeWriter_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nullptr, nullptr, nullptr, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */);
  hoc_register_prop_size(_mechtype, 0, 2);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
 add_nrn_artcell(_mechtype, 0);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 SpikeWriter /mnt/mydata/cerebellum/mod/SpikeWriter.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int write();
 
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#include <stdlib.h>

#ifndef DISABLE_MPI
#include <mpi.h>
#endif

#ifndef NRN_VERSION_GTEQ_8_2_0
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
#endif

extern int nrnmpi_myid;
#endif  // CORENEURON_BUILD
 
static int  write (  ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
        double *time = NULL, *gid = NULL;
        int num_spikes = 0;
        char *filePath = NULL;

        // first vector is time of spikes
        if (ifarg(1)) {
            IvocVect *v = vector_arg(1);
            time = vector_vec(v);
            num_spikes = vector_capacity(v);
        }

        // second vector is associated gids
        if (ifarg(2)) {
            IvocVect *v = vector_arg(2);
            gid = vector_vec(v);
        }

        // third argument is file path
        if(ifarg(3)) {
            filePath = hoc_gargstr(3);
        } else  {
#ifndef DISABLE_MPI
            if(nrnmpi_myid == 0) {
                fprintf(stderr, " Error : No spike file path provided, can't write spikes! \n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
#else
            fprintf(stderr, " Error : No spike file path provided, can't write spikes! \n");
            exit(-1);
#endif
        }

        // rank 0 write extra string at the begining as "/scatter"
        unsigned num_entries = nrnmpi_myid == 0 ? (num_spikes + 1) : num_spikes;

        // each spike record in the file is max 48 chars
        const int spike_record_length = 48;

        // amount of data for recording spikes +  zero termination
        unsigned num_bytes = (sizeof(char) * num_entries * spike_record_length) + 1;

        char *spike_data = (char *) malloc(num_bytes);

        if(spike_data == NULL) {
            fprintf(stderr, "Error : Memory allocation failed for spike buffer I/O!\n");
#ifndef DISABLE_MPI
            MPI_Abort(MPI_COMM_WORLD, 1);
#else
            exit(-1);
#endif
        }

        strcpy(spike_data, "");

        if(nrnmpi_myid == 0) {
            strcat(spike_data, "/scatter\n");
        }

        // populate buffer with all entries
        int i;
        for(i = 0; i < num_spikes; i++) {
            char str[spike_record_length];
            int nstr = snprintf(str, spike_record_length, "%.3f\t%d\n", time[i], (int)gid[i]);
            if (nstr >= spike_record_length) {
                fprintf(stderr, "Error : Record written is larger than spike record buffer\n");
                free(spike_data);
#ifndef DISABLE_MPI
                MPI_Abort(MPI_COMM_WORLD, 1);
#else
                exit(-1);
#endif
            }
            strcat(spike_data, str);
        }

        // calculate offset into global file. note that we don't write
        // num_bytes but only "populated" characters
        unsigned long num_chars = strlen(spike_data);

#ifndef DISABLE_MPI
        unsigned long offset = 0;

        MPI_Exscan(&num_chars, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (nrnmpi_myid == 0) {
            offset = 0;
        }

        // write to file using parallel mpi i/o Remove it first in case it exists.
        // Must delete because MPI_File_open does not have a Truncate mode
        MPI_File fh;
        MPI_Status status;
        MPI_File_delete(filePath, MPI_INFO_NULL);
        MPI_File_open(MPI_COMM_WORLD, filePath, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        MPI_File_write_at_all(fh, offset, spike_data, num_chars, MPI_BYTE, &status);

        MPI_File_close(&fh);
#else
        FILE *spike_file = fopen(filePath, "w");
        fprintf(spike_file, spike_data);
        fclose(spike_file);
#endif

        free(spike_data);
#endif  // CORENEURON_BUILD
  return 0; }
 
static double _hoc_write(void* _vptr) {
 double _r;
    auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _setdata(_p);
 _r = 1.;
 write (  );
 return(_r);
}

static void initmodel() {
  int _i; double _save;_ninits++;
{

}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
 v=_v;
{
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/mnt/mydata/cerebellum/mod/SpikeWriter.mod";
    const char* nmodl_file_text = 
  "COMMENT\n"
  "/**\n"
  " * @file SpikeWriter.mod\n"
  " * @brief Interface to write spikes in parallel using MPI-IO\n"
  " */\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    ARTIFICIAL_CELL SpikeWriter\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "#ifndef CORENEURON_BUILD\n"
  "#include <stdlib.h>\n"
  "\n"
  "#ifndef DISABLE_MPI\n"
  "#include <mpi.h>\n"
  "#endif\n"
  "\n"
  "#ifndef NRN_VERSION_GTEQ_8_2_0\n"
  "extern double* vector_vec();\n"
  "extern int vector_capacity();\n"
  "extern void* vector_arg();\n"
  "#endif\n"
  "\n"
  "extern int nrnmpi_myid;\n"
  "#endif  // CORENEURON_BUILD\n"
  "ENDVERBATIM\n"
  "\n"
  ": write_spikes to file in parallel\n"
  "PROCEDURE write() {\n"
  "\n"
  "    VERBATIM\n"
  "#ifndef CORENEURON_BUILD\n"
  "        double *time = NULL, *gid = NULL;\n"
  "        int num_spikes = 0;\n"
  "        char *filePath = NULL;\n"
  "\n"
  "        // first vector is time of spikes\n"
  "        if (ifarg(1)) {\n"
  "            IvocVect *v = vector_arg(1);\n"
  "            time = vector_vec(v);\n"
  "            num_spikes = vector_capacity(v);\n"
  "        }\n"
  "\n"
  "        // second vector is associated gids\n"
  "        if (ifarg(2)) {\n"
  "            IvocVect *v = vector_arg(2);\n"
  "            gid = vector_vec(v);\n"
  "        }\n"
  "\n"
  "        // third argument is file path\n"
  "        if(ifarg(3)) {\n"
  "            filePath = hoc_gargstr(3);\n"
  "        } else  {\n"
  "#ifndef DISABLE_MPI\n"
  "            if(nrnmpi_myid == 0) {\n"
  "                fprintf(stderr, \" Error : No spike file path provided, can't write spikes! \\n\");\n"
  "                MPI_Abort(MPI_COMM_WORLD, 1);\n"
  "            }\n"
  "#else\n"
  "            fprintf(stderr, \" Error : No spike file path provided, can't write spikes! \\n\");\n"
  "            exit(-1);\n"
  "#endif\n"
  "        }\n"
  "\n"
  "        // rank 0 write extra string at the begining as \"/scatter\"\n"
  "        unsigned num_entries = nrnmpi_myid == 0 ? (num_spikes + 1) : num_spikes;\n"
  "\n"
  "        // each spike record in the file is max 48 chars\n"
  "        const int spike_record_length = 48;\n"
  "\n"
  "        // amount of data for recording spikes +  zero termination\n"
  "        unsigned num_bytes = (sizeof(char) * num_entries * spike_record_length) + 1;\n"
  "\n"
  "        char *spike_data = (char *) malloc(num_bytes);\n"
  "\n"
  "        if(spike_data == NULL) {\n"
  "            fprintf(stderr, \"Error : Memory allocation failed for spike buffer I/O!\\n\");\n"
  "#ifndef DISABLE_MPI\n"
  "            MPI_Abort(MPI_COMM_WORLD, 1);\n"
  "#else\n"
  "            exit(-1);\n"
  "#endif\n"
  "        }\n"
  "\n"
  "        strcpy(spike_data, \"\");\n"
  "\n"
  "        if(nrnmpi_myid == 0) {\n"
  "            strcat(spike_data, \"/scatter\\n\");\n"
  "        }\n"
  "\n"
  "        // populate buffer with all entries\n"
  "        int i;\n"
  "        for(i = 0; i < num_spikes; i++) {\n"
  "            char str[spike_record_length];\n"
  "            int nstr = snprintf(str, spike_record_length, \"%.3f\\t%d\\n\", time[i], (int)gid[i]);\n"
  "            if (nstr >= spike_record_length) {\n"
  "                fprintf(stderr, \"Error : Record written is larger than spike record buffer\\n\");\n"
  "                free(spike_data);\n"
  "#ifndef DISABLE_MPI\n"
  "                MPI_Abort(MPI_COMM_WORLD, 1);\n"
  "#else\n"
  "                exit(-1);\n"
  "#endif\n"
  "            }\n"
  "            strcat(spike_data, str);\n"
  "        }\n"
  "\n"
  "        // calculate offset into global file. note that we don't write\n"
  "        // num_bytes but only \"populated\" characters\n"
  "        unsigned long num_chars = strlen(spike_data);\n"
  "\n"
  "#ifndef DISABLE_MPI\n"
  "        unsigned long offset = 0;\n"
  "\n"
  "        MPI_Exscan(&num_chars, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);\n"
  "\n"
  "        if (nrnmpi_myid == 0) {\n"
  "            offset = 0;\n"
  "        }\n"
  "\n"
  "        // write to file using parallel mpi i/o Remove it first in case it exists.\n"
  "        // Must delete because MPI_File_open does not have a Truncate mode\n"
  "        MPI_File fh;\n"
  "        MPI_Status status;\n"
  "        MPI_File_delete(filePath, MPI_INFO_NULL);\n"
  "        MPI_File_open(MPI_COMM_WORLD, filePath, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);\n"
  "        MPI_File_write_at_all(fh, offset, spike_data, num_chars, MPI_BYTE, &status);\n"
  "\n"
  "        MPI_File_close(&fh);\n"
  "#else\n"
  "        FILE *spike_file = fopen(filePath, \"w\");\n"
  "        fprintf(spike_file, spike_data);\n"
  "        fclose(spike_file);\n"
  "#endif\n"
  "\n"
  "        free(spike_data);\n"
  "#endif  // CORENEURON_BUILD\n"
  "    ENDVERBATIM\n"
  "}\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
