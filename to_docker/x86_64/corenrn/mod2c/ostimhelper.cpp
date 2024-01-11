/*********************************************************
Model Name      : OdorStimHelper
Filename        : ostimhelper.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Thu Jan 11 15:15:20 2024
Backend         : C++ (api-compatibility)
NMODL Compiler  : 0.6 [6f6db3b6 2023-10-20 15:10:33 +0200]
*********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <coreneuron/gpu/nrn_acc_manager.hpp>
#include <coreneuron/mechanism/mech/mod2c_core_thread.hpp>
#include <coreneuron/mechanism/register_mech.hpp>
#include <coreneuron/nrnconf.h>
#include <coreneuron/nrniv/nrniv_decl.h>
#include <coreneuron/sim/multicore.hpp>
#include <coreneuron/sim/scopmath/newton_thread.hpp>
#include <coreneuron/utils/ivocvect.hpp>
#include <coreneuron/utils/nrnoc_aux.hpp>
#include <coreneuron/utils/randoms/nrnran123.h>


namespace coreneuron {
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "OdorStimHelper",
        "start",
        "dur",
        "invl_min",
        "invl_max",
        0,
        0,
        0,
        "space",
        0
    };


    /** all global variables */
    struct OdorStimHelper_Store {
        int point_type{};
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<OdorStimHelper_Store>);
    static_assert(std::is_trivially_move_constructible_v<OdorStimHelper_Store>);
    static_assert(std::is_trivially_copy_assignable_v<OdorStimHelper_Store>);
    static_assert(std::is_trivially_move_assignable_v<OdorStimHelper_Store>);
    static_assert(std::is_trivially_destructible_v<OdorStimHelper_Store>);
    OdorStimHelper_Store OdorStimHelper_global;


    /** all mechanism instance variables and global variables */
    struct OdorStimHelper_Instance  {
        const double* start{};
        const double* dur{};
        const double* invl_min{};
        const double* invl_max{};
        double* v_unused{};
        double* tsave{};
        const double* node_area{};
        void** point_process{};
        void** space{};
        void** tqitem{};
        OdorStimHelper_Store* global{&OdorStimHelper_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {nullptr, nullptr}
    };


    /** connect global (array) variables to hoc -- */
    static DoubVec hoc_vector_double[] = {
        {nullptr, nullptr, 0}
    };


    static inline int first_pointer_var_index() {
        return 2;
    }


    static inline int num_net_receive_args() {
        return 1;
    }


    static inline int float_variables_size() {
        return 6;
    }


    static inline int int_variables_size() {
        return 4;
    }


    static inline int get_mech_type() {
        return OdorStimHelper_global.mech_type;
    }


    static inline Memb_list* get_memb_list(NrnThread* nt) {
        if (!nt->_ml_list) {
            return nullptr;
        }
        return nt->_ml_list[get_mech_type()];
    }


    static inline void* mem_alloc(size_t num, size_t size, size_t alignment = 16) {
        void* ptr;
        posix_memalign(&ptr, alignment, num*size);
        memset(ptr, 0, size);
        return ptr;
    }


    static inline void mem_free(void* ptr) {
        free(ptr);
    }


    static inline void coreneuron_abort() {
        abort();
    }

    // Allocate instance structure
    static void nrn_private_constructor_OdorStimHelper(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new OdorStimHelper_Instance{};
        assert(inst->global == &OdorStimHelper_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(OdorStimHelper_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_OdorStimHelper(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<OdorStimHelper_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &OdorStimHelper_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(OdorStimHelper_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<OdorStimHelper_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &OdorStimHelper_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(OdorStimHelper_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->start = ml->data+0*pnodecount;
        inst->dur = ml->data+1*pnodecount;
        inst->invl_min = ml->data+2*pnodecount;
        inst->invl_max = ml->data+3*pnodecount;
        inst->v_unused = ml->data+4*pnodecount;
        inst->tsave = ml->data+5*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = nt->_vdata;
        inst->space = nt->_vdata;
        inst->tqitem = nt->_vdata;
    }



    static void nrn_alloc_OdorStimHelper(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_OdorStimHelper(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<OdorStimHelper_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_OdorStimHelper(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<OdorStimHelper_Instance*>(ml->instance);

        #endif
    }


    inline double invl_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double uniform_pick_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int initstream_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int noiseFromRandom123_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


static void bbcore_write(double* x, int* d, int* xx, int *offset, int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
#if !NRNBBCORE
  /* error if using the legacy normrand */
  if (!nt->_vdata[indexes[2*pnodecount + id]]) {
    fprintf(stderr, "OStimHelper: noiseFromRandom123(1,2,3) not called.\n");
    assert(0);
  }
  if (d) {
    uint32_t* di = ((uint32_t*)d) + *offset;
    nrnran123_State** pv = (nrnran123_State**)(&nt->_vdata[indexes[2*pnodecount + id]]);
    nrnran123_getids3(*pv, di, di+1, di+2);
  }
#endif
  *offset += 3;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
  uint32_t* di = ((uint32_t*)d) + *offset;
  nrnran123_State** pv = (nrnran123_State**)(&nt->_vdata[indexes[2*pnodecount + id]]);
#if !NRNBBCORE
  if(*pv) {
    nrnran123_deletestream(*pv);
  }
#endif
  *pv = nrnran123_newstream3(di[0], di[1], di[2]);
  *offset += 3;
}


namespace coreneuron {


    inline int initstream_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_initstream = 0;
          if (inst->space[indexes[2*pnodecount + id]]) {
            nrnran123_setseq((nrnran123_State*)inst->space[indexes[2*pnodecount + id]], 0, 0);
          }

        return ret_initstream;
    }


    inline int noiseFromRandom123_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_noiseFromRandom123 = 0;
        #if !NRNBBCORE
         {
          nrnran123_State** pv = (nrnran123_State**)(&inst->space[indexes[2*pnodecount + id]]);
          if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
          }
          *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));
         }
        #endif

        return ret_noiseFromRandom123;
    }


    inline double invl_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_invl = 0.0;
        double uniform_pick_in_0;
        {
              if (inst->space[indexes[2*pnodecount + id]]) {
                uniform_pick_in_0 = nrnran123_dblpick((nrnran123_State*)inst->space[indexes[2*pnodecount + id]]);
              }else{
                uniform_pick_in_0 = 0.5;
              }

        }
        ret_invl = (inst->invl_max[id] - inst->invl_min[id]) * uniform_pick_in_0 + inst->invl_min[id];
        return ret_invl;
    }


    inline double uniform_pick_OdorStimHelper(int id, int pnodecount, OdorStimHelper_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_uniform_pick = 0.0;
          if (inst->space[indexes[2*pnodecount + id]]) {
            ret_uniform_pick = nrnran123_dblpick((nrnran123_State*)inst->space[indexes[2*pnodecount + id]]);
          }else{
            ret_uniform_pick = 0.5;
          }

        return ret_uniform_pick;
    }


    static inline void net_receive_OdorStimHelper(Point_process* pnt, int weight_index, double flag) {
        int tid = pnt->_tid;
        int id = pnt->_i_instance;
        double v = 0;
        NrnThread* nt = nrn_threads + tid;
        Memb_list* ml = nt->_ml_list[pnt->_type];
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        double* data = ml->data;
        double* weights = nt->weights;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<OdorStimHelper_Instance*>(ml->instance);

        double t = nt->_t;
        inst->tsave[id] = t;
        {
            double tinv;
            if (flag == 1.0) {
                double invl_in_0;
                net_event(pnt, t);
                {
                    double uniform_pick_in_0;
                    {
                          if (inst->space[indexes[2*pnodecount + id]]) {
                            uniform_pick_in_0 = nrnran123_dblpick((nrnran123_State*)inst->space[indexes[2*pnodecount + id]]);
                          }else{
                            uniform_pick_in_0 = 0.5;
                          }

                    }
                    invl_in_0 = (inst->invl_max[id] - inst->invl_min[id]) * uniform_pick_in_0 + inst->invl_min[id];
                }
                tinv = invl_in_0;
                if (t + tinv <= inst->dur[id]) {
                    artcell_net_send(&inst->tqitem[indexes[3*pnodecount + id]], weight_index, pnt, nt->_t+tinv, 1.0);
                }
            }
        }
    }


    /** initialize channel */
    void nrn_init_OdorStimHelper(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<OdorStimHelper_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                inst->tsave[id] = -1e20;
                double v = 0.0;
                {
                      if (inst->space[indexes[2*pnodecount + id]]) {
                        nrnran123_setseq((nrnran123_State*)inst->space[indexes[2*pnodecount + id]], 0, 0);
                      }

                }
                if (nt->_t <= inst->start[id]) {
                    artcell_net_send(&inst->tqitem[indexes[3*pnodecount + id]], 0, (Point_process*)inst->point_process[indexes[1*pnodecount + id]], nt->_t+inst->start[id], 1.0);
                }
            }
        }
    }


    /** register channel with the simulator */
    void _ostimhelper_reg() {

        int mech_type = nrn_get_mechtype("OdorStimHelper");
        OdorStimHelper_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_OdorStimHelper, nullptr, nullptr, nullptr, nrn_init_OdorStimHelper, nrn_private_constructor_OdorStimHelper, nrn_private_destructor_OdorStimHelper, first_pointer_var_index(), nullptr, nullptr, 1);

        hoc_reg_bbcore_read(mech_type, bbcore_read);
        hoc_reg_bbcore_write(mech_type, bbcore_write);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "bbcorepointer");
        hoc_register_dparam_semantics(mech_type, 3, "netsend");
        add_nrn_has_net_event(mech_type);
        add_nrn_artcell(mech_type, 3);
        set_pnt_receive(mech_type, net_receive_OdorStimHelper, nullptr, num_net_receive_args());
        hoc_register_net_send_buffering(mech_type);
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
