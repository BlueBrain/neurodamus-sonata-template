/*********************************************************
Model Name      : VecStim
Filename        : VecStim.mod
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
        "VecStim",
        "ping",
        0,
        "index",
        "etime",
        0,
        0,
        0
    };


    /** all global variables */
    struct VecStim_Store {
        int point_type{};
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<VecStim_Store>);
    static_assert(std::is_trivially_move_constructible_v<VecStim_Store>);
    static_assert(std::is_trivially_copy_assignable_v<VecStim_Store>);
    static_assert(std::is_trivially_move_assignable_v<VecStim_Store>);
    static_assert(std::is_trivially_destructible_v<VecStim_Store>);
    VecStim_Store VecStim_global;


    /** all mechanism instance variables and global variables */
    struct VecStim_Instance  {
        const double* ping{};
        double* index{};
        double* etime{};
        double* space{};
        double* v_unused{};
        double* tsave{};
        const double* node_area{};
        void** point_process{};
        void** tqitem{};
        VecStim_Store* global{&VecStim_global};
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
        return -1;
    }


    static inline int num_net_receive_args() {
        return 1;
    }


    static inline int float_variables_size() {
        return 6;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return VecStim_global.mech_type;
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
    static void nrn_private_constructor_VecStim(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new VecStim_Instance{};
        assert(inst->global == &VecStim_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(VecStim_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_VecStim(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<VecStim_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &VecStim_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(VecStim_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<VecStim_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &VecStim_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(VecStim_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->ping = ml->data+0*pnodecount;
        inst->index = ml->data+1*pnodecount;
        inst->etime = ml->data+2*pnodecount;
        inst->space = ml->data+3*pnodecount;
        inst->v_unused = ml->data+4*pnodecount;
        inst->tsave = ml->data+5*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = nt->_vdata;
        inst->tqitem = nt->_vdata;
    }



    static void nrn_alloc_VecStim(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_VecStim(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<VecStim_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_VecStim(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<VecStim_Instance*>(ml->instance);

        #endif
    }


    inline double element_VecStim(int id, int pnodecount, VecStim_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int restartEvent_VecStim(int id, int pnodecount, VecStim_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int play_VecStim(int id, int pnodecount, VecStim_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


#ifdef STIM_DEBUG
# define debug_printf(...) printf(__VA_ARGS__)
#else
# define debug_printf(...)
#endif
#if defined(NRN_VERSION_GTEQ)
#if NRN_VERSION_GTEQ(9,0,0)
#define NRN_VERSION_GTEQ_9_0_0
#endif
#endif


#ifndef NRN_VERSION_GTEQ_8_2_0
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
#endif


namespace coreneuron {


    inline int restartEvent_VecStim(int id, int pnodecount, VecStim_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_restartEvent = 0;
        inst->index[id] = 0.0;
        #ifndef CORENEURON_BUILD
            while (element_VecStim(id, pnodecount, inst, data, indexes, thread, nt, v) && inst->etime[id] < nt->_t) {}  
            if (inst->index[id] > 0) {
                debug_printf("[VecStim] restartEvent(): index=%d, etime=%g, t=%g\n", (int)inst->index[id] - 1, inst->etime[id], nt->_t);
                #if defined(NRN_VERSION_GTEQ_9_0_0)
                artcell_net_send(&nt->_vdata[indexes[2*pnodecount + id]], (double*)0, indexes[1].get<Point_process*>(), inst->etime[id], 1.0);
                #else
                artcell_net_send(&nt->_vdata[indexes[2*pnodecount + id]], (double*)0, (Point_process*)indexes[1]._pvoid, inst->etime[id], 1.0);
                #endif
            }
        #endif

        return ret_restartEvent;
    }


    inline int play_VecStim(int id, int pnodecount, VecStim_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_play = 0;
            #ifndef CORENEURON_BUILD
            void** vv;
            vv = (void**)(&inst->space[id]);
            *vv = NULL;
            if (ifarg(1)) {
                *vv = vector_arg(1);
            }
            inst->index[id] = -2;
            #endif

        return ret_play;
    }


    inline double element_VecStim(int id, int pnodecount, VecStim_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_element = 0.0;
            const int i = (int)inst->index[id];
            IvocVect* const vv = *((IvocVect**)(&inst->space[id]));
            int size; double* px;
            if (i < 0 || vv == NULL)
                return 0;
            size = vector_capacity(vv);
            px = vector_vec(vv);
            if (i < size) {
                inst->etime[id] = px[i];
                inst->index[id] += 1.;
                debug_printf("[VecStim] element(): index=%d, etime=%g, t=%g\n", (int)inst->index[id] - 1, inst->etime[id], nt->_t);
                return inst->index[id];
            }
            inst->index[id] = -1;
            return 0;

        return ret_element;
    }


    static inline void net_receive_VecStim(Point_process* pnt, int weight_index, double flag) {
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
        auto* const inst = static_cast<VecStim_Instance*>(ml->instance);

        double t = nt->_t;
        inst->tsave[id] = t;
        {
            if (flag == 1.0) {
                        debug_printf("[VecStim] net_event(): index=%d, etime=%g, t=%g\n", (int)inst->index[id] - 1, inst->etime[id], t);

                net_event(pnt, t);
                if (element_VecStim(id, pnodecount, inst, data, indexes, thread, nt, v) > 0.0) {
                    if (inst->etime[id] < t) {
                        printf("[VecStim] WARNING: spike time (%g ms) before current time (%g ms)\n", inst->etime[id], t);
                    } else {
                        artcell_net_send(&inst->tqitem[indexes[2*pnodecount + id]], weight_index, pnt, nt->_t+inst->etime[id] - t, 1.0);
                    }
                }
            } else if (flag == 2.0) {
                if (inst->index[id] ==  -2.0) {
                    printf("[VecStim] Detected new time vector.\n");
                    {
                        inst->index[id] = 0.0;
                        #ifndef CORENEURON_BUILD
                            while (element_VecStim(id, pnodecount, inst, data, indexes, thread, nt, v) && inst->etime[id] < t) {}  
                            if (inst->index[id] > 0) {
                                debug_printf("[VecStim] restartEvent(): index=%d, etime=%g, t=%g\n", (int)inst->index[id] - 1, inst->etime[id], t);
                                #if defined(NRN_VERSION_GTEQ_9_0_0)
                                artcell_net_send(&nt->_vdata[indexes[2*pnodecount + id]], (double*)0, indexes[1].get<Point_process*>(), inst->etime[id], 1.0);
                                #else
                                artcell_net_send(&nt->_vdata[indexes[2*pnodecount + id]], (double*)0, (Point_process*)indexes[1]._pvoid, inst->etime[id], 1.0);
                                #endif
                            }
                        #endif

                    }
                }
                artcell_net_send(&inst->tqitem[indexes[2*pnodecount + id]], weight_index, pnt, nt->_t+inst->ping[id], 2.0);
            }
        }
    }


    /** initialize channel */
    void nrn_init_VecStim(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<VecStim_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                inst->tsave[id] = -1e20;
                double v = 0.0;
                 #ifndef CORENEURON_BUILD

                inst->index[id] = 0.0;
                if (element_VecStim(id, pnodecount, inst, data, indexes, thread, nt, v)) {
                    artcell_net_send(&inst->tqitem[indexes[2*pnodecount + id]], 0, (Point_process*)inst->point_process[indexes[1*pnodecount + id]], nt->_t+inst->etime[id] - nt->_t, 1.0);
                }
                if (inst->ping[id] > 0.0) {
                    artcell_net_send(&inst->tqitem[indexes[2*pnodecount + id]], 0, (Point_process*)inst->point_process[indexes[1*pnodecount + id]], nt->_t+inst->ping[id], 2.0);
                }
                 #endif

            }
        }
    }


    /** register channel with the simulator */
    void _VecStim_reg() {

        int mech_type = nrn_get_mechtype("VecStim");
        VecStim_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_VecStim, nullptr, nullptr, nullptr, nrn_init_VecStim, nrn_private_constructor_VecStim, nrn_private_destructor_VecStim, first_pointer_var_index(), nullptr, nullptr, 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "netsend");
        add_nrn_has_net_event(mech_type);
        add_nrn_artcell(mech_type, 2);
        set_pnt_receive(mech_type, net_receive_VecStim, nullptr, num_net_receive_args());
        hoc_register_net_send_buffering(mech_type);
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
