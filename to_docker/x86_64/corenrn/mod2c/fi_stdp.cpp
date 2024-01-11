/*********************************************************
Model Name      : FastInhibSTDP
Filename        : fi_stdp.mod
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
        "FastInhibSTDP",
        "tau1",
        "tau2",
        "gmax",
        "e",
        "wmax",
        "aLTP",
        "aLTD",
        "mgid",
        "ggid",
        "srcgid",
        0,
        "i",
        "interval",
        "tlast_pre",
        "tlast_post",
        "M",
        "P",
        "deltaw",
        "wsyn",
        0,
        "A",
        "B",
        0,
        0
    };


    /** all global variables */
    struct FastInhibSTDP_Store {
        int point_type{};
        double A0{};
        double B0{};
        int reset{};
        int mech_type{};
        double tauLTP{20};
        double tauLTD{20};
        double on{1};
        int slist1[2]{18, 19};
        int dlist1[2]{22, 23};
    };
    static_assert(std::is_trivially_copy_constructible_v<FastInhibSTDP_Store>);
    static_assert(std::is_trivially_move_constructible_v<FastInhibSTDP_Store>);
    static_assert(std::is_trivially_copy_assignable_v<FastInhibSTDP_Store>);
    static_assert(std::is_trivially_move_assignable_v<FastInhibSTDP_Store>);
    static_assert(std::is_trivially_destructible_v<FastInhibSTDP_Store>);
    FastInhibSTDP_Store FastInhibSTDP_global;


    /** all mechanism instance variables and global variables */
    struct FastInhibSTDP_Instance  {
        double* tau1{};
        const double* tau2{};
        const double* gmax{};
        const double* e{};
        const double* wmax{};
        const double* aLTP{};
        const double* aLTD{};
        const double* mgid{};
        const double* ggid{};
        const double* srcgid{};
        double* i{};
        double* interval{};
        double* tlast_pre{};
        double* tlast_post{};
        double* M{};
        double* P{};
        double* deltaw{};
        double* wsyn{};
        double* A{};
        double* B{};
        double* g{};
        double* factor{};
        double* DA{};
        double* DB{};
        double* v_unused{};
        double* g_unused{};
        double* tsave{};
        const double* node_area{};
        const int* point_process{};
        FastInhibSTDP_Store* global{&FastInhibSTDP_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"tauLTP_FastInhibSTDP", &FastInhibSTDP_global.tauLTP},
        {"tauLTD_FastInhibSTDP", &FastInhibSTDP_global.tauLTD},
        {"on_FastInhibSTDP", &FastInhibSTDP_global.on},
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
        return 27;
    }


    static inline int int_variables_size() {
        return 2;
    }


    static inline int get_mech_type() {
        return FastInhibSTDP_global.mech_type;
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
    static void nrn_private_constructor_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new FastInhibSTDP_Instance{};
        assert(inst->global == &FastInhibSTDP_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(FastInhibSTDP_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &FastInhibSTDP_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(FastInhibSTDP_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &FastInhibSTDP_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(FastInhibSTDP_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->tau1 = ml->data+0*pnodecount;
        inst->tau2 = ml->data+1*pnodecount;
        inst->gmax = ml->data+2*pnodecount;
        inst->e = ml->data+3*pnodecount;
        inst->wmax = ml->data+4*pnodecount;
        inst->aLTP = ml->data+5*pnodecount;
        inst->aLTD = ml->data+6*pnodecount;
        inst->mgid = ml->data+7*pnodecount;
        inst->ggid = ml->data+8*pnodecount;
        inst->srcgid = ml->data+9*pnodecount;
        inst->i = ml->data+10*pnodecount;
        inst->interval = ml->data+11*pnodecount;
        inst->tlast_pre = ml->data+12*pnodecount;
        inst->tlast_post = ml->data+13*pnodecount;
        inst->M = ml->data+14*pnodecount;
        inst->P = ml->data+15*pnodecount;
        inst->deltaw = ml->data+16*pnodecount;
        inst->wsyn = ml->data+17*pnodecount;
        inst->A = ml->data+18*pnodecount;
        inst->B = ml->data+19*pnodecount;
        inst->g = ml->data+20*pnodecount;
        inst->factor = ml->data+21*pnodecount;
        inst->DA = ml->data+22*pnodecount;
        inst->DB = ml->data+23*pnodecount;
        inst->v_unused = ml->data+24*pnodecount;
        inst->g_unused = ml->data+25*pnodecount;
        inst->tsave = ml->data+26*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = ml->pdata;
    }



    static void nrn_alloc_FastInhibSTDP(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);

        #endif
    }


    static inline void net_receive_kernel_FastInhibSTDP(double t, Point_process* pnt, FastInhibSTDP_Instance* inst, NrnThread* nt, Memb_list* ml, int weight_index, double flag) {
        int tid = pnt->_tid;
        int id = pnt->_i_instance;
        double v = 0;
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        double* data = ml->data;
        double* weights = nt->weights;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        double* w = weights + weight_index + 0;
        inst->tsave[id] = t;
        {
            if ((*w) >= 0.0) {
                inst->P[id] = inst->P[id] * exp((inst->tlast_pre[id] - t) / inst->global->tauLTP) + inst->aLTP[id];
                inst->interval[id] = inst->tlast_post[id] - t;
                inst->tlast_pre[id] = t;
                inst->deltaw[id] = inst->wmax[id] * inst->M[id] * exp(inst->interval[id] / inst->global->tauLTD);
            } else {
                inst->M[id] = inst->M[id] * exp((inst->tlast_post[id] - t) / inst->global->tauLTD) - inst->aLTD[id];
                inst->interval[id] = t - inst->tlast_pre[id];
                inst->tlast_post[id] = t;
                inst->deltaw[id] = inst->wmax[id] * inst->P[id] * exp( -inst->interval[id] / inst->global->tauLTP);
            }
            if (inst->global->on) {
                inst->wsyn[id] = inst->wsyn[id] + inst->deltaw[id];
                if (inst->wsyn[id] > inst->wmax[id]) {
                    inst->wsyn[id] = inst->wmax[id];
                }
                if (inst->wsyn[id] < 0.0) {
                    inst->wsyn[id] = 0.0;
                }
            }
            inst->A[id] = inst->A[id] + inst->wsyn[id] * inst->factor[id];
            inst->B[id] = inst->B[id] + inst->wsyn[id] * inst->factor[id];
        }
    }


    static void net_receive_FastInhibSTDP(Point_process* pnt, int weight_index, double flag) {
        NrnThread* nt = nrn_threads + pnt->_tid;
        Memb_list* ml = get_memb_list(nt);
        NetReceiveBuffer_t* nrb = ml->_net_receive_buffer;
        if (nrb->_cnt >= nrb->_size) {
            realloc_net_receive_buffer(nt, ml);
        }
        int id = nrb->_cnt;
        nrb->_pnt_index[id] = pnt-nt->pntprocs;
        nrb->_weight_index[id] = weight_index;
        nrb->_nrb_t[id] = nt->_t;
        nrb->_nrb_flag[id] = flag;
        nrb->_cnt++;
    }


    void net_buf_receive_FastInhibSTDP(NrnThread* nt) {
        Memb_list* ml = get_memb_list(nt);
        if (!ml) {
            return;
        }

        NetReceiveBuffer_t* nrb = ml->_net_receive_buffer;
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);
        int count = nrb->_displ_cnt;
        #pragma omp simd
        #pragma ivdep
        for (int i = 0; i < count; i++) {
            int start = nrb->_displ[i];
            int end = nrb->_displ[i+1];
            for (int j = start; j < end; j++) {
                int index = nrb->_nrb_index[j];
                int offset = nrb->_pnt_index[index];
                double t = nrb->_nrb_t[index];
                int weight_index = nrb->_weight_index[index];
                double flag = nrb->_nrb_flag[index];
                Point_process* point_process = nt->pntprocs + offset;
                net_receive_kernel_FastInhibSTDP(t, point_process, inst, nt, ml, weight_index, flag);
            }
        }
        nrb->_displ_cnt = 0;
        nrb->_cnt = 0;
    }


    /** initialize channel */
    void nrn_init_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                inst->tsave[id] = -1e20;
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->A[id] = inst->global->A0;
                inst->B[id] = inst->global->B0;
                double tp;
                if (inst->tau1[id] / inst->tau2[id] > .9999) {
                    inst->tau1[id] = .9999 * inst->tau2[id];
                }
                inst->A[id] = 0.0;
                inst->B[id] = 0.0;
                tp = (inst->tau1[id] * inst->tau2[id]) / (inst->tau2[id] - inst->tau1[id]) * log(inst->tau2[id] / inst->tau1[id]);
                inst->factor[id] =  -exp( -tp / inst->tau1[id]) + exp( -tp / inst->tau2[id]);
                inst->factor[id] = 1.0 / inst->factor[id];
                inst->interval[id] = 0.0;
                inst->tlast_pre[id] = 0.0;
                inst->tlast_post[id] = 0.0;
                inst->M[id] = 0.0;
                inst->P[id] = 0.0;
                inst->deltaw[id] = 0.0;
                inst->wsyn[id] = 0.0;
            }
        }
    }


    inline double nrn_current_FastInhibSTDP(int id, int pnodecount, FastInhibSTDP_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->g[id] = (inst->B[id] - inst->A[id]) * inst->gmax[id];
        inst->i[id] = inst->g[id] * (v - inst->e[id]);
        current += inst->i[id];
        return current;
    }


    /** update current */
    void nrn_cur_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        double* shadow_rhs = nt->_shadow_rhs;
        double* shadow_d = nt->_shadow_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            double g = nrn_current_FastInhibSTDP(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double rhs = nrn_current_FastInhibSTDP(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            double mfactor = 1.e2/inst->node_area[indexes[0*pnodecount + id]];
            g = g*mfactor;
            rhs = rhs*mfactor;
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            shadow_rhs[id] = rhs;
            shadow_d[id] = g;
        }
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            vec_rhs[node_id] -= shadow_rhs[id];
            vec_d[node_id] += shadow_d[id];
        }
    }


    /** update state */
    void nrn_state_FastInhibSTDP(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<FastInhibSTDP_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->A[id] = inst->A[id] + (1.0 - exp(nt->_dt * (( -1.0) / inst->tau1[id]))) * ( -(0.0) / (( -1.0) / inst->tau1[id]) - inst->A[id]);
            inst->B[id] = inst->B[id] + (1.0 - exp(nt->_dt * (( -1.0) / inst->tau2[id]))) * ( -(0.0) / (( -1.0) / inst->tau2[id]) - inst->B[id]);
        }
    }


    /** register channel with the simulator */
    void _fi_stdp_reg() {

        int mech_type = nrn_get_mechtype("FastInhibSTDP");
        FastInhibSTDP_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_FastInhibSTDP, nrn_cur_FastInhibSTDP, nullptr, nrn_state_FastInhibSTDP, nrn_init_FastInhibSTDP, nrn_private_constructor_FastInhibSTDP, nrn_private_destructor_FastInhibSTDP, first_pointer_var_index(), nullptr, nullptr, 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_net_receive_buffering(net_buf_receive_FastInhibSTDP, mech_type);
        set_pnt_receive(mech_type, net_receive_FastInhibSTDP, nullptr, num_net_receive_args());
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
