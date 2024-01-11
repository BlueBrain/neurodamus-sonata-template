/*********************************************************
Model Name      : Kv2_2_0010
Filename        : Kv22.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Thu Jan 11 15:15:19 2024
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
        "Kv2_2_0010",
        "gKv2_2bar_Kv2_2_0010",
        "BBiD_Kv2_2_0010",
        0,
        "ik_Kv2_2_0010",
        "gKv2_2_Kv2_2_0010",
        0,
        "m_Kv2_2_0010",
        "h_Kv2_2_0010",
        0,
        0
    };


    /** all global variables */
    struct Kv2_2_0010_Store {
        int k_type{};
        double m0{};
        double h0{};
        int reset{};
        int mech_type{};
        int slist1[2]{4, 5};
        int dlist1[2]{11, 12};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kv2_2_0010_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kv2_2_0010_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kv2_2_0010_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kv2_2_0010_Store>);
    static_assert(std::is_trivially_destructible_v<Kv2_2_0010_Store>);
    Kv2_2_0010_Store Kv2_2_0010_global;


    /** all mechanism instance variables and global variables */
    struct Kv2_2_0010_Instance  {
        const double* gKv2_2bar{};
        const double* BBiD{};
        double* ik{};
        double* gKv2_2{};
        double* m{};
        double* h{};
        double* ek{};
        double* mInf{};
        double* mTau{};
        double* hInf{};
        double* hTau{};
        double* Dm{};
        double* Dh{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Kv2_2_0010_Store* global{&Kv2_2_0010_global};
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


    static inline int float_variables_size() {
        return 15;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return Kv2_2_0010_global.mech_type;
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
    static void nrn_private_constructor_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kv2_2_0010_Instance{};
        assert(inst->global == &Kv2_2_0010_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kv2_2_0010_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv2_2_0010_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv2_2_0010_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv2_2_0010_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv2_2_0010_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gKv2_2bar = ml->data+0*pnodecount;
        inst->BBiD = ml->data+1*pnodecount;
        inst->ik = ml->data+2*pnodecount;
        inst->gKv2_2 = ml->data+3*pnodecount;
        inst->m = ml->data+4*pnodecount;
        inst->h = ml->data+5*pnodecount;
        inst->ek = ml->data+6*pnodecount;
        inst->mInf = ml->data+7*pnodecount;
        inst->mTau = ml->data+8*pnodecount;
        inst->hInf = ml->data+9*pnodecount;
        inst->hTau = ml->data+10*pnodecount;
        inst->Dm = ml->data+11*pnodecount;
        inst->Dh = ml->data+12*pnodecount;
        inst->v_unused = ml->data+13*pnodecount;
        inst->g_unused = ml->data+14*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Kv2_2_0010(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);

        #endif
    }


    inline int rates_Kv2_2_0010(int id, int pnodecount, Kv2_2_0010_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);


    inline int rates_Kv2_2_0010(int id, int pnodecount, Kv2_2_0010_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_rates = 0;
        inst->mInf[id] = 1.0 / (1.0 + exp(((v - (5.0)) / ( -12.0))));
        inst->mTau[id] = 130.0 / (1.0 + exp(((v - ( -46.560)) / ( -44.140))));
        inst->hInf[id] = 1.0 / (1.0 + exp(((v - ( -16.300)) / (4.800))));
        inst->hTau[id] = 10000.0 / (1.0 + exp(((v - ( -46.560)) / ( -44.140))));
        return ret_rates;
    }


    /** initialize channel */
    void nrn_init_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
                inst->m[id] = inst->global->m0;
                inst->h[id] = inst->global->h0;
                {
                    inst->mInf[id] = 1.0 / (1.0 + exp(((v - (5.0)) / ( -12.0))));
                    inst->mTau[id] = 130.0 / (1.0 + exp(((v - ( -46.560)) / ( -44.140))));
                    inst->hInf[id] = 1.0 / (1.0 + exp(((v - ( -16.300)) / (4.800))));
                    inst->hTau[id] = 10000.0 / (1.0 + exp(((v - ( -46.560)) / ( -44.140))));
                }
                inst->m[id] = inst->mInf[id];
                inst->h[id] = inst->hInf[id];
            }
        }
    }


    inline double nrn_current_Kv2_2_0010(int id, int pnodecount, Kv2_2_0010_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->gKv2_2[id] = inst->gKv2_2bar[id] * inst->m[id] * inst->h[id];
        inst->ik[id] = inst->gKv2_2[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_Kv2_2_0010(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kv2_2_0010(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            inst->ion_dikdv[indexes[2*pnodecount + id]] += (dik-inst->ik[id])/0.001;
            inst->ion_ik[indexes[1*pnodecount + id]] += inst->ik[id];
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            vec_rhs[node_id] -= rhs;
            vec_d[node_id] += g;
        }
    }


    /** update state */
    void nrn_state_Kv2_2_0010(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv2_2_0010_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            {
                inst->mInf[id] = 1.0 / (1.0 + exp(((v - (5.0)) / ( -12.0))));
                inst->mTau[id] = 130.0 / (1.0 + exp(((v - ( -46.560)) / ( -44.140))));
                inst->hInf[id] = 1.0 / (1.0 + exp(((v - ( -16.300)) / (4.800))));
                inst->hTau[id] = 10000.0 / (1.0 + exp(((v - ( -46.560)) / ( -44.140))));
            }
            inst->m[id] = inst->m[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->mTau[id]))) * ( -(((inst->mInf[id])) / inst->mTau[id]) / (((( -1.0))) / inst->mTau[id]) - inst->m[id]);
            inst->h[id] = inst->h[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->hTau[id]))) * ( -(((inst->hInf[id])) / inst->hTau[id]) / (((( -1.0))) / inst->hTau[id]) - inst->h[id]);
        }
    }


    /** register channel with the simulator */
    void _Kv22_reg() {

        int mech_type = nrn_get_mechtype("Kv2_2_0010");
        Kv2_2_0010_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kv2_2_0010, nrn_cur_Kv2_2_0010, nullptr, nrn_state_Kv2_2_0010, nrn_init_Kv2_2_0010, nrn_private_constructor_Kv2_2_0010, nrn_private_destructor_Kv2_2_0010, first_pointer_var_index(), 1);
        Kv2_2_0010_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
