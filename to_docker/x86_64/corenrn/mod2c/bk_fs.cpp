/*********************************************************
Model Name      : bk_fs
Filename        : bk_fs.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Wed Sep 27 13:48:16 2023
Backend         : C (api-compatibility)
NMODL Compiler  : 0.6 [2ce4a2b9 2023-06-05 16:57:21 +0200]
*********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <coreneuron/nrnconf.h>
#include <coreneuron/sim/multicore.hpp>
#include <coreneuron/mechanism/register_mech.hpp>
#include <coreneuron/gpu/nrn_acc_manager.hpp>
#include <coreneuron/utils/randoms/nrnran123.h>
#include <coreneuron/nrniv/nrniv_decl.h>
#include <coreneuron/utils/ivocvect.hpp>
#include <coreneuron/utils/nrnoc_aux.hpp>
#include <coreneuron/mechanism/mech/mod2c_core_thread.hpp>
#include <coreneuron/sim/scopmath/newton_thread.hpp>


namespace coreneuron {


    /** constants used in nmodl from UNITS */
    static const double FARADAY = 0x1.81f0fae775425p+6;
    static const double R = 0x1.0a1013e8990bep+3;
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "bk_fs",
        "gbar_bk_fs",
        0,
        "ik_bk_fs",
        0,
        "o_bk_fs",
        0,
        0
    };


    /** all global variables */
    struct bk_fs_Store {
        int ca_type{};
        int k_type{};
        double o0{};
        int reset{};
        int mech_type{};
        double k1{0.003};
        double k4{0.009};
        double d1{0.84};
        double d4{1};
        double q{1};
        int slist1[1]{2};
        int dlist1[1]{7};
    };
    static_assert(std::is_trivially_copy_constructible_v<bk_fs_Store>);
    static_assert(std::is_trivially_move_constructible_v<bk_fs_Store>);
    static_assert(std::is_trivially_copy_assignable_v<bk_fs_Store>);
    static_assert(std::is_trivially_move_assignable_v<bk_fs_Store>);
    static_assert(std::is_trivially_destructible_v<bk_fs_Store>);
    bk_fs_Store bk_fs_global;


    /** all mechanism instance variables and global variables */
    struct bk_fs_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gbar{};
        double* ik{};
        double* o{};
        double* cai{};
        double* ek{};
        double* oinf{};
        double* otau{};
        double* Do{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_cai{};
        const double* ion_cao{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        bk_fs_Store* global{&bk_fs_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"k1_bk_fs", &bk_fs_global.k1},
        {"k4_bk_fs", &bk_fs_global.k4},
        {"d1_bk_fs", &bk_fs_global.d1},
        {"d4_bk_fs", &bk_fs_global.d4},
        {"q_bk_fs", &bk_fs_global.q},
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
        return 10;
    }


    static inline int int_variables_size() {
        return 5;
    }


    static inline int get_mech_type() {
        return bk_fs_global.mech_type;
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
    static void nrn_private_constructor_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new bk_fs_Instance{};
        assert(inst->global == &bk_fs_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(bk_fs_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &bk_fs_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(bk_fs_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &bk_fs_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(bk_fs_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gbar = ml->data+0*pnodecount;
        inst->ik = ml->data+1*pnodecount;
        inst->o = ml->data+2*pnodecount;
        inst->cai = ml->data+3*pnodecount;
        inst->ek = ml->data+4*pnodecount;
        inst->oinf = ml->data+5*pnodecount;
        inst->otau = ml->data+6*pnodecount;
        inst->Do = ml->data+7*pnodecount;
        inst->v_unused = ml->data+8*pnodecount;
        inst->g_unused = ml->data+9*pnodecount;
        inst->ion_cai = nt->_data;
        inst->ion_cao = nt->_data;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_bk_fs(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);

        #endif
    }


    inline int rate_bk_fs(int id, int pnodecount, bk_fs_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double ca);


    inline int rate_bk_fs(int id, int pnodecount, bk_fs_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double ca) {
        int ret_rate = 0;
        double a, b, sum, z;
        z = 1e-3 * 2.0 * FARADAY / (R * (*(inst->celsius) + 273.15));
        a = 0.48 * ca / (ca + inst->global->k1 * exp( -z * inst->global->d1 * arg_v));
        b = 0.28 / (1.0 + ca / (inst->global->k4 * exp( -z * inst->global->d4 * arg_v)));
        sum = a + b;
        inst->oinf[id] = a / sum;
        inst->otau[id] = 1.0 / sum;
        return ret_rate;
    }


    /** initialize channel */
    void nrn_init_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma ivdep
            #pragma omp simd
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->cai[id] = inst->ion_cai[indexes[0*pnodecount + id]];
                inst->ek[id] = inst->ion_ek[indexes[2*pnodecount + id]];
                inst->o[id] = inst->global->o0;
                {
                    double a, b, sum, z, v_in_1, ca_in_1;
                    v_in_1 = v;
                    ca_in_1 = inst->cai[id];
                    z = 1e-3 * 2.0 * FARADAY / (R * (*(inst->celsius) + 273.15));
                    a = 0.48 * ca_in_1 / (ca_in_1 + inst->global->k1 * exp( -z * inst->global->d1 * v_in_1));
                    b = 0.28 / (1.0 + ca_in_1 / (inst->global->k4 * exp( -z * inst->global->d4 * v_in_1)));
                    sum = a + b;
                    inst->oinf[id] = a / sum;
                    inst->otau[id] = 1.0 / sum;
                }
                inst->o[id] = inst->oinf[id];
            }
        }
    }


    inline double nrn_current_bk_fs(int id, int pnodecount, bk_fs_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->ik[id] = inst->gbar[id] * inst->o[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);

        #pragma ivdep
        #pragma omp simd
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->cai[id] = inst->ion_cai[indexes[0*pnodecount + id]];
            inst->ek[id] = inst->ion_ek[indexes[2*pnodecount + id]];
            double g = nrn_current_bk_fs(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_bk_fs(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            inst->ion_dikdv[indexes[4*pnodecount + id]] += (dik-inst->ik[id])/0.001;
            inst->ion_ik[indexes[3*pnodecount + id]] += inst->ik[id];
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            vec_rhs[node_id] -= rhs;
            vec_d[node_id] += g;
        }
    }


    /** update state */
    void nrn_state_bk_fs(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<bk_fs_Instance*>(ml->instance);

        #pragma ivdep
        #pragma omp simd
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->cai[id] = inst->ion_cai[indexes[0*pnodecount + id]];
            inst->ek[id] = inst->ion_ek[indexes[2*pnodecount + id]];
            {
                double a, b, sum, z, v_in_0, ca_in_0;
                v_in_0 = v;
                ca_in_0 = inst->cai[id];
                z = 1e-3 * 2.0 * FARADAY / (R * (*(inst->celsius) + 273.15));
                a = 0.48 * ca_in_0 / (ca_in_0 + inst->global->k1 * exp( -z * inst->global->d1 * v_in_0));
                b = 0.28 / (1.0 + ca_in_0 / (inst->global->k4 * exp( -z * inst->global->d4 * v_in_0)));
                sum = a + b;
                inst->oinf[id] = a / sum;
                inst->otau[id] = 1.0 / sum;
            }
            inst->o[id] = inst->o[id] + (1.0 - exp(nt->_dt * ((((( -1.0))) / inst->otau[id]) * (inst->global->q)))) * ( -((((inst->oinf[id])) / inst->otau[id]) * (inst->global->q)) / ((((( -1.0))) / inst->otau[id]) * (inst->global->q)) - inst->o[id]);
        }
    }


    /** register channel with the simulator */
    void _bk_fs_reg() {

        int mech_type = nrn_get_mechtype("bk_fs");
        bk_fs_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_bk_fs, nrn_cur_bk_fs, nullptr, nrn_state_bk_fs, nrn_init_bk_fs, nrn_private_constructor_bk_fs, nrn_private_destructor_bk_fs, first_pointer_var_index(), 1);
        bk_fs_global.ca_type = nrn_get_mechtype("ca_ion");
        bk_fs_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 1, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_dparam_semantics(mech_type, 3, "k_ion");
        hoc_register_dparam_semantics(mech_type, 4, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
