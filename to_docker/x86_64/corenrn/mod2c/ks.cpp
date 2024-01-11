/*********************************************************
Model Name      : ks
Filename        : ks.mod
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
        "ks",
        "gksbar_ks",
        0,
        "pinf_ks",
        "qinf_ks",
        "taup_ks",
        "tauq_ks",
        "gks_ks",
        0,
        "p_ks",
        "q_ks",
        0,
        0
    };


    /** all global variables */
    struct ks_Store {
        int k_type{};
        double p0{};
        double q0{};
        int reset{};
        int mech_type{};
        double q10{2};
        double vhalfp{-35};
        double tp{6};
        double kp{2};
        double vhalfq{-50};
        double kq{6.6};
        double to{10};
        double tk{6.85};
        double a0q{100};
        double tvh{-50};
        int slist1[2]{6, 7};
        int dlist1[2]{9, 10};
    };
    static_assert(std::is_trivially_copy_constructible_v<ks_Store>);
    static_assert(std::is_trivially_move_constructible_v<ks_Store>);
    static_assert(std::is_trivially_copy_assignable_v<ks_Store>);
    static_assert(std::is_trivially_move_assignable_v<ks_Store>);
    static_assert(std::is_trivially_destructible_v<ks_Store>);
    ks_Store ks_global;


    /** all mechanism instance variables and global variables */
    struct ks_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gksbar{};
        double* pinf{};
        double* qinf{};
        double* taup{};
        double* tauq{};
        double* gks{};
        double* p{};
        double* q{};
        double* ek{};
        double* Dp{};
        double* Dq{};
        double* ik{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        ks_Store* global{&ks_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"q10_ks", &ks_global.q10},
        {"vhalfp_ks", &ks_global.vhalfp},
        {"tp_ks", &ks_global.tp},
        {"kp_ks", &ks_global.kp},
        {"vhalfq_ks", &ks_global.vhalfq},
        {"kq_ks", &ks_global.kq},
        {"to_ks", &ks_global.to},
        {"tk_ks", &ks_global.tk},
        {"a0q_ks", &ks_global.a0q},
        {"tvh_ks", &ks_global.tvh},
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
        return 14;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return ks_global.mech_type;
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
    static void nrn_private_constructor_ks(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new ks_Instance{};
        assert(inst->global == &ks_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(ks_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_ks(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<ks_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &ks_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(ks_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<ks_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &ks_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(ks_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gksbar = ml->data+0*pnodecount;
        inst->pinf = ml->data+1*pnodecount;
        inst->qinf = ml->data+2*pnodecount;
        inst->taup = ml->data+3*pnodecount;
        inst->tauq = ml->data+4*pnodecount;
        inst->gks = ml->data+5*pnodecount;
        inst->p = ml->data+6*pnodecount;
        inst->q = ml->data+7*pnodecount;
        inst->ek = ml->data+8*pnodecount;
        inst->Dp = ml->data+9*pnodecount;
        inst->Dq = ml->data+10*pnodecount;
        inst->ik = ml->data+11*pnodecount;
        inst->v_unused = ml->data+12*pnodecount;
        inst->g_unused = ml->data+13*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_ks(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_ks(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<ks_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_ks(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<ks_Instance*>(ml->instance);

        #endif
    }


    inline int rates_ks(int id, int pnodecount, ks_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int rates_ks(int id, int pnodecount, ks_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_rates = 0;
        double qt;
        qt = pow(inst->global->q10, ((*(inst->celsius) - 25.0) / 10.0));
        inst->pinf[id] = (1.0 / (1.0 + exp( -(arg_v - inst->global->vhalfp) / inst->global->kp)));
        inst->qinf[id] = (1.0 / (1.0 + exp((arg_v - inst->global->vhalfq) / inst->global->kq)));
        inst->taup[id] = inst->global->tp / qt;
        inst->tauq[id] = inst->global->to + inst->global->a0q / (1.0 + exp( -(arg_v - inst->global->tvh) / inst->global->tk)) / qt;
        return ret_rates;
    }


    /** initialize channel */
    void nrn_init_ks(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<ks_Instance*>(ml->instance);

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
                inst->p[id] = inst->global->p0;
                inst->q[id] = inst->global->q0;
                {
                    double qt, v_in_0;
                    v_in_0 = v;
                    qt = pow(inst->global->q10, ((*(inst->celsius) - 25.0) / 10.0));
                    inst->pinf[id] = (1.0 / (1.0 + exp( -(v_in_0 - inst->global->vhalfp) / inst->global->kp)));
                    inst->qinf[id] = (1.0 / (1.0 + exp((v_in_0 - inst->global->vhalfq) / inst->global->kq)));
                    inst->taup[id] = inst->global->tp / qt;
                    inst->tauq[id] = inst->global->to + inst->global->a0q / (1.0 + exp( -(v_in_0 - inst->global->tvh) / inst->global->tk)) / qt;
                }
                inst->p[id] = inst->pinf[id];
                inst->q[id] = inst->qinf[id];
            }
        }
    }


    inline double nrn_current_ks(int id, int pnodecount, ks_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->gks[id] = inst->gksbar[id] * inst->p[id] * inst->q[id];
        inst->ik[id] = inst->gks[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_ks(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<ks_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_ks(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_ks(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_ks(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<ks_Instance*>(ml->instance);

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
                double qt, v_in_1;
                v_in_1 = v;
                qt = pow(inst->global->q10, ((*(inst->celsius) - 25.0) / 10.0));
                inst->pinf[id] = (1.0 / (1.0 + exp( -(v_in_1 - inst->global->vhalfp) / inst->global->kp)));
                inst->qinf[id] = (1.0 / (1.0 + exp((v_in_1 - inst->global->vhalfq) / inst->global->kq)));
                inst->taup[id] = inst->global->tp / qt;
                inst->tauq[id] = inst->global->to + inst->global->a0q / (1.0 + exp( -(v_in_1 - inst->global->tvh) / inst->global->tk)) / qt;
            }
            inst->p[id] = inst->p[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->taup[id]))) * ( -(((inst->pinf[id])) / inst->taup[id]) / (((( -1.0))) / inst->taup[id]) - inst->p[id]);
            inst->q[id] = inst->q[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->tauq[id]))) * ( -(((inst->qinf[id])) / inst->tauq[id]) / (((( -1.0))) / inst->tauq[id]) - inst->q[id]);
        }
    }


    /** register channel with the simulator */
    void _ks_reg() {

        int mech_type = nrn_get_mechtype("ks");
        ks_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_ks, nrn_cur_ks, nullptr, nrn_state_ks, nrn_init_ks, nrn_private_constructor_ks, nrn_private_destructor_ks, first_pointer_var_index(), 1);
        ks_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
