/*********************************************************
Model Name      : kdrmt
Filename        : kdrmt.mod
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
        "kdrmt",
        "gbar_kdrmt",
        "vhalfm_kdrmt",
        "q10_kdrmt",
        "alpm_kdrmt",
        "betm_kdrmt",
        0,
        "minf_kdrmt",
        "mtau_kdrmt",
        0,
        "m_kdrmt",
        0,
        0
    };


    /** all global variables */
    struct kdrmt_Store {
        int k_type{};
        double m0{};
        int reset{};
        int mech_type{};
        double a0m{0.0035};
        double zetam{0.055};
        double gmm{0.5};
        int slist1[1]{7};
        int dlist1[1]{10};
    };
    static_assert(std::is_trivially_copy_constructible_v<kdrmt_Store>);
    static_assert(std::is_trivially_move_constructible_v<kdrmt_Store>);
    static_assert(std::is_trivially_copy_assignable_v<kdrmt_Store>);
    static_assert(std::is_trivially_move_assignable_v<kdrmt_Store>);
    static_assert(std::is_trivially_destructible_v<kdrmt_Store>);
    kdrmt_Store kdrmt_global;


    /** all mechanism instance variables and global variables */
    struct kdrmt_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gbar{};
        const double* vhalfm{};
        const double* q10{};
        double* alpm{};
        double* betm{};
        double* minf{};
        double* mtau{};
        double* m{};
        double* ek{};
        double* ik{};
        double* Dm{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        kdrmt_Store* global{&kdrmt_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"a0m_kdrmt", &kdrmt_global.a0m},
        {"zetam_kdrmt", &kdrmt_global.zetam},
        {"gmm_kdrmt", &kdrmt_global.gmm},
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
        return 13;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return kdrmt_global.mech_type;
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
    static void nrn_private_constructor_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new kdrmt_Instance{};
        assert(inst->global == &kdrmt_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(kdrmt_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &kdrmt_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(kdrmt_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &kdrmt_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(kdrmt_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gbar = ml->data+0*pnodecount;
        inst->vhalfm = ml->data+1*pnodecount;
        inst->q10 = ml->data+2*pnodecount;
        inst->alpm = ml->data+3*pnodecount;
        inst->betm = ml->data+4*pnodecount;
        inst->minf = ml->data+5*pnodecount;
        inst->mtau = ml->data+6*pnodecount;
        inst->m = ml->data+7*pnodecount;
        inst->ek = ml->data+8*pnodecount;
        inst->ik = ml->data+9*pnodecount;
        inst->Dm = ml->data+10*pnodecount;
        inst->v_unused = ml->data+11*pnodecount;
        inst->g_unused = ml->data+12*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_kdrmt(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);

        #endif
    }


    inline int trates_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline int falpm_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline int fbetm_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int trates_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_trates = 0;
        double qt;
        qt = pow(inst->q10[id], ((*(inst->celsius) - 24.0) / 10.0));
        inst->minf[id] = 1.0 / (1.0 + exp( -(arg_v - 21.0) / 10.0));
        {
            double v_in_0;
            v_in_0 = arg_v;
            inst->alpm[id] = exp(inst->global->zetam * (v_in_0 - inst->vhalfm[id]));
        }
        {
            double v_in_1;
            v_in_1 = arg_v;
            inst->betm[id] = exp(inst->global->zetam * inst->global->gmm * (v_in_1 - inst->vhalfm[id]));
        }
        inst->mtau[id] = inst->betm[id] / (qt * inst->global->a0m * (1.0 + inst->alpm[id]));
        return ret_trates;
    }


    inline int falpm_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_falpm = 0;
        inst->alpm[id] = exp(inst->global->zetam * (arg_v - inst->vhalfm[id]));
        return ret_falpm;
    }


    inline int fbetm_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_fbetm = 0;
        inst->betm[id] = exp(inst->global->zetam * inst->global->gmm * (arg_v - inst->vhalfm[id]));
        return ret_fbetm;
    }


    /** initialize channel */
    void nrn_init_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);

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
                {
                    double qt, v_in_2;
                    v_in_2 = v;
                    qt = pow(inst->q10[id], ((*(inst->celsius) - 24.0) / 10.0));
                    inst->minf[id] = 1.0 / (1.0 + exp( -(v_in_2 - 21.0) / 10.0));
                    {
                        double v_in_0;
                        v_in_0 = v_in_2;
                        inst->alpm[id] = exp(inst->global->zetam * (v_in_0 - inst->vhalfm[id]));
                    }
                    {
                        double v_in_1;
                        v_in_1 = v_in_2;
                        inst->betm[id] = exp(inst->global->zetam * inst->global->gmm * (v_in_1 - inst->vhalfm[id]));
                    }
                    inst->mtau[id] = inst->betm[id] / (qt * inst->global->a0m * (1.0 + inst->alpm[id]));
                }
                inst->m[id] = inst->minf[id];
            }
        }
    }


    inline double nrn_current_kdrmt(int id, int pnodecount, kdrmt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->ik[id] = inst->gbar[id] * inst->m[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_kdrmt(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_kdrmt(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_kdrmt(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kdrmt_Instance*>(ml->instance);

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
                double qt, v_in_3;
                v_in_3 = v;
                qt = pow(inst->q10[id], ((*(inst->celsius) - 24.0) / 10.0));
                inst->minf[id] = 1.0 / (1.0 + exp( -(v_in_3 - 21.0) / 10.0));
                {
                    double v_in_0;
                    v_in_0 = v_in_3;
                    inst->alpm[id] = exp(inst->global->zetam * (v_in_0 - inst->vhalfm[id]));
                }
                {
                    double v_in_1;
                    v_in_1 = v_in_3;
                    inst->betm[id] = exp(inst->global->zetam * inst->global->gmm * (v_in_1 - inst->vhalfm[id]));
                }
                inst->mtau[id] = inst->betm[id] / (qt * inst->global->a0m * (1.0 + inst->alpm[id]));
            }
            inst->m[id] = inst->m[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->mtau[id]))) * ( -(((inst->minf[id])) / inst->mtau[id]) / (((( -1.0))) / inst->mtau[id]) - inst->m[id]);
        }
    }


    /** register channel with the simulator */
    void _kdrmt_reg() {

        int mech_type = nrn_get_mechtype("kdrmt");
        kdrmt_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_kdrmt, nrn_cur_kdrmt, nullptr, nrn_state_kdrmt, nrn_init_kdrmt, nrn_private_constructor_kdrmt, nrn_private_destructor_kdrmt, first_pointer_var_index(), 1);
        kdrmt_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
