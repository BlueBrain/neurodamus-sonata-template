/*********************************************************
Model Name      : kamt
Filename        : kamt.mod
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
        "kamt",
        "gbar_kamt",
        "q10_kamt",
        0,
        "minf_kamt",
        "mtau_kamt",
        "hinf_kamt",
        "htau_kamt",
        0,
        "m_kamt",
        "h_kamt",
        0,
        0
    };


    /** all global variables */
    struct kamt_Store {
        int k_type{};
        double m0{};
        double h0{};
        int reset{};
        int mech_type{};
        double a0m{0.04};
        double vhalfm{-45};
        double zetam{0.1};
        double gmm{0.75};
        double a0h{0.018};
        double vhalfh{-70};
        double zetah{0.2};
        double gmh{0.99};
        double sha{9.9};
        double shi{5.7};
        int slist1[2]{6, 7};
        int dlist1[2]{10, 11};
    };
    static_assert(std::is_trivially_copy_constructible_v<kamt_Store>);
    static_assert(std::is_trivially_move_constructible_v<kamt_Store>);
    static_assert(std::is_trivially_copy_assignable_v<kamt_Store>);
    static_assert(std::is_trivially_move_assignable_v<kamt_Store>);
    static_assert(std::is_trivially_destructible_v<kamt_Store>);
    kamt_Store kamt_global;


    /** all mechanism instance variables and global variables */
    struct kamt_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gbar{};
        const double* q10{};
        double* minf{};
        double* mtau{};
        double* hinf{};
        double* htau{};
        double* m{};
        double* h{};
        double* ek{};
        double* ik{};
        double* Dm{};
        double* Dh{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        kamt_Store* global{&kamt_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"a0m_kamt", &kamt_global.a0m},
        {"vhalfm_kamt", &kamt_global.vhalfm},
        {"zetam_kamt", &kamt_global.zetam},
        {"gmm_kamt", &kamt_global.gmm},
        {"a0h_kamt", &kamt_global.a0h},
        {"vhalfh_kamt", &kamt_global.vhalfh},
        {"zetah_kamt", &kamt_global.zetah},
        {"gmh_kamt", &kamt_global.gmh},
        {"sha_kamt", &kamt_global.sha},
        {"shi_kamt", &kamt_global.shi},
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
        return kamt_global.mech_type;
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
    static void nrn_private_constructor_kamt(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new kamt_Instance{};
        assert(inst->global == &kamt_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(kamt_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_kamt(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &kamt_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(kamt_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &kamt_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(kamt_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gbar = ml->data+0*pnodecount;
        inst->q10 = ml->data+1*pnodecount;
        inst->minf = ml->data+2*pnodecount;
        inst->mtau = ml->data+3*pnodecount;
        inst->hinf = ml->data+4*pnodecount;
        inst->htau = ml->data+5*pnodecount;
        inst->m = ml->data+6*pnodecount;
        inst->h = ml->data+7*pnodecount;
        inst->ek = ml->data+8*pnodecount;
        inst->ik = ml->data+9*pnodecount;
        inst->Dm = ml->data+10*pnodecount;
        inst->Dh = ml->data+11*pnodecount;
        inst->v_unused = ml->data+12*pnodecount;
        inst->g_unused = ml->data+13*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_kamt(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_kamt(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_kamt(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);

        #endif
    }


    inline double alpm_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double betm_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double alph_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double beth_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline int trates_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int trates_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_trates = 0;
        double qt, betm_in_0, alpm_in_0, beth_in_0, alph_in_0;
        qt = pow(inst->q10[id], ((*(inst->celsius) - 24.0) / 10.0));
        inst->minf[id] = 1.0 / (1.0 + exp( -(arg_v - inst->global->sha - 7.6) / 14.0));
        {
            double v_in_0;
            v_in_0 = arg_v;
            betm_in_0 = exp(inst->global->zetam * inst->global->gmm * (v_in_0 - inst->global->vhalfm));
        }
        {
            double v_in_1;
            v_in_1 = arg_v;
            alpm_in_0 = exp(inst->global->zetam * (v_in_1 - inst->global->vhalfm));
        }
        inst->mtau[id] = betm_in_0 / (qt * inst->global->a0m * (1.0 + alpm_in_0));
        inst->hinf[id] = 1.0 / (1.0 + exp((arg_v - inst->global->shi + 47.4) / 6.0));
        {
            double v_in_2;
            v_in_2 = arg_v;
            beth_in_0 = exp(inst->global->zetah * inst->global->gmh * (v_in_2 - inst->global->vhalfh));
        }
        {
            double v_in_3;
            v_in_3 = arg_v;
            alph_in_0 = exp(inst->global->zetah * (v_in_3 - inst->global->vhalfh));
        }
        inst->htau[id] = beth_in_0 / (qt * inst->global->a0h * (1.0 + alph_in_0));
        return ret_trates;
    }


    inline double alpm_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alpm = 0.0;
        ret_alpm = exp(inst->global->zetam * (arg_v - inst->global->vhalfm));
        return ret_alpm;
    }


    inline double betm_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_betm = 0.0;
        ret_betm = exp(inst->global->zetam * inst->global->gmm * (arg_v - inst->global->vhalfm));
        return ret_betm;
    }


    inline double alph_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alph = 0.0;
        ret_alph = exp(inst->global->zetah * (arg_v - inst->global->vhalfh));
        return ret_alph;
    }


    inline double beth_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_beth = 0.0;
        ret_beth = exp(inst->global->zetah * inst->global->gmh * (arg_v - inst->global->vhalfh));
        return ret_beth;
    }


    /** initialize channel */
    void nrn_init_kamt(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);

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
                    double qt, betm_in_0, alpm_in_0, beth_in_0, alph_in_0, v_in_4;
                    v_in_4 = v;
                    qt = pow(inst->q10[id], ((*(inst->celsius) - 24.0) / 10.0));
                    inst->minf[id] = 1.0 / (1.0 + exp( -(v_in_4 - inst->global->sha - 7.6) / 14.0));
                    {
                        double v_in_0;
                        v_in_0 = v_in_4;
                        betm_in_0 = exp(inst->global->zetam * inst->global->gmm * (v_in_0 - inst->global->vhalfm));
                    }
                    {
                        double v_in_1;
                        v_in_1 = v_in_4;
                        alpm_in_0 = exp(inst->global->zetam * (v_in_1 - inst->global->vhalfm));
                    }
                    inst->mtau[id] = betm_in_0 / (qt * inst->global->a0m * (1.0 + alpm_in_0));
                    inst->hinf[id] = 1.0 / (1.0 + exp((v_in_4 - inst->global->shi + 47.4) / 6.0));
                    {
                        double v_in_2;
                        v_in_2 = v_in_4;
                        beth_in_0 = exp(inst->global->zetah * inst->global->gmh * (v_in_2 - inst->global->vhalfh));
                    }
                    {
                        double v_in_3;
                        v_in_3 = v_in_4;
                        alph_in_0 = exp(inst->global->zetah * (v_in_3 - inst->global->vhalfh));
                    }
                    inst->htau[id] = beth_in_0 / (qt * inst->global->a0h * (1.0 + alph_in_0));
                }
                inst->m[id] = inst->minf[id];
                inst->h[id] = inst->hinf[id];
            }
        }
    }


    inline double nrn_current_kamt(int id, int pnodecount, kamt_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->ik[id] = inst->gbar[id] * inst->m[id] * inst->h[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_kamt(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_kamt(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_kamt(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_kamt(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<kamt_Instance*>(ml->instance);

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
                double qt, betm_in_0, alpm_in_0, beth_in_0, alph_in_0, v_in_5;
                v_in_5 = v;
                qt = pow(inst->q10[id], ((*(inst->celsius) - 24.0) / 10.0));
                inst->minf[id] = 1.0 / (1.0 + exp( -(v_in_5 - inst->global->sha - 7.6) / 14.0));
                {
                    double v_in_0;
                    v_in_0 = v_in_5;
                    betm_in_0 = exp(inst->global->zetam * inst->global->gmm * (v_in_0 - inst->global->vhalfm));
                }
                {
                    double v_in_1;
                    v_in_1 = v_in_5;
                    alpm_in_0 = exp(inst->global->zetam * (v_in_1 - inst->global->vhalfm));
                }
                inst->mtau[id] = betm_in_0 / (qt * inst->global->a0m * (1.0 + alpm_in_0));
                inst->hinf[id] = 1.0 / (1.0 + exp((v_in_5 - inst->global->shi + 47.4) / 6.0));
                {
                    double v_in_2;
                    v_in_2 = v_in_5;
                    beth_in_0 = exp(inst->global->zetah * inst->global->gmh * (v_in_2 - inst->global->vhalfh));
                }
                {
                    double v_in_3;
                    v_in_3 = v_in_5;
                    alph_in_0 = exp(inst->global->zetah * (v_in_3 - inst->global->vhalfh));
                }
                inst->htau[id] = beth_in_0 / (qt * inst->global->a0h * (1.0 + alph_in_0));
            }
            inst->m[id] = inst->m[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->mtau[id]))) * ( -(((inst->minf[id])) / inst->mtau[id]) / (((( -1.0))) / inst->mtau[id]) - inst->m[id]);
            inst->h[id] = inst->h[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->htau[id]))) * ( -(((inst->hinf[id])) / inst->htau[id]) / (((( -1.0))) / inst->htau[id]) - inst->h[id]);
        }
    }


    /** register channel with the simulator */
    void _kamt_reg() {

        int mech_type = nrn_get_mechtype("kamt");
        kamt_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_kamt, nrn_cur_kamt, nullptr, nrn_state_kamt, nrn_init_kamt, nrn_private_constructor_kamt, nrn_private_destructor_kamt, first_pointer_var_index(), 1);
        kamt_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
