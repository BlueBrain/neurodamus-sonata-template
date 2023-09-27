/*********************************************************
Model Name      : Im_ms
Filename        : Im_ms.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Wed Sep 27 13:48:15 2023
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
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "Im_ms",
        "gbar_Im_ms",
        "modACh_Im_ms",
        "maxModACh_Im_ms",
        "levelACh_Im_ms",
        0,
        "ik_Im_ms",
        "gIm_Im_ms",
        0,
        "m_Im_ms",
        0,
        0
    };


    /** all global variables */
    struct Im_ms_Store {
        int k_type{};
        double m0{};
        int reset{};
        int mech_type{};
        int slist1[1]{6};
        int dlist1[1]{12};
    };
    static_assert(std::is_trivially_copy_constructible_v<Im_ms_Store>);
    static_assert(std::is_trivially_move_constructible_v<Im_ms_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Im_ms_Store>);
    static_assert(std::is_trivially_move_assignable_v<Im_ms_Store>);
    static_assert(std::is_trivially_destructible_v<Im_ms_Store>);
    Im_ms_Store Im_ms_global;


    /** all mechanism instance variables and global variables */
    struct Im_ms_Instance  {
        const double* gbar{};
        const double* modACh{};
        const double* maxModACh{};
        const double* levelACh{};
        double* ik{};
        double* gIm{};
        double* m{};
        double* ek{};
        double* mInf{};
        double* mTau{};
        double* mAlpha{};
        double* mBeta{};
        double* Dm{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Im_ms_Store* global{&Im_ms_global};
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
        return Im_ms_global.mech_type;
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
    static void nrn_private_constructor_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Im_ms_Instance{};
        assert(inst->global == &Im_ms_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Im_ms_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Im_ms_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Im_ms_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Im_ms_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Im_ms_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gbar = ml->data+0*pnodecount;
        inst->modACh = ml->data+1*pnodecount;
        inst->maxModACh = ml->data+2*pnodecount;
        inst->levelACh = ml->data+3*pnodecount;
        inst->ik = ml->data+4*pnodecount;
        inst->gIm = ml->data+5*pnodecount;
        inst->m = ml->data+6*pnodecount;
        inst->ek = ml->data+7*pnodecount;
        inst->mInf = ml->data+8*pnodecount;
        inst->mTau = ml->data+9*pnodecount;
        inst->mAlpha = ml->data+10*pnodecount;
        inst->mBeta = ml->data+11*pnodecount;
        inst->Dm = ml->data+12*pnodecount;
        inst->v_unused = ml->data+13*pnodecount;
        inst->g_unused = ml->data+14*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Im_ms(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);

        #endif
    }


    inline double modulationACh_Im_ms(int id, int pnodecount, Im_ms_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int rates_Im_ms(int id, int pnodecount, Im_ms_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);


    inline int rates_Im_ms(int id, int pnodecount, Im_ms_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_rates = 0;
        double qt;
        qt = pow(2.3, ((34.0 - 21.0) / 10.0));
        inst->mAlpha[id] = 3.3e-3 * exp(2.5 * 0.04 * (v -  -35.0));
        inst->mBeta[id] = 3.3e-3 * exp( -2.5 * 0.04 * (v -  -35.0));
        inst->mInf[id] = inst->mAlpha[id] / (inst->mAlpha[id] + inst->mBeta[id]);
        inst->mTau[id] = (1.0 / (inst->mAlpha[id] + inst->mBeta[id])) / qt;
        return ret_rates;
    }


    inline double modulationACh_Im_ms(int id, int pnodecount, Im_ms_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_modulationACh = 0.0;
        ret_modulationACh = 1.0 + inst->modACh[id] * (inst->maxModACh[id] - 1.0) * inst->levelACh[id];
        return ret_modulationACh;
    }


    /** initialize channel */
    void nrn_init_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma ivdep
            #pragma omp simd
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
                inst->m[id] = inst->global->m0;
                {
                    double qt;
                    qt = pow(2.3, ((34.0 - 21.0) / 10.0));
                    inst->mAlpha[id] = 3.3e-3 * exp(2.5 * 0.04 * (v -  -35.0));
                    inst->mBeta[id] = 3.3e-3 * exp( -2.5 * 0.04 * (v -  -35.0));
                    inst->mInf[id] = inst->mAlpha[id] / (inst->mAlpha[id] + inst->mBeta[id]);
                    inst->mTau[id] = (1.0 / (inst->mAlpha[id] + inst->mBeta[id])) / qt;
                }
                inst->m[id] = inst->mInf[id];
            }
        }
    }


    inline double nrn_current_Im_ms(int id, int pnodecount, Im_ms_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double modulationACh_in_0;
        {
            modulationACh_in_0 = 1.0 + inst->modACh[id] * (inst->maxModACh[id] - 1.0) * inst->levelACh[id];
        }
        inst->gIm[id] = inst->gbar[id] * inst->m[id] * modulationACh_in_0;
        inst->ik[id] = inst->gIm[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);

        #pragma ivdep
        #pragma omp simd
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_Im_ms(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Im_ms(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Im_ms(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Im_ms_Instance*>(ml->instance);

        #pragma ivdep
        #pragma omp simd
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            {
                double qt;
                qt = pow(2.3, ((34.0 - 21.0) / 10.0));
                inst->mAlpha[id] = 3.3e-3 * exp(2.5 * 0.04 * (v -  -35.0));
                inst->mBeta[id] = 3.3e-3 * exp( -2.5 * 0.04 * (v -  -35.0));
                inst->mInf[id] = inst->mAlpha[id] / (inst->mAlpha[id] + inst->mBeta[id]);
                inst->mTau[id] = (1.0 / (inst->mAlpha[id] + inst->mBeta[id])) / qt;
            }
            inst->m[id] = inst->m[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->mTau[id]))) * ( -(((inst->mInf[id])) / inst->mTau[id]) / (((( -1.0))) / inst->mTau[id]) - inst->m[id]);
        }
    }


    /** register channel with the simulator */
    void _Im_ms_reg() {

        int mech_type = nrn_get_mechtype("Im_ms");
        Im_ms_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Im_ms, nrn_cur_Im_ms, nullptr, nrn_state_Im_ms, nrn_init_Im_ms, nrn_private_constructor_Im_ms, nrn_private_destructor_Im_ms, first_pointer_var_index(), 1);
        Im_ms_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
