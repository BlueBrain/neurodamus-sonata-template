/*********************************************************
Model Name      : Kv1_1
Filename        : Kv11.mod
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
        "Kv1_1",
        "gateCurrent_Kv1_1",
        "gbar_Kv1_1",
        "gunit_Kv1_1",
        0,
        "ik_Kv1_1",
        "i_Kv1_1",
        "igate_Kv1_1",
        "g_Kv1_1",
        "nc_Kv1_1",
        "ninf_Kv1_1",
        "taun_Kv1_1",
        0,
        "n_Kv1_1",
        0,
        0
    };


    /** all global variables */
    struct Kv1_1_Store {
        int k_type{};
        double n0{};
        int reset{};
        int mech_type{};
        double e0{1.60218e-19};
        double q10{2.7};
        double ca{0.12889};
        double cva{45};
        double cka{-33.9088};
        double cb{0.12889};
        double cvb{45};
        double ckb{12.421};
        double zn{2.7978};
        int slist1[1]{10};
        int dlist1[1]{15};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kv1_1_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kv1_1_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kv1_1_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kv1_1_Store>);
    static_assert(std::is_trivially_destructible_v<Kv1_1_Store>);
    Kv1_1_Store Kv1_1_global;


    /** all mechanism instance variables and global variables */
    struct Kv1_1_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gateCurrent{};
        const double* gbar{};
        const double* gunit{};
        double* ik{};
        double* i{};
        double* igate{};
        double* g{};
        double* nc{};
        double* ninf{};
        double* taun{};
        double* n{};
        double* ek{};
        double* alphan{};
        double* betan{};
        double* qt{};
        double* Dn{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Kv1_1_Store* global{&Kv1_1_global};
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
        return 18;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return Kv1_1_global.mech_type;
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
    static void nrn_private_constructor_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kv1_1_Instance{};
        assert(inst->global == &Kv1_1_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kv1_1_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv1_1_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv1_1_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv1_1_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv1_1_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gateCurrent = ml->data+0*pnodecount;
        inst->gbar = ml->data+1*pnodecount;
        inst->gunit = ml->data+2*pnodecount;
        inst->ik = ml->data+3*pnodecount;
        inst->i = ml->data+4*pnodecount;
        inst->igate = ml->data+5*pnodecount;
        inst->g = ml->data+6*pnodecount;
        inst->nc = ml->data+7*pnodecount;
        inst->ninf = ml->data+8*pnodecount;
        inst->taun = ml->data+9*pnodecount;
        inst->n = ml->data+10*pnodecount;
        inst->ek = ml->data+11*pnodecount;
        inst->alphan = ml->data+12*pnodecount;
        inst->betan = ml->data+13*pnodecount;
        inst->qt = ml->data+14*pnodecount;
        inst->Dn = ml->data+15*pnodecount;
        inst->v_unused = ml->data+16*pnodecount;
        inst->g_unused = ml->data+17*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Kv1_1(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);

        #endif
    }


    inline double alphanfkt_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double betanfkt_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double ngateFlip_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int rates_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int rates_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_rates = 0;
        double alphanfkt_in_0, betanfkt_in_0;
        {
            double v_in_0;
            v_in_0 = arg_v;
            alphanfkt_in_0 = inst->global->ca * exp( -(v_in_0 + inst->global->cva) / inst->global->cka);
        }
        inst->alphan[id] = alphanfkt_in_0;
        {
            double v_in_1;
            v_in_1 = arg_v;
            betanfkt_in_0 = inst->global->cb * exp( -(v_in_1 + inst->global->cvb) / inst->global->ckb);
        }
        inst->betan[id] = betanfkt_in_0;
        inst->ninf[id] = inst->alphan[id] / (inst->alphan[id] + inst->betan[id]);
        inst->taun[id] = 1.0 / (inst->qt[id] * (inst->alphan[id] + inst->betan[id]));
        return ret_rates;
    }


    inline double alphanfkt_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alphanfkt = 0.0;
        ret_alphanfkt = inst->global->ca * exp( -(arg_v + inst->global->cva) / inst->global->cka);
        return ret_alphanfkt;
    }


    inline double betanfkt_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_betanfkt = 0.0;
        ret_betanfkt = inst->global->cb * exp( -(arg_v + inst->global->cvb) / inst->global->ckb);
        return ret_betanfkt;
    }


    inline double ngateFlip_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_ngateFlip = 0.0;
        ret_ngateFlip = (inst->ninf[id] - inst->n[id]) / inst->taun[id];
        return ret_ngateFlip;
    }


    /** initialize channel */
    void nrn_init_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);

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
                inst->n[id] = inst->global->n0;
                inst->nc[id] = (1e12) * inst->gbar[id] / inst->gunit[id];
                inst->qt[id] = pow(inst->global->q10, ((*(inst->celsius) - 22.0) / 10.0));
                {
                    double alphanfkt_in_0, betanfkt_in_0, v_in_2;
                    v_in_2 = v;
                    {
                        double v_in_0;
                        v_in_0 = v_in_2;
                        alphanfkt_in_0 = inst->global->ca * exp( -(v_in_0 + inst->global->cva) / inst->global->cka);
                    }
                    inst->alphan[id] = alphanfkt_in_0;
                    {
                        double v_in_1;
                        v_in_1 = v_in_2;
                        betanfkt_in_0 = inst->global->cb * exp( -(v_in_1 + inst->global->cvb) / inst->global->ckb);
                    }
                    inst->betan[id] = betanfkt_in_0;
                    inst->ninf[id] = inst->alphan[id] / (inst->alphan[id] + inst->betan[id]);
                    inst->taun[id] = 1.0 / (inst->qt[id] * (inst->alphan[id] + inst->betan[id]));
                }
                inst->n[id] = inst->ninf[id];
            }
        }
    }


    inline double nrn_current_Kv1_1(int id, int pnodecount, Kv1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double ngateFlip_in_0;
        inst->g[id] = inst->gbar[id] * pow(inst->n[id], 4.0);
        inst->ik[id] = inst->g[id] * (v - inst->ek[id]);
        {
            ngateFlip_in_0 = (inst->ninf[id] - inst->n[id]) / inst->taun[id];
        }
        inst->igate[id] = inst->nc[id] * (1e6) * inst->global->e0 * 4.0 * inst->global->zn * ngateFlip_in_0;
        if (inst->gateCurrent[id] != 0.0) {
            inst->i[id] = inst->igate[id];
        }
        current += inst->i[id];
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_Kv1_1(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kv1_1(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Kv1_1(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_1_Instance*>(ml->instance);

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
                double alphanfkt_in_0, betanfkt_in_0, v_in_3;
                v_in_3 = v;
                {
                    double v_in_0;
                    v_in_0 = v_in_3;
                    alphanfkt_in_0 = inst->global->ca * exp( -(v_in_0 + inst->global->cva) / inst->global->cka);
                }
                inst->alphan[id] = alphanfkt_in_0;
                {
                    double v_in_1;
                    v_in_1 = v_in_3;
                    betanfkt_in_0 = inst->global->cb * exp( -(v_in_1 + inst->global->cvb) / inst->global->ckb);
                }
                inst->betan[id] = betanfkt_in_0;
                inst->ninf[id] = inst->alphan[id] / (inst->alphan[id] + inst->betan[id]);
                inst->taun[id] = 1.0 / (inst->qt[id] * (inst->alphan[id] + inst->betan[id]));
            }
            inst->n[id] = inst->n[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->taun[id]))) * ( -(((inst->ninf[id])) / inst->taun[id]) / (((( -1.0))) / inst->taun[id]) - inst->n[id]);
        }
    }


    /** register channel with the simulator */
    void _Kv11_reg() {

        int mech_type = nrn_get_mechtype("Kv1_1");
        Kv1_1_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kv1_1, nrn_cur_Kv1_1, nullptr, nrn_state_Kv1_1, nrn_init_Kv1_1, nrn_private_constructor_Kv1_1, nrn_private_destructor_Kv1_1, first_pointer_var_index(), 1);
        Kv1_1_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
