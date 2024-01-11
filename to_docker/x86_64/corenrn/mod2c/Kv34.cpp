/*********************************************************
Model Name      : Kv3_4
Filename        : Kv34.mod
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
        "Kv3_4",
        "gkbar_Kv3_4",
        0,
        "ik_Kv3_4",
        "minf_Kv3_4",
        "mtau_Kv3_4",
        "hinf_Kv3_4",
        "htau_Kv3_4",
        0,
        "m_Kv3_4",
        "h_Kv3_4",
        0,
        0
    };


    /** all global variables */
    struct Kv3_4_Store {
        int k_type{};
        double m0{};
        double h0{};
        int reset{};
        int mech_type{};
        double mivh{-24};
        double mik{15.4};
        double mty0{0.00012851};
        double mtvh1{100.7};
        double mtk1{12.9};
        double mtvh2{-56};
        double mtk2{-23.1};
        double hiy0{0.31};
        double hiA{0.69};
        double hivh{-5.802};
        double hik{11.2};
        double q10{3};
        int slist1[2]{6, 7};
        int dlist1[2]{10, 11};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kv3_4_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kv3_4_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kv3_4_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kv3_4_Store>);
    static_assert(std::is_trivially_destructible_v<Kv3_4_Store>);
    Kv3_4_Store Kv3_4_global;


    /** all mechanism instance variables and global variables */
    struct Kv3_4_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gkbar{};
        double* ik{};
        double* minf{};
        double* mtau{};
        double* hinf{};
        double* htau{};
        double* m{};
        double* h{};
        double* ek{};
        double* qt{};
        double* Dm{};
        double* Dh{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Kv3_4_Store* global{&Kv3_4_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"mivh_Kv3_4", &Kv3_4_global.mivh},
        {"mik_Kv3_4", &Kv3_4_global.mik},
        {"mty0_Kv3_4", &Kv3_4_global.mty0},
        {"mtvh1_Kv3_4", &Kv3_4_global.mtvh1},
        {"mtk1_Kv3_4", &Kv3_4_global.mtk1},
        {"mtvh2_Kv3_4", &Kv3_4_global.mtvh2},
        {"mtk2_Kv3_4", &Kv3_4_global.mtk2},
        {"hiy0_Kv3_4", &Kv3_4_global.hiy0},
        {"hiA_Kv3_4", &Kv3_4_global.hiA},
        {"hivh_Kv3_4", &Kv3_4_global.hivh},
        {"hik_Kv3_4", &Kv3_4_global.hik},
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
        return Kv3_4_global.mech_type;
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
    static void nrn_private_constructor_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kv3_4_Instance{};
        assert(inst->global == &Kv3_4_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kv3_4_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv3_4_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv3_4_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv3_4_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv3_4_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gkbar = ml->data+0*pnodecount;
        inst->ik = ml->data+1*pnodecount;
        inst->minf = ml->data+2*pnodecount;
        inst->mtau = ml->data+3*pnodecount;
        inst->hinf = ml->data+4*pnodecount;
        inst->htau = ml->data+5*pnodecount;
        inst->m = ml->data+6*pnodecount;
        inst->h = ml->data+7*pnodecount;
        inst->ek = ml->data+8*pnodecount;
        inst->qt = ml->data+9*pnodecount;
        inst->Dm = ml->data+10*pnodecount;
        inst->Dh = ml->data+11*pnodecount;
        inst->v_unused = ml->data+12*pnodecount;
        inst->g_unused = ml->data+13*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Kv3_4(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);

        #endif
    }


    inline double mtau_func_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double htau_func_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double Vm);
    inline int rates_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double Vm);


    inline int rates_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double Vm) {
        int ret_rates = 0;
        double v_r_0, mtau_func_in_0, htau_func_in_0;
        v_r_0 = Vm + 11.0;
        inst->minf[id] = 1.0 / (1.0 + exp( -(v_r_0 - inst->global->mivh) / inst->global->mik));
        {
            double v_in_0;
            v_in_0 = v_r_0;
            if (v_in_0 <  -35.0) {
                mtau_func_in_0 = (3.4225e-5 + .00498 * exp( -v_in_0 /  -28.29)) * 3.0;
            } else {
                mtau_func_in_0 = (inst->global->mty0 + 1.0 / (exp((v_in_0 + inst->global->mtvh1) / inst->global->mtk1) + exp((v_in_0 + inst->global->mtvh2) / inst->global->mtk2)));
            }
        }
        inst->mtau[id] = (1000.0) * mtau_func_in_0 / inst->qt[id];
        inst->hinf[id] = inst->global->hiy0 + inst->global->hiA / (1.0 + exp((v_r_0 - inst->global->hivh) / inst->global->hik));
        {
            double Vm_in_0;
            Vm_in_0 = v_r_0;
            if (Vm_in_0 > 0.0) {
                htau_func_in_0 = .0012 + .0023 * exp( -.141 * Vm_in_0);
            } else {
                htau_func_in_0 = 1.2202e-05 + .012 * exp( -pow(((Vm_in_0 - ( -56.3)) / 49.6), 2.0));
            }
        }
        inst->htau[id] = 1000.0 * htau_func_in_0 / inst->qt[id];
        return ret_rates;
    }


    inline double mtau_func_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_mtau_func = 0.0;
        if (arg_v <  -35.0) {
            ret_mtau_func = (3.4225e-5 + .00498 * exp( -arg_v /  -28.29)) * 3.0;
        } else {
            ret_mtau_func = (inst->global->mty0 + 1.0 / (exp((arg_v + inst->global->mtvh1) / inst->global->mtk1) + exp((arg_v + inst->global->mtvh2) / inst->global->mtk2)));
        }
        return ret_mtau_func;
    }


    inline double htau_func_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double Vm) {
        double ret_htau_func = 0.0;
        if (Vm > 0.0) {
            ret_htau_func = .0012 + .0023 * exp( -.141 * Vm);
        } else {
            ret_htau_func = 1.2202e-05 + .012 * exp( -pow(((Vm - ( -56.3)) / 49.6), 2.0));
        }
        return ret_htau_func;
    }


    /** initialize channel */
    void nrn_init_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);

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
                    double v_r_0, mtau_func_in_0, htau_func_in_0, Vm_in_1;
                    Vm_in_1 = v;
                    v_r_0 = Vm_in_1 + 11.0;
                    inst->minf[id] = 1.0 / (1.0 + exp( -(v_r_0 - inst->global->mivh) / inst->global->mik));
                    {
                        double v_in_0;
                        v_in_0 = v_r_0;
                        if (v_in_0 <  -35.0) {
                            mtau_func_in_0 = (3.4225e-5 + .00498 * exp( -v_in_0 /  -28.29)) * 3.0;
                        } else {
                            mtau_func_in_0 = (inst->global->mty0 + 1.0 / (exp((v_in_0 + inst->global->mtvh1) / inst->global->mtk1) + exp((v_in_0 + inst->global->mtvh2) / inst->global->mtk2)));
                        }
                    }
                    inst->mtau[id] = (1000.0) * mtau_func_in_0 / inst->qt[id];
                    inst->hinf[id] = inst->global->hiy0 + inst->global->hiA / (1.0 + exp((v_r_0 - inst->global->hivh) / inst->global->hik));
                    {
                        double Vm_in_0;
                        Vm_in_0 = v_r_0;
                        if (Vm_in_0 > 0.0) {
                            htau_func_in_0 = .0012 + .0023 * exp( -.141 * Vm_in_0);
                        } else {
                            htau_func_in_0 = 1.2202e-05 + .012 * exp( -pow(((Vm_in_0 - ( -56.3)) / 49.6), 2.0));
                        }
                    }
                    inst->htau[id] = 1000.0 * htau_func_in_0 / inst->qt[id];
                }
                inst->m[id] = inst->minf[id];
                inst->h[id] = inst->hinf[id];
                inst->qt[id] = pow(inst->global->q10, ((*(inst->celsius) - 37.0) / 10.0));
            }
        }
    }


    inline double nrn_current_Kv3_4(int id, int pnodecount, Kv3_4_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->ik[id] = inst->gkbar[id] * pow(inst->m[id], 3.0) * inst->h[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_Kv3_4(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kv3_4(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Kv3_4(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv3_4_Instance*>(ml->instance);

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
                double v_r_0, mtau_func_in_0, htau_func_in_0, Vm_in_2;
                Vm_in_2 = v;
                v_r_0 = Vm_in_2 + 11.0;
                inst->minf[id] = 1.0 / (1.0 + exp( -(v_r_0 - inst->global->mivh) / inst->global->mik));
                {
                    double v_in_0;
                    v_in_0 = v_r_0;
                    if (v_in_0 <  -35.0) {
                        mtau_func_in_0 = (3.4225e-5 + .00498 * exp( -v_in_0 /  -28.29)) * 3.0;
                    } else {
                        mtau_func_in_0 = (inst->global->mty0 + 1.0 / (exp((v_in_0 + inst->global->mtvh1) / inst->global->mtk1) + exp((v_in_0 + inst->global->mtvh2) / inst->global->mtk2)));
                    }
                }
                inst->mtau[id] = (1000.0) * mtau_func_in_0 / inst->qt[id];
                inst->hinf[id] = inst->global->hiy0 + inst->global->hiA / (1.0 + exp((v_r_0 - inst->global->hivh) / inst->global->hik));
                {
                    double Vm_in_0;
                    Vm_in_0 = v_r_0;
                    if (Vm_in_0 > 0.0) {
                        htau_func_in_0 = .0012 + .0023 * exp( -.141 * Vm_in_0);
                    } else {
                        htau_func_in_0 = 1.2202e-05 + .012 * exp( -pow(((Vm_in_0 - ( -56.3)) / 49.6), 2.0));
                    }
                }
                inst->htau[id] = 1000.0 * htau_func_in_0 / inst->qt[id];
            }
            inst->m[id] = inst->m[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->mtau[id]))) * ( -(((inst->minf[id])) / inst->mtau[id]) / (((( -1.0))) / inst->mtau[id]) - inst->m[id]);
            inst->h[id] = inst->h[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->htau[id]))) * ( -(((inst->hinf[id])) / inst->htau[id]) / (((( -1.0))) / inst->htau[id]) - inst->h[id]);
        }
    }


    /** register channel with the simulator */
    void _Kv34_reg() {

        int mech_type = nrn_get_mechtype("Kv3_4");
        Kv3_4_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kv3_4, nrn_cur_Kv3_4, nullptr, nrn_state_Kv3_4, nrn_init_Kv3_4, nrn_private_constructor_Kv3_4, nrn_private_destructor_Kv3_4, first_pointer_var_index(), 1);
        Kv3_4_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
