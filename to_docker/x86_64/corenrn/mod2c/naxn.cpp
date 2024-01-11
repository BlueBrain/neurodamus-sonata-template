/*********************************************************
Model Name      : nax
Filename        : naxn.mod
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
        "nax",
        "sh_nax",
        "gbar_nax",
        0,
        "minf_nax",
        "hinf_nax",
        "mtau_nax",
        "htau_nax",
        0,
        "m_nax",
        "h_nax",
        0,
        0
    };


    /** all global variables */
    struct nax_Store {
        int na_type{};
        double m0{};
        double h0{};
        int reset{};
        int mech_type{};
        double thinf{-50};
        double qinf{4};
        double tha{-30};
        double qa{7.2};
        double Ra{0.4};
        double Rb{0.124};
        double thi1{-45};
        double thi2{-45};
        double qd{1.5};
        double qg{1.5};
        double mmin{0.02};
        double hmin{0.5};
        double q10{2};
        double Rg{0.01};
        double Rd{0.03};
        int slist1[2]{6, 7};
        int dlist1[2]{11, 12};
    };
    static_assert(std::is_trivially_copy_constructible_v<nax_Store>);
    static_assert(std::is_trivially_move_constructible_v<nax_Store>);
    static_assert(std::is_trivially_copy_assignable_v<nax_Store>);
    static_assert(std::is_trivially_move_assignable_v<nax_Store>);
    static_assert(std::is_trivially_destructible_v<nax_Store>);
    nax_Store nax_global;


    /** all mechanism instance variables and global variables */
    struct nax_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* sh{};
        const double* gbar{};
        double* minf{};
        double* hinf{};
        double* mtau{};
        double* htau{};
        double* m{};
        double* h{};
        double* ena{};
        double* ina{};
        double* thegna{};
        double* Dm{};
        double* Dh{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ena{};
        double* ion_ina{};
        double* ion_dinadv{};
        nax_Store* global{&nax_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"thinf_nax", &nax_global.thinf},
        {"qinf_nax", &nax_global.qinf},
        {"tha_nax", &nax_global.tha},
        {"qa_nax", &nax_global.qa},
        {"Ra_nax", &nax_global.Ra},
        {"Rb_nax", &nax_global.Rb},
        {"thi1_nax", &nax_global.thi1},
        {"thi2_nax", &nax_global.thi2},
        {"qd_nax", &nax_global.qd},
        {"qg_nax", &nax_global.qg},
        {"mmin_nax", &nax_global.mmin},
        {"hmin_nax", &nax_global.hmin},
        {"q10_nax", &nax_global.q10},
        {"Rg_nax", &nax_global.Rg},
        {"Rd_nax", &nax_global.Rd},
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
        return nax_global.mech_type;
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
    static void nrn_private_constructor_nax(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new nax_Instance{};
        assert(inst->global == &nax_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(nax_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_nax(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<nax_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &nax_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(nax_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<nax_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &nax_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(nax_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->sh = ml->data+0*pnodecount;
        inst->gbar = ml->data+1*pnodecount;
        inst->minf = ml->data+2*pnodecount;
        inst->hinf = ml->data+3*pnodecount;
        inst->mtau = ml->data+4*pnodecount;
        inst->htau = ml->data+5*pnodecount;
        inst->m = ml->data+6*pnodecount;
        inst->h = ml->data+7*pnodecount;
        inst->ena = ml->data+8*pnodecount;
        inst->ina = ml->data+9*pnodecount;
        inst->thegna = ml->data+10*pnodecount;
        inst->Dm = ml->data+11*pnodecount;
        inst->Dh = ml->data+12*pnodecount;
        inst->v_unused = ml->data+13*pnodecount;
        inst->g_unused = ml->data+14*pnodecount;
        inst->ion_ena = nt->_data;
        inst->ion_ina = nt->_data;
        inst->ion_dinadv = nt->_data;
    }



    static void nrn_alloc_nax(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_nax(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<nax_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_nax(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<nax_Instance*>(ml->instance);

        #endif
    }


    inline double trap0_nax(int id, int pnodecount, nax_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double th, double a, double q);
    inline int trates_nax(int id, int pnodecount, nax_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double vm, double sh2);


    inline int trates_nax(int id, int pnodecount, nax_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double vm, double sh2) {
        int ret_trates = 0;
        double a, b, qt, trap0_in_0, trap0_in_1, trap0_in_2, trap0_in_3;
        qt = pow(inst->global->q10, ((*(inst->celsius) - 24.0) / 10.0));
        {
            double v_in_0, th_in_0, a_in_0, q_in_0;
            v_in_0 = vm;
            th_in_0 = inst->global->tha + sh2;
            a_in_0 = inst->global->Ra;
            q_in_0 = inst->global->qa;
            if (fabs(v_in_0 - th_in_0) > 1e-6) {
                trap0_in_0 = a_in_0 * (v_in_0 - th_in_0) / (1.0 - exp( -(v_in_0 - th_in_0) / q_in_0));
            } else {
                trap0_in_0 = a_in_0 * q_in_0;
            }
        }
        a = trap0_in_0;
        {
            double v_in_1, th_in_1, a_in_1, q_in_1;
            v_in_1 =  -vm;
            th_in_1 =  -inst->global->tha - sh2;
            a_in_1 = inst->global->Rb;
            q_in_1 = inst->global->qa;
            if (fabs(v_in_1 - th_in_1) > 1e-6) {
                trap0_in_1 = a_in_1 * (v_in_1 - th_in_1) / (1.0 - exp( -(v_in_1 - th_in_1) / q_in_1));
            } else {
                trap0_in_1 = a_in_1 * q_in_1;
            }
        }
        b = trap0_in_1;
        inst->mtau[id] = 1.0 / (a + b) / qt;
        if (inst->mtau[id] < inst->global->mmin) {
            inst->mtau[id] = inst->global->mmin;
        }
        inst->minf[id] = a / (a + b);
        {
            double v_in_2, th_in_2, a_in_2, q_in_2;
            v_in_2 = vm;
            th_in_2 = inst->global->thi1 + sh2;
            a_in_2 = inst->global->Rd;
            q_in_2 = inst->global->qd;
            if (fabs(v_in_2 - th_in_2) > 1e-6) {
                trap0_in_2 = a_in_2 * (v_in_2 - th_in_2) / (1.0 - exp( -(v_in_2 - th_in_2) / q_in_2));
            } else {
                trap0_in_2 = a_in_2 * q_in_2;
            }
        }
        a = trap0_in_2;
        {
            double v_in_3, th_in_3, a_in_3, q_in_3;
            v_in_3 =  -vm;
            th_in_3 =  -inst->global->thi2 - sh2;
            a_in_3 = inst->global->Rg;
            q_in_3 = inst->global->qg;
            if (fabs(v_in_3 - th_in_3) > 1e-6) {
                trap0_in_3 = a_in_3 * (v_in_3 - th_in_3) / (1.0 - exp( -(v_in_3 - th_in_3) / q_in_3));
            } else {
                trap0_in_3 = a_in_3 * q_in_3;
            }
        }
        b = trap0_in_3;
        inst->htau[id] = 1.0 / (a + b) / qt;
        if (inst->htau[id] < inst->global->hmin) {
            inst->htau[id] = inst->global->hmin;
        }
        inst->hinf[id] = 1.0 / (1.0 + exp((vm - inst->global->thinf - sh2) / inst->global->qinf));
        return ret_trates;
    }


    inline double trap0_nax(int id, int pnodecount, nax_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double th, double a, double q) {
        double ret_trap0 = 0.0;
        if (fabs(arg_v - th) > 1e-6) {
            ret_trap0 = a * (arg_v - th) / (1.0 - exp( -(arg_v - th) / q));
        } else {
            ret_trap0 = a * q;
        }
        return ret_trap0;
    }


    /** initialize channel */
    void nrn_init_nax(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<nax_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->ena[id] = inst->ion_ena[indexes[0*pnodecount + id]];
                inst->m[id] = inst->global->m0;
                inst->h[id] = inst->global->h0;
                {
                    double a, b, qt, trap0_in_0, trap0_in_1, trap0_in_2, trap0_in_3, vm_in_0, sh2_in_0;
                    vm_in_0 = v;
                    sh2_in_0 = inst->sh[id];
                    qt = pow(inst->global->q10, ((*(inst->celsius) - 24.0) / 10.0));
                    {
                        double v_in_0, th_in_0, a_in_0, q_in_0;
                        v_in_0 = vm_in_0;
                        th_in_0 = inst->global->tha + sh2_in_0;
                        a_in_0 = inst->global->Ra;
                        q_in_0 = inst->global->qa;
                        if (fabs(v_in_0 - th_in_0) > 1e-6) {
                            trap0_in_0 = a_in_0 * (v_in_0 - th_in_0) / (1.0 - exp( -(v_in_0 - th_in_0) / q_in_0));
                        } else {
                            trap0_in_0 = a_in_0 * q_in_0;
                        }
                    }
                    a = trap0_in_0;
                    {
                        double v_in_1, th_in_1, a_in_1, q_in_1;
                        v_in_1 =  -vm_in_0;
                        th_in_1 =  -inst->global->tha - sh2_in_0;
                        a_in_1 = inst->global->Rb;
                        q_in_1 = inst->global->qa;
                        if (fabs(v_in_1 - th_in_1) > 1e-6) {
                            trap0_in_1 = a_in_1 * (v_in_1 - th_in_1) / (1.0 - exp( -(v_in_1 - th_in_1) / q_in_1));
                        } else {
                            trap0_in_1 = a_in_1 * q_in_1;
                        }
                    }
                    b = trap0_in_1;
                    inst->mtau[id] = 1.0 / (a + b) / qt;
                    if (inst->mtau[id] < inst->global->mmin) {
                        inst->mtau[id] = inst->global->mmin;
                    }
                    inst->minf[id] = a / (a + b);
                    {
                        double v_in_2, th_in_2, a_in_2, q_in_2;
                        v_in_2 = vm_in_0;
                        th_in_2 = inst->global->thi1 + sh2_in_0;
                        a_in_2 = inst->global->Rd;
                        q_in_2 = inst->global->qd;
                        if (fabs(v_in_2 - th_in_2) > 1e-6) {
                            trap0_in_2 = a_in_2 * (v_in_2 - th_in_2) / (1.0 - exp( -(v_in_2 - th_in_2) / q_in_2));
                        } else {
                            trap0_in_2 = a_in_2 * q_in_2;
                        }
                    }
                    a = trap0_in_2;
                    {
                        double v_in_3, th_in_3, a_in_3, q_in_3;
                        v_in_3 =  -vm_in_0;
                        th_in_3 =  -inst->global->thi2 - sh2_in_0;
                        a_in_3 = inst->global->Rg;
                        q_in_3 = inst->global->qg;
                        if (fabs(v_in_3 - th_in_3) > 1e-6) {
                            trap0_in_3 = a_in_3 * (v_in_3 - th_in_3) / (1.0 - exp( -(v_in_3 - th_in_3) / q_in_3));
                        } else {
                            trap0_in_3 = a_in_3 * q_in_3;
                        }
                    }
                    b = trap0_in_3;
                    inst->htau[id] = 1.0 / (a + b) / qt;
                    if (inst->htau[id] < inst->global->hmin) {
                        inst->htau[id] = inst->global->hmin;
                    }
                    inst->hinf[id] = 1.0 / (1.0 + exp((vm_in_0 - inst->global->thinf - sh2_in_0) / inst->global->qinf));
                }
                inst->m[id] = inst->minf[id];
                inst->h[id] = inst->hinf[id];
            }
        }
    }


    inline double nrn_current_nax(int id, int pnodecount, nax_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->thegna[id] = inst->gbar[id] * inst->m[id] * inst->m[id] * inst->m[id] * inst->h[id];
        inst->ina[id] = inst->thegna[id] * (v - inst->ena[id]);
        current += inst->ina[id];
        return current;
    }


    /** update current */
    void nrn_cur_nax(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<nax_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ena[id] = inst->ion_ena[indexes[0*pnodecount + id]];
            double g = nrn_current_nax(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dina = inst->ina[id];
            double rhs = nrn_current_nax(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            inst->ion_dinadv[indexes[2*pnodecount + id]] += (dina-inst->ina[id])/0.001;
            inst->ion_ina[indexes[1*pnodecount + id]] += inst->ina[id];
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            vec_rhs[node_id] -= rhs;
            vec_d[node_id] += g;
        }
    }


    /** update state */
    void nrn_state_nax(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<nax_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ena[id] = inst->ion_ena[indexes[0*pnodecount + id]];
            {
                double a, b, qt, trap0_in_0, trap0_in_1, trap0_in_2, trap0_in_3, vm_in_1, sh2_in_1;
                vm_in_1 = v;
                sh2_in_1 = inst->sh[id];
                qt = pow(inst->global->q10, ((*(inst->celsius) - 24.0) / 10.0));
                {
                    double v_in_0, th_in_0, a_in_0, q_in_0;
                    v_in_0 = vm_in_1;
                    th_in_0 = inst->global->tha + sh2_in_1;
                    a_in_0 = inst->global->Ra;
                    q_in_0 = inst->global->qa;
                    if (fabs(v_in_0 - th_in_0) > 1e-6) {
                        trap0_in_0 = a_in_0 * (v_in_0 - th_in_0) / (1.0 - exp( -(v_in_0 - th_in_0) / q_in_0));
                    } else {
                        trap0_in_0 = a_in_0 * q_in_0;
                    }
                }
                a = trap0_in_0;
                {
                    double v_in_1, th_in_1, a_in_1, q_in_1;
                    v_in_1 =  -vm_in_1;
                    th_in_1 =  -inst->global->tha - sh2_in_1;
                    a_in_1 = inst->global->Rb;
                    q_in_1 = inst->global->qa;
                    if (fabs(v_in_1 - th_in_1) > 1e-6) {
                        trap0_in_1 = a_in_1 * (v_in_1 - th_in_1) / (1.0 - exp( -(v_in_1 - th_in_1) / q_in_1));
                    } else {
                        trap0_in_1 = a_in_1 * q_in_1;
                    }
                }
                b = trap0_in_1;
                inst->mtau[id] = 1.0 / (a + b) / qt;
                if (inst->mtau[id] < inst->global->mmin) {
                    inst->mtau[id] = inst->global->mmin;
                }
                inst->minf[id] = a / (a + b);
                {
                    double v_in_2, th_in_2, a_in_2, q_in_2;
                    v_in_2 = vm_in_1;
                    th_in_2 = inst->global->thi1 + sh2_in_1;
                    a_in_2 = inst->global->Rd;
                    q_in_2 = inst->global->qd;
                    if (fabs(v_in_2 - th_in_2) > 1e-6) {
                        trap0_in_2 = a_in_2 * (v_in_2 - th_in_2) / (1.0 - exp( -(v_in_2 - th_in_2) / q_in_2));
                    } else {
                        trap0_in_2 = a_in_2 * q_in_2;
                    }
                }
                a = trap0_in_2;
                {
                    double v_in_3, th_in_3, a_in_3, q_in_3;
                    v_in_3 =  -vm_in_1;
                    th_in_3 =  -inst->global->thi2 - sh2_in_1;
                    a_in_3 = inst->global->Rg;
                    q_in_3 = inst->global->qg;
                    if (fabs(v_in_3 - th_in_3) > 1e-6) {
                        trap0_in_3 = a_in_3 * (v_in_3 - th_in_3) / (1.0 - exp( -(v_in_3 - th_in_3) / q_in_3));
                    } else {
                        trap0_in_3 = a_in_3 * q_in_3;
                    }
                }
                b = trap0_in_3;
                inst->htau[id] = 1.0 / (a + b) / qt;
                if (inst->htau[id] < inst->global->hmin) {
                    inst->htau[id] = inst->global->hmin;
                }
                inst->hinf[id] = 1.0 / (1.0 + exp((vm_in_1 - inst->global->thinf - sh2_in_1) / inst->global->qinf));
            }
            inst->m[id] = inst->m[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->mtau[id]))) * ( -(((inst->minf[id])) / inst->mtau[id]) / (((( -1.0))) / inst->mtau[id]) - inst->m[id]);
            inst->h[id] = inst->h[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->htau[id]))) * ( -(((inst->hinf[id])) / inst->htau[id]) / (((( -1.0))) / inst->htau[id]) - inst->h[id]);
        }
    }


    /** register channel with the simulator */
    void _naxn_reg() {

        int mech_type = nrn_get_mechtype("nax");
        nax_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_nax, nrn_cur_nax, nullptr, nrn_state_nax, nrn_init_nax, nrn_private_constructor_nax, nrn_private_destructor_nax, first_pointer_var_index(), 1);
        nax_global.na_type = nrn_get_mechtype("na_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "na_ion");
        hoc_register_dparam_semantics(mech_type, 1, "na_ion");
        hoc_register_dparam_semantics(mech_type, 2, "na_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
