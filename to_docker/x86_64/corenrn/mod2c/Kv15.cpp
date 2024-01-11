/*********************************************************
Model Name      : Kv1_5
Filename        : Kv15.mod
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


    /** constants used in nmodl from UNITS */
    static const double F = 0x1.78e555060882cp+16;
    static const double R = 0x1.0a1013e8990bep+3;
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "Kv1_5",
        "gKur_Kv1_5",
        "Tauact_Kv1_5",
        "Tauinactf_Kv1_5",
        "Tauinacts_Kv1_5",
        "gnonspec_Kv1_5",
        0,
        "ik_Kv1_5",
        "minf_Kv1_5",
        "ninf_Kv1_5",
        "uinf_Kv1_5",
        "mtau_Kv1_5",
        "ntau_Kv1_5",
        "utau_Kv1_5",
        "ino_Kv1_5",
        0,
        "m_Kv1_5",
        "n_Kv1_5",
        "u_Kv1_5",
        0,
        0
    };


    /** all global variables */
    struct Kv1_5_Store {
        int k_type{};
        int na_type{};
        int no_type{};
        double m0{};
        double n0{};
        double u0{};
        int reset{};
        int mech_type{};
        int slist1[3]{13, 14, 15};
        int dlist1[3]{16, 17, 18};
        int slist2[3]{13, 14, 15};
        ThreadDatum ext_call_thread[3]{};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kv1_5_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kv1_5_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kv1_5_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kv1_5_Store>);
    static_assert(std::is_trivially_destructible_v<Kv1_5_Store>);
    Kv1_5_Store Kv1_5_global;


    /** all mechanism instance variables and global variables */
    struct Kv1_5_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gKur{};
        const double* Tauact{};
        const double* Tauinactf{};
        const double* Tauinacts{};
        const double* gnonspec{};
        double* ik{};
        double* minf{};
        double* ninf{};
        double* uinf{};
        double* mtau{};
        double* ntau{};
        double* utau{};
        double* ino{};
        double* m{};
        double* n{};
        double* u{};
        double* Dm{};
        double* Dn{};
        double* Du{};
        double* ek{};
        double* ki{};
        double* ko{};
        double* nai{};
        double* nao{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        const double* ion_ki{};
        const double* ion_ko{};
        double* ion_ik{};
        double* ion_dikdv{};
        const double* ion_nai{};
        const double* ion_nao{};
        double* ion_ino{};
        double* ion_dinodv{};
        Kv1_5_Store* global{&Kv1_5_global};
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


    /** thread specific helper routines for derivimplicit */

    static inline int* deriv1_advance(ThreadDatum* thread) {
        return &(thread[0].i);
    }

    static inline int dith1() {
        return 1;
    }

    static inline void** newtonspace1(ThreadDatum* thread) {
        return &(thread[2]._pvoid);
    }


    static inline int float_variables_size() {
        return 26;
    }


    static inline int int_variables_size() {
        return 9;
    }


    static inline int get_mech_type() {
        return Kv1_5_global.mech_type;
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


    /** thread memory allocation callback */
    static void thread_mem_init(ThreadDatum* thread)  {
        thread[dith1()].pval = nullptr;
    }


    /** thread memory cleanup callback */
    static void thread_mem_cleanup(ThreadDatum* thread)  {
        free(thread[dith1()].pval);
        nrn_destroy_newtonspace(static_cast<NewtonSpace*>(*newtonspace1(thread)));
    }

    // Allocate instance structure
    static void nrn_private_constructor_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kv1_5_Instance{};
        assert(inst->global == &Kv1_5_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kv1_5_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv1_5_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv1_5_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv1_5_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv1_5_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gKur = ml->data+0*pnodecount;
        inst->Tauact = ml->data+1*pnodecount;
        inst->Tauinactf = ml->data+2*pnodecount;
        inst->Tauinacts = ml->data+3*pnodecount;
        inst->gnonspec = ml->data+4*pnodecount;
        inst->ik = ml->data+5*pnodecount;
        inst->minf = ml->data+6*pnodecount;
        inst->ninf = ml->data+7*pnodecount;
        inst->uinf = ml->data+8*pnodecount;
        inst->mtau = ml->data+9*pnodecount;
        inst->ntau = ml->data+10*pnodecount;
        inst->utau = ml->data+11*pnodecount;
        inst->ino = ml->data+12*pnodecount;
        inst->m = ml->data+13*pnodecount;
        inst->n = ml->data+14*pnodecount;
        inst->u = ml->data+15*pnodecount;
        inst->Dm = ml->data+16*pnodecount;
        inst->Dn = ml->data+17*pnodecount;
        inst->Du = ml->data+18*pnodecount;
        inst->ek = ml->data+19*pnodecount;
        inst->ki = ml->data+20*pnodecount;
        inst->ko = ml->data+21*pnodecount;
        inst->nai = ml->data+22*pnodecount;
        inst->nao = ml->data+23*pnodecount;
        inst->v_unused = ml->data+24*pnodecount;
        inst->g_unused = ml->data+25*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ki = nt->_data;
        inst->ion_ko = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
        inst->ion_nai = nt->_data;
        inst->ion_nao = nt->_data;
        inst->ion_ino = nt->_data;
        inst->ion_dinodv = nt->_data;
    }



    static void nrn_alloc_Kv1_5(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);

        #endif
    }


    inline double alp_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double i);
    inline double bet_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double i);
    inline double ce_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double i);
    inline int rates_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int rates_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_rates = 0;
        double a, b, c, alp_in_0, bet_in_0, ce_in_0, alp_in_1, bet_in_1, ce_in_1, ce_in_2;
        {
            double q10, v_in_0, i_in_0;
            v_in_0 = arg_v;
            i_in_0 = 0.0;
            v_in_0 = v_in_0;
            q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
            if (i_in_0 == 0.0) {
                alp_in_0 = q10 * 0.65 / (exp( -(v_in_0 + 10.0) / 8.5) + exp( -(v_in_0 - 30.0) / 59.0));
            } else if (i_in_0 == 1.0) {
                alp_in_0 = 0.001 * q10 / (2.4 + 10.9 * exp( -(v_in_0 + 90.0) / 78.0));
            }
        }
        a = alp_in_0;
        {
            double q10, v_in_1, i_in_1;
            v_in_1 = arg_v;
            i_in_1 = 0.0;
            v_in_1 = v_in_1;
            q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
            if (i_in_1 == 0.0) {
                bet_in_0 = q10 * 0.65 / (2.5 + exp((v_in_1 + 82.0) / 17.0));
            } else if (i_in_1 == 1.0) {
                bet_in_0 = q10 * 0.001 * exp((v_in_1 - 168.0) / 16.0);
            }
        }
        b = bet_in_0;
        {
            double v_in_2, i_in_2;
            v_in_2 = arg_v;
            i_in_2 = 0.0;
            v_in_2 = v_in_2;
            if (i_in_2 == 0.0) {
                ce_in_0 = 1.0 / (1.0 + exp( -(v_in_2 + 30.3) / 9.6));
            } else if (i_in_2 == 1.0) {
                ce_in_0 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_2 + 7.0) / 14.0)));
            } else if (i_in_2 == 2.0) {
                ce_in_0 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_2 + 7.0) / 14.0)));
            }
        }
        c = ce_in_0;
        inst->mtau[id] = 1.0 / (a + b) / 3.0 * inst->Tauact[id];
        inst->minf[id] = c;
        {
            double q10, v_in_3, i_in_3;
            v_in_3 = arg_v;
            i_in_3 = 1.0;
            v_in_3 = v_in_3;
            q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
            if (i_in_3 == 0.0) {
                alp_in_1 = q10 * 0.65 / (exp( -(v_in_3 + 10.0) / 8.5) + exp( -(v_in_3 - 30.0) / 59.0));
            } else if (i_in_3 == 1.0) {
                alp_in_1 = 0.001 * q10 / (2.4 + 10.9 * exp( -(v_in_3 + 90.0) / 78.0));
            }
        }
        a = alp_in_1;
        {
            double q10, v_in_4, i_in_4;
            v_in_4 = arg_v;
            i_in_4 = 1.0;
            v_in_4 = v_in_4;
            q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
            if (i_in_4 == 0.0) {
                bet_in_1 = q10 * 0.65 / (2.5 + exp((v_in_4 + 82.0) / 17.0));
            } else if (i_in_4 == 1.0) {
                bet_in_1 = q10 * 0.001 * exp((v_in_4 - 168.0) / 16.0);
            }
        }
        b = bet_in_1;
        {
            double v_in_5, i_in_5;
            v_in_5 = arg_v;
            i_in_5 = 1.0;
            v_in_5 = v_in_5;
            if (i_in_5 == 0.0) {
                ce_in_1 = 1.0 / (1.0 + exp( -(v_in_5 + 30.3) / 9.6));
            } else if (i_in_5 == 1.0) {
                ce_in_1 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_5 + 7.0) / 14.0)));
            } else if (i_in_5 == 2.0) {
                ce_in_1 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_5 + 7.0) / 14.0)));
            }
        }
        c = ce_in_1;
        inst->ntau[id] = 1.0 / (a + b) / 3.0 * inst->Tauinactf[id];
        inst->ninf[id] = c;
        {
            double v_in_6, i_in_6;
            v_in_6 = arg_v;
            i_in_6 = 2.0;
            v_in_6 = v_in_6;
            if (i_in_6 == 0.0) {
                ce_in_2 = 1.0 / (1.0 + exp( -(v_in_6 + 30.3) / 9.6));
            } else if (i_in_6 == 1.0) {
                ce_in_2 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_6 + 7.0) / 14.0)));
            } else if (i_in_6 == 2.0) {
                ce_in_2 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_6 + 7.0) / 14.0)));
            }
        }
        c = ce_in_2;
        inst->uinf[id] = c;
        inst->utau[id] = 6800.0 * inst->Tauinacts[id];
        return ret_rates;
    }


    inline double alp_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double i) {
        double ret_alp = 0.0;
        double q10;
        arg_v = arg_v;
        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
        if (i == 0.0) {
            ret_alp = q10 * 0.65 / (exp( -(arg_v + 10.0) / 8.5) + exp( -(arg_v - 30.0) / 59.0));
        } else if (i == 1.0) {
            ret_alp = 0.001 * q10 / (2.4 + 10.9 * exp( -(arg_v + 90.0) / 78.0));
        }
        return ret_alp;
    }


    inline double bet_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double i) {
        double ret_bet = 0.0;
        double q10;
        arg_v = arg_v;
        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
        if (i == 0.0) {
            ret_bet = q10 * 0.65 / (2.5 + exp((arg_v + 82.0) / 17.0));
        } else if (i == 1.0) {
            ret_bet = q10 * 0.001 * exp((arg_v - 168.0) / 16.0);
        }
        return ret_bet;
    }


    inline double ce_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double i) {
        double ret_ce = 0.0;
        arg_v = arg_v;
        if (i == 0.0) {
            ret_ce = 1.0 / (1.0 + exp( -(arg_v + 30.3) / 9.6));
        } else if (i == 1.0) {
            ret_ce = 1.0 * (0.25 + 1.0 / (1.35 + exp((arg_v + 7.0) / 14.0)));
        } else if (i == 2.0) {
            ret_ce = 1.0 * (0.1 + 1.0 / (1.1 + exp((arg_v + 7.0) / 14.0)));
        }
        return ret_ce;
    }


    namespace {
        struct _newton_states_Kv1_5 {
            int operator()(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) const {
                auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);
                double* savstate1 = static_cast<double*>(thread[dith1()].pval);
                auto const& slist1 = inst->global->slist1;
                auto const& dlist1 = inst->global->dlist1;
                double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (3*pnodecount);
                {
                    double a, b, c, alp_in_0, bet_in_0, ce_in_0, alp_in_1, bet_in_1, ce_in_1, ce_in_2, v_in_8;
                    v_in_8 = v;
                    {
                        double q10, v_in_0, i_in_0;
                        v_in_0 = v_in_8;
                        i_in_0 = 0.0;
                        v_in_0 = v_in_0;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_0 == 0.0) {
                            alp_in_0 = q10 * 0.65 / (exp( -(v_in_0 + 10.0) / 8.5) + exp( -(v_in_0 - 30.0) / 59.0));
                        } else if (i_in_0 == 1.0) {
                            alp_in_0 = 0.001 * q10 / (2.4 + 10.9 * exp( -(v_in_0 + 90.0) / 78.0));
                        }
                    }
                    a = alp_in_0;
                    {
                        double q10, v_in_1, i_in_1;
                        v_in_1 = v_in_8;
                        i_in_1 = 0.0;
                        v_in_1 = v_in_1;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_1 == 0.0) {
                            bet_in_0 = q10 * 0.65 / (2.5 + exp((v_in_1 + 82.0) / 17.0));
                        } else if (i_in_1 == 1.0) {
                            bet_in_0 = q10 * 0.001 * exp((v_in_1 - 168.0) / 16.0);
                        }
                    }
                    b = bet_in_0;
                    {
                        double v_in_2, i_in_2;
                        v_in_2 = v_in_8;
                        i_in_2 = 0.0;
                        v_in_2 = v_in_2;
                        if (i_in_2 == 0.0) {
                            ce_in_0 = 1.0 / (1.0 + exp( -(v_in_2 + 30.3) / 9.6));
                        } else if (i_in_2 == 1.0) {
                            ce_in_0 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_2 + 7.0) / 14.0)));
                        } else if (i_in_2 == 2.0) {
                            ce_in_0 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_2 + 7.0) / 14.0)));
                        }
                    }
                    c = ce_in_0;
                    inst->mtau[id] = 1.0 / (a + b) / 3.0 * inst->Tauact[id];
                    inst->minf[id] = c;
                    {
                        double q10, v_in_3, i_in_3;
                        v_in_3 = v_in_8;
                        i_in_3 = 1.0;
                        v_in_3 = v_in_3;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_3 == 0.0) {
                            alp_in_1 = q10 * 0.65 / (exp( -(v_in_3 + 10.0) / 8.5) + exp( -(v_in_3 - 30.0) / 59.0));
                        } else if (i_in_3 == 1.0) {
                            alp_in_1 = 0.001 * q10 / (2.4 + 10.9 * exp( -(v_in_3 + 90.0) / 78.0));
                        }
                    }
                    a = alp_in_1;
                    {
                        double q10, v_in_4, i_in_4;
                        v_in_4 = v_in_8;
                        i_in_4 = 1.0;
                        v_in_4 = v_in_4;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_4 == 0.0) {
                            bet_in_1 = q10 * 0.65 / (2.5 + exp((v_in_4 + 82.0) / 17.0));
                        } else if (i_in_4 == 1.0) {
                            bet_in_1 = q10 * 0.001 * exp((v_in_4 - 168.0) / 16.0);
                        }
                    }
                    b = bet_in_1;
                    {
                        double v_in_5, i_in_5;
                        v_in_5 = v_in_8;
                        i_in_5 = 1.0;
                        v_in_5 = v_in_5;
                        if (i_in_5 == 0.0) {
                            ce_in_1 = 1.0 / (1.0 + exp( -(v_in_5 + 30.3) / 9.6));
                        } else if (i_in_5 == 1.0) {
                            ce_in_1 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_5 + 7.0) / 14.0)));
                        } else if (i_in_5 == 2.0) {
                            ce_in_1 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_5 + 7.0) / 14.0)));
                        }
                    }
                    c = ce_in_1;
                    inst->ntau[id] = 1.0 / (a + b) / 3.0 * inst->Tauinactf[id];
                    inst->ninf[id] = c;
                    {
                        double v_in_6, i_in_6;
                        v_in_6 = v_in_8;
                        i_in_6 = 2.0;
                        v_in_6 = v_in_6;
                        if (i_in_6 == 0.0) {
                            ce_in_2 = 1.0 / (1.0 + exp( -(v_in_6 + 30.3) / 9.6));
                        } else if (i_in_6 == 1.0) {
                            ce_in_2 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_6 + 7.0) / 14.0)));
                        } else if (i_in_6 == 2.0) {
                            ce_in_2 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_6 + 7.0) / 14.0)));
                        }
                    }
                    c = ce_in_2;
                    inst->uinf[id] = c;
                    inst->utau[id] = 6800.0 * inst->Tauinacts[id];
                }
                inst->Dm[id] = (inst->minf[id] - inst->m[id]) / inst->mtau[id];
                inst->Dn[id] = (inst->ninf[id] - inst->n[id]) / inst->ntau[id];
                inst->Du[id] = (inst->uinf[id] - inst->u[id]) / inst->utau[id];
                int counter = -1;
                for (int i=0; i<3; i++) {
                    if (*deriv1_advance(thread)) {
                        dlist2[(++counter)*pnodecount+id] = data[dlist1[i]*pnodecount+id]-(data[slist1[i]*pnodecount+id]-savstate1[i*pnodecount+id])/nt->_dt;
                    } else {
                        dlist2[(++counter)*pnodecount+id] = data[slist1[i]*pnodecount+id]-savstate1[i*pnodecount+id];
                    }
                }
                return 0;
            }
        };
    }

    int states_Kv1_5(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);
        double* savstate1 = (double*) thread[dith1()].pval;
        auto const& slist1 = inst->global->slist1;
        auto& slist2 = inst->global->slist2;
        double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (3*pnodecount);
        for (int i=0; i<3; i++) {
            savstate1[i*pnodecount+id] = data[slist1[i]*pnodecount+id];
        }
        int reset = nrn_newton_thread(static_cast<NewtonSpace*>(*newtonspace1(thread)), 3, slist2, _newton_states_Kv1_5{}, dlist2, id, pnodecount, data, indexes, thread, nt, ml, v);
        return reset;
    }




    /** initialize channel */
    void nrn_init_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);


        int& deriv_advance_flag = *deriv1_advance(thread);
        deriv_advance_flag = 0;
        auto ns = newtonspace1(thread);
        auto& th = thread[dith1()];
        if (*ns == nullptr) {
            int vec_size = 2*3*pnodecount*sizeof(double);
            double* vec = makevector(vec_size);
            th.pval = vec;
            *ns = nrn_cons_newtonspace(3, pnodecount);
        }
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
                inst->ki[id] = inst->ion_ki[indexes[1*pnodecount + id]];
                inst->ko[id] = inst->ion_ko[indexes[2*pnodecount + id]];
                inst->nai[id] = inst->ion_nai[indexes[5*pnodecount + id]];
                inst->nao[id] = inst->ion_nao[indexes[6*pnodecount + id]];
                inst->m[id] = inst->global->m0;
                inst->n[id] = inst->global->n0;
                inst->u[id] = inst->global->u0;
                {
                    double a, b, c, alp_in_0, bet_in_0, ce_in_0, alp_in_1, bet_in_1, ce_in_1, ce_in_2, v_in_7;
                    v_in_7 = v;
                    {
                        double q10, v_in_0, i_in_0;
                        v_in_0 = v_in_7;
                        i_in_0 = 0.0;
                        v_in_0 = v_in_0;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_0 == 0.0) {
                            alp_in_0 = q10 * 0.65 / (exp( -(v_in_0 + 10.0) / 8.5) + exp( -(v_in_0 - 30.0) / 59.0));
                        } else if (i_in_0 == 1.0) {
                            alp_in_0 = 0.001 * q10 / (2.4 + 10.9 * exp( -(v_in_0 + 90.0) / 78.0));
                        }
                    }
                    a = alp_in_0;
                    {
                        double q10, v_in_1, i_in_1;
                        v_in_1 = v_in_7;
                        i_in_1 = 0.0;
                        v_in_1 = v_in_1;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_1 == 0.0) {
                            bet_in_0 = q10 * 0.65 / (2.5 + exp((v_in_1 + 82.0) / 17.0));
                        } else if (i_in_1 == 1.0) {
                            bet_in_0 = q10 * 0.001 * exp((v_in_1 - 168.0) / 16.0);
                        }
                    }
                    b = bet_in_0;
                    {
                        double v_in_2, i_in_2;
                        v_in_2 = v_in_7;
                        i_in_2 = 0.0;
                        v_in_2 = v_in_2;
                        if (i_in_2 == 0.0) {
                            ce_in_0 = 1.0 / (1.0 + exp( -(v_in_2 + 30.3) / 9.6));
                        } else if (i_in_2 == 1.0) {
                            ce_in_0 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_2 + 7.0) / 14.0)));
                        } else if (i_in_2 == 2.0) {
                            ce_in_0 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_2 + 7.0) / 14.0)));
                        }
                    }
                    c = ce_in_0;
                    inst->mtau[id] = 1.0 / (a + b) / 3.0 * inst->Tauact[id];
                    inst->minf[id] = c;
                    {
                        double q10, v_in_3, i_in_3;
                        v_in_3 = v_in_7;
                        i_in_3 = 1.0;
                        v_in_3 = v_in_3;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_3 == 0.0) {
                            alp_in_1 = q10 * 0.65 / (exp( -(v_in_3 + 10.0) / 8.5) + exp( -(v_in_3 - 30.0) / 59.0));
                        } else if (i_in_3 == 1.0) {
                            alp_in_1 = 0.001 * q10 / (2.4 + 10.9 * exp( -(v_in_3 + 90.0) / 78.0));
                        }
                    }
                    a = alp_in_1;
                    {
                        double q10, v_in_4, i_in_4;
                        v_in_4 = v_in_7;
                        i_in_4 = 1.0;
                        v_in_4 = v_in_4;
                        q10 = pow(2.2, ((*(inst->celsius) - 37.0) / 10.0));
                        if (i_in_4 == 0.0) {
                            bet_in_1 = q10 * 0.65 / (2.5 + exp((v_in_4 + 82.0) / 17.0));
                        } else if (i_in_4 == 1.0) {
                            bet_in_1 = q10 * 0.001 * exp((v_in_4 - 168.0) / 16.0);
                        }
                    }
                    b = bet_in_1;
                    {
                        double v_in_5, i_in_5;
                        v_in_5 = v_in_7;
                        i_in_5 = 1.0;
                        v_in_5 = v_in_5;
                        if (i_in_5 == 0.0) {
                            ce_in_1 = 1.0 / (1.0 + exp( -(v_in_5 + 30.3) / 9.6));
                        } else if (i_in_5 == 1.0) {
                            ce_in_1 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_5 + 7.0) / 14.0)));
                        } else if (i_in_5 == 2.0) {
                            ce_in_1 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_5 + 7.0) / 14.0)));
                        }
                    }
                    c = ce_in_1;
                    inst->ntau[id] = 1.0 / (a + b) / 3.0 * inst->Tauinactf[id];
                    inst->ninf[id] = c;
                    {
                        double v_in_6, i_in_6;
                        v_in_6 = v_in_7;
                        i_in_6 = 2.0;
                        v_in_6 = v_in_6;
                        if (i_in_6 == 0.0) {
                            ce_in_2 = 1.0 / (1.0 + exp( -(v_in_6 + 30.3) / 9.6));
                        } else if (i_in_6 == 1.0) {
                            ce_in_2 = 1.0 * (0.25 + 1.0 / (1.35 + exp((v_in_6 + 7.0) / 14.0)));
                        } else if (i_in_6 == 2.0) {
                            ce_in_2 = 1.0 * (0.1 + 1.0 / (1.1 + exp((v_in_6 + 7.0) / 14.0)));
                        }
                    }
                    c = ce_in_2;
                    inst->uinf[id] = c;
                    inst->utau[id] = 6800.0 * inst->Tauinacts[id];
                }
                inst->m[id] = inst->minf[id];
                inst->n[id] = inst->ninf[id];
                inst->u[id] = inst->uinf[id];
            }
        }
        deriv_advance_flag = 1;
    }


    inline double nrn_current_Kv1_5(int id, int pnodecount, Kv1_5_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double z;
        z = (R * (*(inst->celsius) + 273.15)) / F;
        inst->ik[id] = inst->gKur[id] * (0.1 + 1.0 / (1.0 + exp( -(v - 15.0) / 13.0))) * inst->m[id] * inst->m[id] * inst->m[id] * inst->n[id] * inst->u[id] * (v - inst->ek[id]);
        inst->ino[id] = inst->gnonspec[id] * (0.1 + 1.0 / (1.0 + exp( -(v - 15.0) / 13.0))) * inst->m[id] * inst->m[id] * inst->m[id] * inst->n[id] * inst->u[id] * (v - z * log((inst->nao[id] + inst->ko[id]) / (inst->nai[id] + inst->ki[id])));
        current += inst->ik[id];
        current += inst->ino[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            inst->ki[id] = inst->ion_ki[indexes[1*pnodecount + id]];
            inst->ko[id] = inst->ion_ko[indexes[2*pnodecount + id]];
            inst->nai[id] = inst->ion_nai[indexes[5*pnodecount + id]];
            inst->nao[id] = inst->ion_nao[indexes[6*pnodecount + id]];
            double g = nrn_current_Kv1_5(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double dino = inst->ino[id];
            double rhs = nrn_current_Kv1_5(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            inst->ion_dikdv[indexes[4*pnodecount + id]] += (dik-inst->ik[id])/0.001;
            inst->ion_dinodv[indexes[8*pnodecount + id]] += (dino-inst->ino[id])/0.001;
            inst->ion_ik[indexes[3*pnodecount + id]] += inst->ik[id];
            inst->ion_ino[indexes[7*pnodecount + id]] += inst->ino[id];
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            vec_rhs[node_id] -= rhs;
            vec_d[node_id] += g;
        }
    }


    /** update state */
    void nrn_state_Kv1_5(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv1_5_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            inst->ki[id] = inst->ion_ki[indexes[1*pnodecount + id]];
            inst->ko[id] = inst->ion_ko[indexes[2*pnodecount + id]];
            inst->nai[id] = inst->ion_nai[indexes[5*pnodecount + id]];
            inst->nao[id] = inst->ion_nao[indexes[6*pnodecount + id]];
            states_Kv1_5(id, pnodecount, data, indexes, thread, nt, ml, v);
        }
    }


    /** register channel with the simulator */
    void _Kv15_reg() {

        int mech_type = nrn_get_mechtype("Kv1_5");
        Kv1_5_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kv1_5, nrn_cur_Kv1_5, nullptr, nrn_state_Kv1_5, nrn_init_Kv1_5, nrn_private_constructor_Kv1_5, nrn_private_destructor_Kv1_5, first_pointer_var_index(), 4);
        Kv1_5_global.k_type = nrn_get_mechtype("k_ion");
        Kv1_5_global.na_type = nrn_get_mechtype("na_ion");
        Kv1_5_global.no_type = nrn_get_mechtype("no_ion");

        thread_mem_init(Kv1_5_global.ext_call_thread);
        _nrn_thread_reg0(mech_type, thread_mem_cleanup);
        _nrn_thread_reg1(mech_type, thread_mem_init);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_dparam_semantics(mech_type, 3, "k_ion");
        hoc_register_dparam_semantics(mech_type, 4, "k_ion");
        hoc_register_dparam_semantics(mech_type, 5, "na_ion");
        hoc_register_dparam_semantics(mech_type, 6, "na_ion");
        hoc_register_dparam_semantics(mech_type, 7, "no_ion");
        hoc_register_dparam_semantics(mech_type, 8, "no_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
