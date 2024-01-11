/*********************************************************
Model Name      : Kv4_3
Filename        : Kv43.mod
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
        "Kv4_3",
        "Aalpha_a_Kv4_3",
        "Kalpha_a_Kv4_3",
        "V0alpha_a_Kv4_3",
        "Abeta_a_Kv4_3",
        "Kbeta_a_Kv4_3",
        "V0beta_a_Kv4_3",
        "Aalpha_b_Kv4_3",
        "Kalpha_b_Kv4_3",
        "V0alpha_b_Kv4_3",
        "Abeta_b_Kv4_3",
        "Kbeta_b_Kv4_3",
        "V0beta_b_Kv4_3",
        "V0_ainf_Kv4_3",
        "K_ainf_Kv4_3",
        "V0_binf_Kv4_3",
        "K_binf_Kv4_3",
        "gkbar_Kv4_3",
        0,
        "ik_Kv4_3",
        "a_inf_Kv4_3",
        "b_inf_Kv4_3",
        "tau_a_Kv4_3",
        "tau_b_Kv4_3",
        "g_Kv4_3",
        "alpha_a_Kv4_3",
        "beta_a_Kv4_3",
        "alpha_b_Kv4_3",
        "beta_b_Kv4_3",
        0,
        "a_Kv4_3",
        "b_Kv4_3",
        0,
        0
    };


    /** all global variables */
    struct Kv4_3_Store {
        int k_type{};
        double a0{};
        double b0{};
        int reset{};
        int mech_type{};
        int slist1[2]{27, 28};
        int dlist1[2]{29, 30};
        int slist2[2]{27, 28};
        double usetable{1};
        double tmin_rate{};
        double mfac_rate{};
        double t_a_inf[13001]{};
        double t_tau_a[13001]{};
        double t_b_inf[13001]{};
        double t_tau_b[13001]{};
        ThreadDatum ext_call_thread[3]{};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kv4_3_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kv4_3_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kv4_3_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kv4_3_Store>);
    static_assert(std::is_trivially_destructible_v<Kv4_3_Store>);
    Kv4_3_Store Kv4_3_global;


    /** all mechanism instance variables and global variables */
    struct Kv4_3_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* Aalpha_a{};
        const double* Kalpha_a{};
        const double* V0alpha_a{};
        const double* Abeta_a{};
        const double* Kbeta_a{};
        const double* V0beta_a{};
        const double* Aalpha_b{};
        const double* Kalpha_b{};
        const double* V0alpha_b{};
        const double* Abeta_b{};
        const double* Kbeta_b{};
        const double* V0beta_b{};
        const double* V0_ainf{};
        const double* K_ainf{};
        const double* V0_binf{};
        const double* K_binf{};
        const double* gkbar{};
        double* ik{};
        double* a_inf{};
        double* b_inf{};
        double* tau_a{};
        double* tau_b{};
        double* g{};
        double* alpha_a{};
        double* beta_a{};
        double* alpha_b{};
        double* beta_b{};
        double* a{};
        double* b{};
        double* Da{};
        double* Db{};
        double* ek{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Kv4_3_Store* global{&Kv4_3_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"usetable_Kv4_3", &Kv4_3_global.usetable},
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
        return 34;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return Kv4_3_global.mech_type;
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
    static void nrn_private_constructor_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kv4_3_Instance{};
        assert(inst->global == &Kv4_3_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kv4_3_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv4_3_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv4_3_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kv4_3_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kv4_3_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->Aalpha_a = ml->data+0*pnodecount;
        inst->Kalpha_a = ml->data+1*pnodecount;
        inst->V0alpha_a = ml->data+2*pnodecount;
        inst->Abeta_a = ml->data+3*pnodecount;
        inst->Kbeta_a = ml->data+4*pnodecount;
        inst->V0beta_a = ml->data+5*pnodecount;
        inst->Aalpha_b = ml->data+6*pnodecount;
        inst->Kalpha_b = ml->data+7*pnodecount;
        inst->V0alpha_b = ml->data+8*pnodecount;
        inst->Abeta_b = ml->data+9*pnodecount;
        inst->Kbeta_b = ml->data+10*pnodecount;
        inst->V0beta_b = ml->data+11*pnodecount;
        inst->V0_ainf = ml->data+12*pnodecount;
        inst->K_ainf = ml->data+13*pnodecount;
        inst->V0_binf = ml->data+14*pnodecount;
        inst->K_binf = ml->data+15*pnodecount;
        inst->gkbar = ml->data+16*pnodecount;
        inst->ik = ml->data+17*pnodecount;
        inst->a_inf = ml->data+18*pnodecount;
        inst->b_inf = ml->data+19*pnodecount;
        inst->tau_a = ml->data+20*pnodecount;
        inst->tau_b = ml->data+21*pnodecount;
        inst->g = ml->data+22*pnodecount;
        inst->alpha_a = ml->data+23*pnodecount;
        inst->beta_a = ml->data+24*pnodecount;
        inst->alpha_b = ml->data+25*pnodecount;
        inst->beta_b = ml->data+26*pnodecount;
        inst->a = ml->data+27*pnodecount;
        inst->b = ml->data+28*pnodecount;
        inst->Da = ml->data+29*pnodecount;
        inst->Db = ml->data+30*pnodecount;
        inst->ek = ml->data+31*pnodecount;
        inst->v_unused = ml->data+32*pnodecount;
        inst->g_unused = ml->data+33*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Kv4_3(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);

        #endif
    }


    inline double alp_a_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double bet_a_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double alp_b_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double bet_b_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double linoid_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double x, double y);
    inline double sigm_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double x, double y);
    inline int rate_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int f_rate_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_f_rate = 0;
        double a_a, b_a, a_b, b_b, alp_a_in_0, bet_a_in_0, alp_b_in_0, bet_b_in_0;
        {
            double Q10, sigm_in_0, v_in_0;
            v_in_0 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            {
                double x_in_0, y_in_0;
                x_in_0 = v_in_0 - inst->V0alpha_a[id];
                y_in_0 = inst->Kalpha_a[id];
                sigm_in_0 = 1.0 / (exp(x_in_0 / y_in_0) + 1.0);
            }
            alp_a_in_0 = Q10 * inst->Aalpha_a[id] * sigm_in_0;
        }
        a_a = alp_a_in_0;
        {
            double Q10, v_in_1;
            v_in_1 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            bet_a_in_0 = Q10 * inst->Abeta_a[id] / (exp((v_in_1 - inst->V0beta_a[id]) / inst->Kbeta_a[id]));
        }
        b_a = bet_a_in_0;
        {
            double Q10, sigm_in_1, v_in_2;
            v_in_2 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            {
                double x_in_1, y_in_1;
                x_in_1 = v_in_2 - inst->V0alpha_b[id];
                y_in_1 = inst->Kalpha_b[id];
                sigm_in_1 = 1.0 / (exp(x_in_1 / y_in_1) + 1.0);
            }
            alp_b_in_0 = Q10 * inst->Aalpha_b[id] * sigm_in_1;
        }
        a_b = alp_b_in_0;
        {
            double Q10, sigm_in_2, v_in_3;
            v_in_3 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            {
                double x_in_2, y_in_2;
                x_in_2 = v_in_3 - inst->V0beta_b[id];
                y_in_2 = inst->Kbeta_b[id];
                sigm_in_2 = 1.0 / (exp(x_in_2 / y_in_2) + 1.0);
            }
            bet_b_in_0 = Q10 * inst->Abeta_b[id] * sigm_in_2;
        }
        b_b = bet_b_in_0;
        inst->a_inf[id] = 1.0 / (1.0 + exp((arg_v - inst->V0_ainf[id]) / inst->K_ainf[id]));
        inst->tau_a[id] = 1.0 / (a_a + b_a);
        inst->b_inf[id] = 1.0 / (1.0 + exp((arg_v - inst->V0_binf[id]) / inst->K_binf[id]));
        inst->tau_b[id] = 1.0 / (a_b + b_b);
        return ret_f_rate;
    }


    void check_rate_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        if (inst->global->usetable == 0) {
            return;
        }
        static bool make_table = true;
        static double save_Aalpha_a;
        static double save_Kalpha_a;
        static double save_V0alpha_a;
        static double save_Abeta_a;
        static double save_Kbeta_a;
        static double save_V0beta_a;
        static double save_Aalpha_b;
        static double save_Kalpha_b;
        static double save_V0alpha_b;
        static double save_Abeta_b;
        static double save_Kbeta_b;
        static double save_V0beta_b;
        static double save_celsius;
        if (save_Aalpha_a != inst->Aalpha_a[id]) {
            make_table = true;
        }
        if (save_Kalpha_a != inst->Kalpha_a[id]) {
            make_table = true;
        }
        if (save_V0alpha_a != inst->V0alpha_a[id]) {
            make_table = true;
        }
        if (save_Abeta_a != inst->Abeta_a[id]) {
            make_table = true;
        }
        if (save_Kbeta_a != inst->Kbeta_a[id]) {
            make_table = true;
        }
        if (save_V0beta_a != inst->V0beta_a[id]) {
            make_table = true;
        }
        if (save_Aalpha_b != inst->Aalpha_b[id]) {
            make_table = true;
        }
        if (save_Kalpha_b != inst->Kalpha_b[id]) {
            make_table = true;
        }
        if (save_V0alpha_b != inst->V0alpha_b[id]) {
            make_table = true;
        }
        if (save_Abeta_b != inst->Abeta_b[id]) {
            make_table = true;
        }
        if (save_Kbeta_b != inst->Kbeta_b[id]) {
            make_table = true;
        }
        if (save_V0beta_b != inst->V0beta_b[id]) {
            make_table = true;
        }
        if (save_celsius != *(inst->celsius)) {
            make_table = true;
        }
        if (make_table) {
            make_table = false;
            inst->global->tmin_rate =  -100.0;
            double tmax = 30.0;
            double dx = (tmax-inst->global->tmin_rate) / 13000.;
            inst->global->mfac_rate = 1./dx;
            double x = inst->global->tmin_rate;
            for (std::size_t i = 0; i < 13001; x += dx, i++) {
                f_rate_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v, x);
                inst->global->t_a_inf[i] = inst->a_inf[id];
                inst->global->t_tau_a[i] = inst->tau_a[id];
                inst->global->t_b_inf[i] = inst->b_inf[id];
                inst->global->t_tau_b[i] = inst->tau_b[id];
            }
            save_Aalpha_a = inst->Aalpha_a[id];
            save_Kalpha_a = inst->Kalpha_a[id];
            save_V0alpha_a = inst->V0alpha_a[id];
            save_Abeta_a = inst->Abeta_a[id];
            save_Kbeta_a = inst->Kbeta_a[id];
            save_V0beta_a = inst->V0beta_a[id];
            save_Aalpha_b = inst->Aalpha_b[id];
            save_Kalpha_b = inst->Kalpha_b[id];
            save_V0alpha_b = inst->V0alpha_b[id];
            save_Abeta_b = inst->Abeta_b[id];
            save_Kbeta_b = inst->Kbeta_b[id];
            save_V0beta_b = inst->V0beta_b[id];
            save_celsius = *(inst->celsius);
        }
    }


    inline int rate_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v){
        if (inst->global->usetable == 0) {
            f_rate_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v, arg_v);
            return 0;
        }
        double xi = inst->global->mfac_rate * (arg_v - inst->global->tmin_rate);
        if (isnan(xi)) {
            inst->a_inf[id] = xi;
            inst->tau_a[id] = xi;
            inst->b_inf[id] = xi;
            inst->tau_b[id] = xi;
            return 0;
        }
        if (xi <= 0. || xi >= 13000.) {
            int index = (xi <= 0.) ? 0 : 13000;
            inst->a_inf[id] = inst->global->t_a_inf[index];
            inst->tau_a[id] = inst->global->t_tau_a[index];
            inst->b_inf[id] = inst->global->t_b_inf[index];
            inst->tau_b[id] = inst->global->t_tau_b[index];
            return 0;
        }
        int i = int(xi);
        double theta = xi - double(i);
        inst->a_inf[id] = inst->global->t_a_inf[i] + theta*(inst->global->t_a_inf[i+1]-inst->global->t_a_inf[i]);
        inst->tau_a[id] = inst->global->t_tau_a[i] + theta*(inst->global->t_tau_a[i+1]-inst->global->t_tau_a[i]);
        inst->b_inf[id] = inst->global->t_b_inf[i] + theta*(inst->global->t_b_inf[i+1]-inst->global->t_b_inf[i]);
        inst->tau_b[id] = inst->global->t_tau_b[i] + theta*(inst->global->t_tau_b[i+1]-inst->global->t_tau_b[i]);
        return 0;
    }


    inline double alp_a_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alp_a = 0.0;
        double Q10, sigm_in_0;
        Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
        {
            double x_in_0, y_in_0;
            x_in_0 = arg_v - inst->V0alpha_a[id];
            y_in_0 = inst->Kalpha_a[id];
            sigm_in_0 = 1.0 / (exp(x_in_0 / y_in_0) + 1.0);
        }
        ret_alp_a = Q10 * inst->Aalpha_a[id] * sigm_in_0;
        return ret_alp_a;
    }


    inline double bet_a_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_bet_a = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
        ret_bet_a = Q10 * inst->Abeta_a[id] / (exp((arg_v - inst->V0beta_a[id]) / inst->Kbeta_a[id]));
        return ret_bet_a;
    }


    inline double alp_b_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alp_b = 0.0;
        double Q10, sigm_in_1;
        Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
        {
            double x_in_1, y_in_1;
            x_in_1 = arg_v - inst->V0alpha_b[id];
            y_in_1 = inst->Kalpha_b[id];
            sigm_in_1 = 1.0 / (exp(x_in_1 / y_in_1) + 1.0);
        }
        ret_alp_b = Q10 * inst->Aalpha_b[id] * sigm_in_1;
        return ret_alp_b;
    }


    inline double bet_b_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_bet_b = 0.0;
        double Q10, sigm_in_2;
        Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
        {
            double x_in_2, y_in_2;
            x_in_2 = arg_v - inst->V0beta_b[id];
            y_in_2 = inst->Kbeta_b[id];
            sigm_in_2 = 1.0 / (exp(x_in_2 / y_in_2) + 1.0);
        }
        ret_bet_b = Q10 * inst->Abeta_b[id] * sigm_in_2;
        return ret_bet_b;
    }


    inline double linoid_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double x, double y) {
        double ret_linoid = 0.0;
        if (fabs(x / y) < 1e-6) {
            ret_linoid = y * (1.0 - x / y / 2.0);
        } else {
            ret_linoid = x / (exp(x / y) - 1.0);
        }
        return ret_linoid;
    }


    inline double sigm_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double x, double y) {
        double ret_sigm = 0.0;
        ret_sigm = 1.0 / (exp(x / y) + 1.0);
        return ret_sigm;
    }


    namespace {
        struct _newton_states_Kv4_3 {
            int operator()(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) const {
                auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);
                double* savstate1 = static_cast<double*>(thread[dith1()].pval);
                auto const& slist1 = inst->global->slist1;
                auto const& dlist1 = inst->global->dlist1;
                double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (2*pnodecount);
                rate_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->Da[id] = (inst->a_inf[id] - inst->a[id]) / inst->tau_a[id];
                inst->Db[id] = (inst->b_inf[id] - inst->b[id]) / inst->tau_b[id];
                int counter = -1;
                for (int i=0; i<2; i++) {
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

    int states_Kv4_3(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);
        double* savstate1 = (double*) thread[dith1()].pval;
        auto const& slist1 = inst->global->slist1;
        auto& slist2 = inst->global->slist2;
        double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (2*pnodecount);
        for (int i=0; i<2; i++) {
            savstate1[i*pnodecount+id] = data[slist1[i]*pnodecount+id];
        }
        int reset = nrn_newton_thread(static_cast<NewtonSpace*>(*newtonspace1(thread)), 2, slist2, _newton_states_Kv4_3{}, dlist2, id, pnodecount, data, indexes, thread, nt, ml, v);
        return reset;
    }




    /** initialize channel */
    void nrn_init_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);


        int& deriv_advance_flag = *deriv1_advance(thread);
        deriv_advance_flag = 0;
        auto ns = newtonspace1(thread);
        auto& th = thread[dith1()];
        if (*ns == nullptr) {
            int vec_size = 2*2*pnodecount*sizeof(double);
            double* vec = makevector(vec_size);
            th.pval = vec;
            *ns = nrn_cons_newtonspace(2, pnodecount);
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
                inst->a[id] = inst->global->a0;
                inst->b[id] = inst->global->b0;
                rate_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->a[id] = inst->a_inf[id];
                inst->b[id] = inst->b_inf[id];
            }
        }
        deriv_advance_flag = 1;
    }


    inline double nrn_current_Kv4_3(int id, int pnodecount, Kv4_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double alp_a_in_1, bet_a_in_1, alp_b_in_1, bet_b_in_1;
        inst->g[id] = inst->gkbar[id] * inst->a[id] * inst->a[id] * inst->a[id] * inst->b[id];
        inst->ik[id] = inst->g[id] * (v - inst->ek[id]);
        {
            double Q10, sigm_in_0, v_in_4;
            v_in_4 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            {
                double x_in_0, y_in_0;
                x_in_0 = v_in_4 - inst->V0alpha_a[id];
                y_in_0 = inst->Kalpha_a[id];
                sigm_in_0 = 1.0 / (exp(x_in_0 / y_in_0) + 1.0);
            }
            alp_a_in_1 = Q10 * inst->Aalpha_a[id] * sigm_in_0;
        }
        inst->alpha_a[id] = alp_a_in_1;
        {
            double Q10, v_in_5;
            v_in_5 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            bet_a_in_1 = Q10 * inst->Abeta_a[id] / (exp((v_in_5 - inst->V0beta_a[id]) / inst->Kbeta_a[id]));
        }
        inst->beta_a[id] = bet_a_in_1;
        {
            double Q10, sigm_in_1, v_in_6;
            v_in_6 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            {
                double x_in_1, y_in_1;
                x_in_1 = v_in_6 - inst->V0alpha_b[id];
                y_in_1 = inst->Kalpha_b[id];
                sigm_in_1 = 1.0 / (exp(x_in_1 / y_in_1) + 1.0);
            }
            alp_b_in_1 = Q10 * inst->Aalpha_b[id] * sigm_in_1;
        }
        inst->alpha_b[id] = alp_b_in_1;
        {
            double Q10, sigm_in_2, v_in_7;
            v_in_7 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 25.5) / 10.0));
            {
                double x_in_2, y_in_2;
                x_in_2 = v_in_7 - inst->V0beta_b[id];
                y_in_2 = inst->Kbeta_b[id];
                sigm_in_2 = 1.0 / (exp(x_in_2 / y_in_2) + 1.0);
            }
            bet_b_in_1 = Q10 * inst->Abeta_b[id] * sigm_in_2;
        }
        inst->beta_b[id] = bet_b_in_1;
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Kv4_3(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            states_Kv4_3(id, pnodecount, data, indexes, thread, nt, ml, v);
        }
    }


    static void check_table_thread_Kv4_3 (int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, int tml_id) {
        setup_instance(nt, ml);
        auto* const inst = static_cast<Kv4_3_Instance*>(ml->instance);
        double v = 0;
        check_rate_Kv4_3(id, pnodecount, inst, data, indexes, thread, nt, v);
    }


    /** register channel with the simulator */
    void _Kv43_reg() {

        int mech_type = nrn_get_mechtype("Kv4_3");
        Kv4_3_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kv4_3, nrn_cur_Kv4_3, nullptr, nrn_state_Kv4_3, nrn_init_Kv4_3, nrn_private_constructor_Kv4_3, nrn_private_destructor_Kv4_3, first_pointer_var_index(), 4);
        Kv4_3_global.k_type = nrn_get_mechtype("k_ion");

        thread_mem_init(Kv4_3_global.ext_call_thread);
        _nrn_thread_reg0(mech_type, thread_mem_cleanup);
        _nrn_thread_reg1(mech_type, thread_mem_init);
        _nrn_thread_table_reg(mech_type, check_table_thread_Kv4_3);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
