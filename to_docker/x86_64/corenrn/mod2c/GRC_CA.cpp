/*********************************************************
Model Name      : GRC_CA
Filename        : GRC_CA.mod
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
        "GRC_CA",
        "Aalpha_s_GRC_CA",
        "Kalpha_s_GRC_CA",
        "V0alpha_s_GRC_CA",
        "Abeta_s_GRC_CA",
        "Kbeta_s_GRC_CA",
        "V0beta_s_GRC_CA",
        "Aalpha_u_GRC_CA",
        "Kalpha_u_GRC_CA",
        "V0alpha_u_GRC_CA",
        "Abeta_u_GRC_CA",
        "Kbeta_u_GRC_CA",
        "V0beta_u_GRC_CA",
        "gcabar_GRC_CA",
        0,
        "ica_GRC_CA",
        "s_inf_GRC_CA",
        "u_inf_GRC_CA",
        "tau_s_GRC_CA",
        "tau_u_GRC_CA",
        "g_GRC_CA",
        "alpha_s_GRC_CA",
        "beta_s_GRC_CA",
        "alpha_u_GRC_CA",
        "beta_u_GRC_CA",
        0,
        "s_GRC_CA",
        "u_GRC_CA",
        0,
        0
    };


    /** all global variables */
    struct GRC_CA_Store {
        int ca_type{};
        double s0{};
        double u0{};
        int reset{};
        int mech_type{};
        int slist1[2]{23, 24};
        int dlist1[2]{26, 27};
        int slist2[2]{23, 24};
        double usetable{1};
        double tmin_rate{};
        double mfac_rate{};
        double t_s_inf[13001]{};
        double t_tau_s[13001]{};
        double t_u_inf[13001]{};
        double t_tau_u[13001]{};
        ThreadDatum ext_call_thread[3]{};
    };
    static_assert(std::is_trivially_copy_constructible_v<GRC_CA_Store>);
    static_assert(std::is_trivially_move_constructible_v<GRC_CA_Store>);
    static_assert(std::is_trivially_copy_assignable_v<GRC_CA_Store>);
    static_assert(std::is_trivially_move_assignable_v<GRC_CA_Store>);
    static_assert(std::is_trivially_destructible_v<GRC_CA_Store>);
    GRC_CA_Store GRC_CA_global;


    /** all mechanism instance variables and global variables */
    struct GRC_CA_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* Aalpha_s{};
        const double* Kalpha_s{};
        const double* V0alpha_s{};
        const double* Abeta_s{};
        const double* Kbeta_s{};
        const double* V0beta_s{};
        const double* Aalpha_u{};
        const double* Kalpha_u{};
        const double* V0alpha_u{};
        const double* Abeta_u{};
        const double* Kbeta_u{};
        const double* V0beta_u{};
        const double* gcabar{};
        double* ica{};
        double* s_inf{};
        double* u_inf{};
        double* tau_s{};
        double* tau_u{};
        double* g{};
        double* alpha_s{};
        double* beta_s{};
        double* alpha_u{};
        double* beta_u{};
        double* s{};
        double* u{};
        double* eca{};
        double* Ds{};
        double* Du{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_eca{};
        double* ion_ica{};
        double* ion_dicadv{};
        GRC_CA_Store* global{&GRC_CA_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"usetable_GRC_CA", &GRC_CA_global.usetable},
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
        return 30;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return GRC_CA_global.mech_type;
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
    static void nrn_private_constructor_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new GRC_CA_Instance{};
        assert(inst->global == &GRC_CA_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(GRC_CA_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &GRC_CA_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(GRC_CA_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &GRC_CA_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(GRC_CA_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->Aalpha_s = ml->data+0*pnodecount;
        inst->Kalpha_s = ml->data+1*pnodecount;
        inst->V0alpha_s = ml->data+2*pnodecount;
        inst->Abeta_s = ml->data+3*pnodecount;
        inst->Kbeta_s = ml->data+4*pnodecount;
        inst->V0beta_s = ml->data+5*pnodecount;
        inst->Aalpha_u = ml->data+6*pnodecount;
        inst->Kalpha_u = ml->data+7*pnodecount;
        inst->V0alpha_u = ml->data+8*pnodecount;
        inst->Abeta_u = ml->data+9*pnodecount;
        inst->Kbeta_u = ml->data+10*pnodecount;
        inst->V0beta_u = ml->data+11*pnodecount;
        inst->gcabar = ml->data+12*pnodecount;
        inst->ica = ml->data+13*pnodecount;
        inst->s_inf = ml->data+14*pnodecount;
        inst->u_inf = ml->data+15*pnodecount;
        inst->tau_s = ml->data+16*pnodecount;
        inst->tau_u = ml->data+17*pnodecount;
        inst->g = ml->data+18*pnodecount;
        inst->alpha_s = ml->data+19*pnodecount;
        inst->beta_s = ml->data+20*pnodecount;
        inst->alpha_u = ml->data+21*pnodecount;
        inst->beta_u = ml->data+22*pnodecount;
        inst->s = ml->data+23*pnodecount;
        inst->u = ml->data+24*pnodecount;
        inst->eca = ml->data+25*pnodecount;
        inst->Ds = ml->data+26*pnodecount;
        inst->Du = ml->data+27*pnodecount;
        inst->v_unused = ml->data+28*pnodecount;
        inst->g_unused = ml->data+29*pnodecount;
        inst->ion_eca = nt->_data;
        inst->ion_ica = nt->_data;
        inst->ion_dicadv = nt->_data;
    }



    static void nrn_alloc_GRC_CA(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);

        #endif
    }


    inline double alp_s_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double bet_s_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double alp_u_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double bet_u_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline int rate_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int f_rate_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_f_rate = 0;
        double a_s, b_s, a_u, b_u, alp_s_in_0, bet_s_in_0, alp_u_in_0, bet_u_in_0;
        {
            double Q10, v_in_0;
            v_in_0 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            alp_s_in_0 = Q10 * inst->Aalpha_s[id] * exp((v_in_0 - inst->V0alpha_s[id]) / inst->Kalpha_s[id]);
        }
        a_s = alp_s_in_0;
        {
            double Q10, v_in_1;
            v_in_1 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            bet_s_in_0 = Q10 * inst->Abeta_s[id] * exp((v_in_1 - inst->V0beta_s[id]) / inst->Kbeta_s[id]);
        }
        b_s = bet_s_in_0;
        {
            double Q10, v_in_2;
            v_in_2 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            alp_u_in_0 = Q10 * inst->Aalpha_u[id] * exp((v_in_2 - inst->V0alpha_u[id]) / inst->Kalpha_u[id]);
        }
        a_u = alp_u_in_0;
        {
            double Q10, v_in_3;
            v_in_3 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            bet_u_in_0 = Q10 * inst->Abeta_u[id] * exp((v_in_3 - inst->V0beta_u[id]) / inst->Kbeta_u[id]);
        }
        b_u = bet_u_in_0;
        inst->s_inf[id] = a_s / (a_s + b_s);
        inst->tau_s[id] = 1.0 / (a_s + b_s);
        inst->u_inf[id] = a_u / (a_u + b_u);
        inst->tau_u[id] = 1.0 / (a_u + b_u);
        return ret_f_rate;
    }


    void check_rate_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        if (inst->global->usetable == 0) {
            return;
        }
        static bool make_table = true;
        static double save_Aalpha_s;
        static double save_Kalpha_s;
        static double save_V0alpha_s;
        static double save_Abeta_s;
        static double save_Kbeta_s;
        static double save_V0beta_s;
        static double save_Aalpha_u;
        static double save_Kalpha_u;
        static double save_V0alpha_u;
        static double save_Abeta_u;
        static double save_Kbeta_u;
        static double save_V0beta_u;
        static double save_celsius;
        if (save_Aalpha_s != inst->Aalpha_s[id]) {
            make_table = true;
        }
        if (save_Kalpha_s != inst->Kalpha_s[id]) {
            make_table = true;
        }
        if (save_V0alpha_s != inst->V0alpha_s[id]) {
            make_table = true;
        }
        if (save_Abeta_s != inst->Abeta_s[id]) {
            make_table = true;
        }
        if (save_Kbeta_s != inst->Kbeta_s[id]) {
            make_table = true;
        }
        if (save_V0beta_s != inst->V0beta_s[id]) {
            make_table = true;
        }
        if (save_Aalpha_u != inst->Aalpha_u[id]) {
            make_table = true;
        }
        if (save_Kalpha_u != inst->Kalpha_u[id]) {
            make_table = true;
        }
        if (save_V0alpha_u != inst->V0alpha_u[id]) {
            make_table = true;
        }
        if (save_Abeta_u != inst->Abeta_u[id]) {
            make_table = true;
        }
        if (save_Kbeta_u != inst->Kbeta_u[id]) {
            make_table = true;
        }
        if (save_V0beta_u != inst->V0beta_u[id]) {
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
                f_rate_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v, x);
                inst->global->t_s_inf[i] = inst->s_inf[id];
                inst->global->t_tau_s[i] = inst->tau_s[id];
                inst->global->t_u_inf[i] = inst->u_inf[id];
                inst->global->t_tau_u[i] = inst->tau_u[id];
            }
            save_Aalpha_s = inst->Aalpha_s[id];
            save_Kalpha_s = inst->Kalpha_s[id];
            save_V0alpha_s = inst->V0alpha_s[id];
            save_Abeta_s = inst->Abeta_s[id];
            save_Kbeta_s = inst->Kbeta_s[id];
            save_V0beta_s = inst->V0beta_s[id];
            save_Aalpha_u = inst->Aalpha_u[id];
            save_Kalpha_u = inst->Kalpha_u[id];
            save_V0alpha_u = inst->V0alpha_u[id];
            save_Abeta_u = inst->Abeta_u[id];
            save_Kbeta_u = inst->Kbeta_u[id];
            save_V0beta_u = inst->V0beta_u[id];
            save_celsius = *(inst->celsius);
        }
    }


    inline int rate_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v){
        if (inst->global->usetable == 0) {
            f_rate_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v, arg_v);
            return 0;
        }
        double xi = inst->global->mfac_rate * (arg_v - inst->global->tmin_rate);
        if (isnan(xi)) {
            inst->s_inf[id] = xi;
            inst->tau_s[id] = xi;
            inst->u_inf[id] = xi;
            inst->tau_u[id] = xi;
            return 0;
        }
        if (xi <= 0. || xi >= 13000.) {
            int index = (xi <= 0.) ? 0 : 13000;
            inst->s_inf[id] = inst->global->t_s_inf[index];
            inst->tau_s[id] = inst->global->t_tau_s[index];
            inst->u_inf[id] = inst->global->t_u_inf[index];
            inst->tau_u[id] = inst->global->t_tau_u[index];
            return 0;
        }
        int i = int(xi);
        double theta = xi - double(i);
        inst->s_inf[id] = inst->global->t_s_inf[i] + theta*(inst->global->t_s_inf[i+1]-inst->global->t_s_inf[i]);
        inst->tau_s[id] = inst->global->t_tau_s[i] + theta*(inst->global->t_tau_s[i+1]-inst->global->t_tau_s[i]);
        inst->u_inf[id] = inst->global->t_u_inf[i] + theta*(inst->global->t_u_inf[i+1]-inst->global->t_u_inf[i]);
        inst->tau_u[id] = inst->global->t_tau_u[i] + theta*(inst->global->t_tau_u[i+1]-inst->global->t_tau_u[i]);
        return 0;
    }


    inline double alp_s_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alp_s = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
        ret_alp_s = Q10 * inst->Aalpha_s[id] * exp((arg_v - inst->V0alpha_s[id]) / inst->Kalpha_s[id]);
        return ret_alp_s;
    }


    inline double bet_s_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_bet_s = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
        ret_bet_s = Q10 * inst->Abeta_s[id] * exp((arg_v - inst->V0beta_s[id]) / inst->Kbeta_s[id]);
        return ret_bet_s;
    }


    inline double alp_u_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alp_u = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
        ret_alp_u = Q10 * inst->Aalpha_u[id] * exp((arg_v - inst->V0alpha_u[id]) / inst->Kalpha_u[id]);
        return ret_alp_u;
    }


    inline double bet_u_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_bet_u = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
        ret_bet_u = Q10 * inst->Abeta_u[id] * exp((arg_v - inst->V0beta_u[id]) / inst->Kbeta_u[id]);
        return ret_bet_u;
    }


    namespace {
        struct _newton_states_GRC_CA {
            int operator()(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) const {
                auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);
                double* savstate1 = static_cast<double*>(thread[dith1()].pval);
                auto const& slist1 = inst->global->slist1;
                auto const& dlist1 = inst->global->dlist1;
                double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (2*pnodecount);
                rate_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->Ds[id] = (inst->s_inf[id] - inst->s[id]) / inst->tau_s[id];
                inst->Du[id] = (inst->u_inf[id] - inst->u[id]) / inst->tau_u[id];
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

    int states_GRC_CA(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);
        double* savstate1 = (double*) thread[dith1()].pval;
        auto const& slist1 = inst->global->slist1;
        auto& slist2 = inst->global->slist2;
        double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (2*pnodecount);
        for (int i=0; i<2; i++) {
            savstate1[i*pnodecount+id] = data[slist1[i]*pnodecount+id];
        }
        int reset = nrn_newton_thread(static_cast<NewtonSpace*>(*newtonspace1(thread)), 2, slist2, _newton_states_GRC_CA{}, dlist2, id, pnodecount, data, indexes, thread, nt, ml, v);
        return reset;
    }




    /** initialize channel */
    void nrn_init_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);


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
                inst->eca[id] = inst->ion_eca[indexes[0*pnodecount + id]];
                inst->s[id] = inst->global->s0;
                inst->u[id] = inst->global->u0;
                rate_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->s[id] = inst->s_inf[id];
                inst->u[id] = inst->u_inf[id];
            }
        }
        deriv_advance_flag = 1;
    }


    inline double nrn_current_GRC_CA(int id, int pnodecount, GRC_CA_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double alp_s_in_1, bet_s_in_1, alp_u_in_1, bet_u_in_1;
        inst->g[id] = inst->gcabar[id] * inst->s[id] * inst->s[id] * inst->u[id];
        inst->ica[id] = inst->g[id] * (v - inst->eca[id]);
        {
            double Q10, v_in_4;
            v_in_4 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            alp_s_in_1 = Q10 * inst->Aalpha_s[id] * exp((v_in_4 - inst->V0alpha_s[id]) / inst->Kalpha_s[id]);
        }
        inst->alpha_s[id] = alp_s_in_1;
        {
            double Q10, v_in_5;
            v_in_5 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            bet_s_in_1 = Q10 * inst->Abeta_s[id] * exp((v_in_5 - inst->V0beta_s[id]) / inst->Kbeta_s[id]);
        }
        inst->beta_s[id] = bet_s_in_1;
        {
            double Q10, v_in_6;
            v_in_6 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            alp_u_in_1 = Q10 * inst->Aalpha_u[id] * exp((v_in_6 - inst->V0alpha_u[id]) / inst->Kalpha_u[id]);
        }
        inst->alpha_u[id] = alp_u_in_1;
        {
            double Q10, v_in_7;
            v_in_7 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            bet_u_in_1 = Q10 * inst->Abeta_u[id] * exp((v_in_7 - inst->V0beta_u[id]) / inst->Kbeta_u[id]);
        }
        inst->beta_u[id] = bet_u_in_1;
        current += inst->ica[id];
        return current;
    }


    /** update current */
    void nrn_cur_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->eca[id] = inst->ion_eca[indexes[0*pnodecount + id]];
            double g = nrn_current_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dica = inst->ica[id];
            double rhs = nrn_current_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            inst->ion_dicadv[indexes[2*pnodecount + id]] += (dica-inst->ica[id])/0.001;
            inst->ion_ica[indexes[1*pnodecount + id]] += inst->ica[id];
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            vec_rhs[node_id] -= rhs;
            vec_d[node_id] += g;
        }
    }


    /** update state */
    void nrn_state_GRC_CA(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->eca[id] = inst->ion_eca[indexes[0*pnodecount + id]];
            states_GRC_CA(id, pnodecount, data, indexes, thread, nt, ml, v);
        }
    }


    static void check_table_thread_GRC_CA (int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, int tml_id) {
        setup_instance(nt, ml);
        auto* const inst = static_cast<GRC_CA_Instance*>(ml->instance);
        double v = 0;
        check_rate_GRC_CA(id, pnodecount, inst, data, indexes, thread, nt, v);
    }


    /** register channel with the simulator */
    void _GRC_CA_reg() {

        int mech_type = nrn_get_mechtype("GRC_CA");
        GRC_CA_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_GRC_CA, nrn_cur_GRC_CA, nullptr, nrn_state_GRC_CA, nrn_init_GRC_CA, nrn_private_constructor_GRC_CA, nrn_private_destructor_GRC_CA, first_pointer_var_index(), 4);
        GRC_CA_global.ca_type = nrn_get_mechtype("ca_ion");

        thread_mem_init(GRC_CA_global.ext_call_thread);
        _nrn_thread_reg0(mech_type, thread_mem_cleanup);
        _nrn_thread_reg1(mech_type, thread_mem_init);
        _nrn_thread_table_reg(mech_type, check_table_thread_GRC_CA);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 1, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 2, "ca_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
