/*********************************************************
Model Name      : GRC_KM
Filename        : GRC_KM.mod
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
        "GRC_KM",
        "Aalpha_n_GRC_KM",
        "Kalpha_n_GRC_KM",
        "V0alpha_n_GRC_KM",
        "Abeta_n_GRC_KM",
        "Kbeta_n_GRC_KM",
        "V0beta_n_GRC_KM",
        "V0_ninf_GRC_KM",
        "B_ninf_GRC_KM",
        "gkbar_GRC_KM",
        0,
        "ik_GRC_KM",
        "n_inf_GRC_KM",
        "tau_n_GRC_KM",
        "g_GRC_KM",
        "alpha_n_GRC_KM",
        "beta_n_GRC_KM",
        0,
        "n_GRC_KM",
        0,
        0
    };


    /** all global variables */
    struct GRC_KM_Store {
        int k_type{};
        double n0{};
        int reset{};
        int mech_type{};
        int slist1[1]{15};
        int dlist1[1]{17};
        int slist2[1]{15};
        double usetable{1};
        double tmin_rate{};
        double mfac_rate{};
        double t_n_inf[13001]{};
        double t_tau_n[13001]{};
        ThreadDatum ext_call_thread[3]{};
    };
    static_assert(std::is_trivially_copy_constructible_v<GRC_KM_Store>);
    static_assert(std::is_trivially_move_constructible_v<GRC_KM_Store>);
    static_assert(std::is_trivially_copy_assignable_v<GRC_KM_Store>);
    static_assert(std::is_trivially_move_assignable_v<GRC_KM_Store>);
    static_assert(std::is_trivially_destructible_v<GRC_KM_Store>);
    GRC_KM_Store GRC_KM_global;


    /** all mechanism instance variables and global variables */
    struct GRC_KM_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* Aalpha_n{};
        const double* Kalpha_n{};
        const double* V0alpha_n{};
        const double* Abeta_n{};
        const double* Kbeta_n{};
        const double* V0beta_n{};
        const double* V0_ninf{};
        const double* B_ninf{};
        const double* gkbar{};
        double* ik{};
        double* n_inf{};
        double* tau_n{};
        double* g{};
        double* alpha_n{};
        double* beta_n{};
        double* n{};
        double* ek{};
        double* Dn{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        GRC_KM_Store* global{&GRC_KM_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"usetable_GRC_KM", &GRC_KM_global.usetable},
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
        return 20;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return GRC_KM_global.mech_type;
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
    static void nrn_private_constructor_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new GRC_KM_Instance{};
        assert(inst->global == &GRC_KM_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(GRC_KM_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &GRC_KM_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(GRC_KM_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &GRC_KM_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(GRC_KM_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->Aalpha_n = ml->data+0*pnodecount;
        inst->Kalpha_n = ml->data+1*pnodecount;
        inst->V0alpha_n = ml->data+2*pnodecount;
        inst->Abeta_n = ml->data+3*pnodecount;
        inst->Kbeta_n = ml->data+4*pnodecount;
        inst->V0beta_n = ml->data+5*pnodecount;
        inst->V0_ninf = ml->data+6*pnodecount;
        inst->B_ninf = ml->data+7*pnodecount;
        inst->gkbar = ml->data+8*pnodecount;
        inst->ik = ml->data+9*pnodecount;
        inst->n_inf = ml->data+10*pnodecount;
        inst->tau_n = ml->data+11*pnodecount;
        inst->g = ml->data+12*pnodecount;
        inst->alpha_n = ml->data+13*pnodecount;
        inst->beta_n = ml->data+14*pnodecount;
        inst->n = ml->data+15*pnodecount;
        inst->ek = ml->data+16*pnodecount;
        inst->Dn = ml->data+17*pnodecount;
        inst->v_unused = ml->data+18*pnodecount;
        inst->g_unused = ml->data+19*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_GRC_KM(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);

        #endif
    }


    inline double alp_n_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double bet_n_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline int rate_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int f_rate_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_f_rate = 0;
        double a_n, b_n, alp_n_in_0, bet_n_in_0;
        {
            double Q10, v_in_0;
            v_in_0 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 22.0) / 10.0));
            alp_n_in_0 = Q10 * inst->Aalpha_n[id] * exp((v_in_0 - inst->V0alpha_n[id]) / inst->Kalpha_n[id]);
        }
        a_n = alp_n_in_0;
        {
            double Q10, v_in_1;
            v_in_1 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 22.0) / 10.0));
            bet_n_in_0 = Q10 * inst->Abeta_n[id] * exp((v_in_1 - inst->V0beta_n[id]) / inst->Kbeta_n[id]);
        }
        b_n = bet_n_in_0;
        inst->tau_n[id] = 1.0 / (a_n + b_n);
        inst->n_inf[id] = 1.0 / (1.0 + exp( -(arg_v - inst->V0_ninf[id]) / inst->B_ninf[id]));
        return ret_f_rate;
    }


    void check_rate_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        if (inst->global->usetable == 0) {
            return;
        }
        static bool make_table = true;
        static double save_Aalpha_n;
        static double save_Kalpha_n;
        static double save_V0alpha_n;
        static double save_Abeta_n;
        static double save_Kbeta_n;
        static double save_V0beta_n;
        static double save_V0_ninf;
        static double save_B_ninf;
        static double save_celsius;
        if (save_Aalpha_n != inst->Aalpha_n[id]) {
            make_table = true;
        }
        if (save_Kalpha_n != inst->Kalpha_n[id]) {
            make_table = true;
        }
        if (save_V0alpha_n != inst->V0alpha_n[id]) {
            make_table = true;
        }
        if (save_Abeta_n != inst->Abeta_n[id]) {
            make_table = true;
        }
        if (save_Kbeta_n != inst->Kbeta_n[id]) {
            make_table = true;
        }
        if (save_V0beta_n != inst->V0beta_n[id]) {
            make_table = true;
        }
        if (save_V0_ninf != inst->V0_ninf[id]) {
            make_table = true;
        }
        if (save_B_ninf != inst->B_ninf[id]) {
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
                f_rate_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v, x);
                inst->global->t_n_inf[i] = inst->n_inf[id];
                inst->global->t_tau_n[i] = inst->tau_n[id];
            }
            save_Aalpha_n = inst->Aalpha_n[id];
            save_Kalpha_n = inst->Kalpha_n[id];
            save_V0alpha_n = inst->V0alpha_n[id];
            save_Abeta_n = inst->Abeta_n[id];
            save_Kbeta_n = inst->Kbeta_n[id];
            save_V0beta_n = inst->V0beta_n[id];
            save_V0_ninf = inst->V0_ninf[id];
            save_B_ninf = inst->B_ninf[id];
            save_celsius = *(inst->celsius);
        }
    }


    inline int rate_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v){
        if (inst->global->usetable == 0) {
            f_rate_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v, arg_v);
            return 0;
        }
        double xi = inst->global->mfac_rate * (arg_v - inst->global->tmin_rate);
        if (isnan(xi)) {
            inst->n_inf[id] = xi;
            inst->tau_n[id] = xi;
            return 0;
        }
        if (xi <= 0. || xi >= 13000.) {
            int index = (xi <= 0.) ? 0 : 13000;
            inst->n_inf[id] = inst->global->t_n_inf[index];
            inst->tau_n[id] = inst->global->t_tau_n[index];
            return 0;
        }
        int i = int(xi);
        double theta = xi - double(i);
        inst->n_inf[id] = inst->global->t_n_inf[i] + theta*(inst->global->t_n_inf[i+1]-inst->global->t_n_inf[i]);
        inst->tau_n[id] = inst->global->t_tau_n[i] + theta*(inst->global->t_tau_n[i+1]-inst->global->t_tau_n[i]);
        return 0;
    }


    inline double alp_n_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alp_n = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 22.0) / 10.0));
        ret_alp_n = Q10 * inst->Aalpha_n[id] * exp((arg_v - inst->V0alpha_n[id]) / inst->Kalpha_n[id]);
        return ret_alp_n;
    }


    inline double bet_n_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_bet_n = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 22.0) / 10.0));
        ret_bet_n = Q10 * inst->Abeta_n[id] * exp((arg_v - inst->V0beta_n[id]) / inst->Kbeta_n[id]);
        return ret_bet_n;
    }


    namespace {
        struct _newton_states_GRC_KM {
            int operator()(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) const {
                auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);
                double* savstate1 = static_cast<double*>(thread[dith1()].pval);
                auto const& slist1 = inst->global->slist1;
                auto const& dlist1 = inst->global->dlist1;
                double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (1*pnodecount);
                rate_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->Dn[id] = (inst->n_inf[id] - inst->n[id]) / inst->tau_n[id];
                int counter = -1;
                for (int i=0; i<1; i++) {
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

    int states_GRC_KM(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);
        double* savstate1 = (double*) thread[dith1()].pval;
        auto const& slist1 = inst->global->slist1;
        auto& slist2 = inst->global->slist2;
        double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (1*pnodecount);
        for (int i=0; i<1; i++) {
            savstate1[i*pnodecount+id] = data[slist1[i]*pnodecount+id];
        }
        int reset = nrn_newton_thread(static_cast<NewtonSpace*>(*newtonspace1(thread)), 1, slist2, _newton_states_GRC_KM{}, dlist2, id, pnodecount, data, indexes, thread, nt, ml, v);
        return reset;
    }




    /** initialize channel */
    void nrn_init_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);


        int& deriv_advance_flag = *deriv1_advance(thread);
        deriv_advance_flag = 0;
        auto ns = newtonspace1(thread);
        auto& th = thread[dith1()];
        if (*ns == nullptr) {
            int vec_size = 2*1*pnodecount*sizeof(double);
            double* vec = makevector(vec_size);
            th.pval = vec;
            *ns = nrn_cons_newtonspace(1, pnodecount);
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
                inst->n[id] = inst->global->n0;
                rate_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->n[id] = inst->n_inf[id];
            }
        }
        deriv_advance_flag = 1;
    }


    inline double nrn_current_GRC_KM(int id, int pnodecount, GRC_KM_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double alp_n_in_1, bet_n_in_1;
        inst->g[id] = inst->gkbar[id] * inst->n[id];
        inst->ik[id] = inst->g[id] * (v - inst->ek[id]);
        {
            double Q10, v_in_2;
            v_in_2 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 22.0) / 10.0));
            alp_n_in_1 = Q10 * inst->Aalpha_n[id] * exp((v_in_2 - inst->V0alpha_n[id]) / inst->Kalpha_n[id]);
        }
        inst->alpha_n[id] = alp_n_in_1;
        {
            double Q10, v_in_3;
            v_in_3 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 22.0) / 10.0));
            bet_n_in_1 = Q10 * inst->Abeta_n[id] * exp((v_in_3 - inst->V0beta_n[id]) / inst->Kbeta_n[id]);
        }
        inst->beta_n[id] = bet_n_in_1;
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_GRC_KM(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            states_GRC_KM(id, pnodecount, data, indexes, thread, nt, ml, v);
        }
    }


    static void check_table_thread_GRC_KM (int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, int tml_id) {
        setup_instance(nt, ml);
        auto* const inst = static_cast<GRC_KM_Instance*>(ml->instance);
        double v = 0;
        check_rate_GRC_KM(id, pnodecount, inst, data, indexes, thread, nt, v);
    }


    /** register channel with the simulator */
    void _GRC_KM_reg() {

        int mech_type = nrn_get_mechtype("GRC_KM");
        GRC_KM_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_GRC_KM, nrn_cur_GRC_KM, nullptr, nrn_state_GRC_KM, nrn_init_GRC_KM, nrn_private_constructor_GRC_KM, nrn_private_destructor_GRC_KM, first_pointer_var_index(), 4);
        GRC_KM_global.k_type = nrn_get_mechtype("k_ion");

        thread_mem_init(GRC_KM_global.ext_call_thread);
        _nrn_thread_reg0(mech_type, thread_mem_cleanup);
        _nrn_thread_reg1(mech_type, thread_mem_init);
        _nrn_thread_table_reg(mech_type, check_table_thread_GRC_KM);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
