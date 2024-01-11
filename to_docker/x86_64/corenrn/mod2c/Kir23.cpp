/*********************************************************
Model Name      : Kir2_3
Filename        : Kir23.mod
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
        "Kir2_3",
        "Aalpha_d_Kir2_3",
        "Kalpha_d_Kir2_3",
        "V0alpha_d_Kir2_3",
        "Abeta_d_Kir2_3",
        "Kbeta_d_Kir2_3",
        "V0beta_d_Kir2_3",
        "gkbar_Kir2_3",
        0,
        "ik_Kir2_3",
        "d_inf_Kir2_3",
        "tau_d_Kir2_3",
        "g_Kir2_3",
        "alpha_d_Kir2_3",
        "beta_d_Kir2_3",
        0,
        "d_Kir2_3",
        0,
        0
    };


    /** all global variables */
    struct Kir2_3_Store {
        int k_type{};
        double d0{};
        int reset{};
        int mech_type{};
        int slist1[1]{13};
        int dlist1[1]{15};
        int slist2[1]{13};
        double usetable{1};
        double tmin_rate{};
        double mfac_rate{};
        double t_d_inf[201]{};
        double t_tau_d[201]{};
        ThreadDatum ext_call_thread[3]{};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kir2_3_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kir2_3_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kir2_3_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kir2_3_Store>);
    static_assert(std::is_trivially_destructible_v<Kir2_3_Store>);
    Kir2_3_Store Kir2_3_global;


    /** all mechanism instance variables and global variables */
    struct Kir2_3_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* Aalpha_d{};
        const double* Kalpha_d{};
        const double* V0alpha_d{};
        const double* Abeta_d{};
        const double* Kbeta_d{};
        const double* V0beta_d{};
        const double* gkbar{};
        double* ik{};
        double* d_inf{};
        double* tau_d{};
        double* g{};
        double* alpha_d{};
        double* beta_d{};
        double* d{};
        double* ek{};
        double* Dd{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Kir2_3_Store* global{&Kir2_3_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"usetable_Kir2_3", &Kir2_3_global.usetable},
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
        return 18;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return Kir2_3_global.mech_type;
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
    static void nrn_private_constructor_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kir2_3_Instance{};
        assert(inst->global == &Kir2_3_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kir2_3_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kir2_3_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kir2_3_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kir2_3_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kir2_3_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->Aalpha_d = ml->data+0*pnodecount;
        inst->Kalpha_d = ml->data+1*pnodecount;
        inst->V0alpha_d = ml->data+2*pnodecount;
        inst->Abeta_d = ml->data+3*pnodecount;
        inst->Kbeta_d = ml->data+4*pnodecount;
        inst->V0beta_d = ml->data+5*pnodecount;
        inst->gkbar = ml->data+6*pnodecount;
        inst->ik = ml->data+7*pnodecount;
        inst->d_inf = ml->data+8*pnodecount;
        inst->tau_d = ml->data+9*pnodecount;
        inst->g = ml->data+10*pnodecount;
        inst->alpha_d = ml->data+11*pnodecount;
        inst->beta_d = ml->data+12*pnodecount;
        inst->d = ml->data+13*pnodecount;
        inst->ek = ml->data+14*pnodecount;
        inst->Dd = ml->data+15*pnodecount;
        inst->v_unused = ml->data+16*pnodecount;
        inst->g_unused = ml->data+17*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Kir2_3(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);

        #endif
    }


    inline double alp_d_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double bet_d_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline int rate_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);


    inline int f_rate_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        int ret_f_rate = 0;
        double a_d, b_d, alp_d_in_0, bet_d_in_0;
        {
            double Q10, v_in_0;
            v_in_0 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            alp_d_in_0 = Q10 * inst->Aalpha_d[id] * exp((v_in_0 - inst->V0alpha_d[id]) / inst->Kalpha_d[id]);
        }
        a_d = alp_d_in_0;
        {
            double Q10, v_in_1;
            v_in_1 = arg_v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            bet_d_in_0 = Q10 * inst->Abeta_d[id] * exp((v_in_1 - inst->V0beta_d[id]) / inst->Kbeta_d[id]);
        }
        b_d = bet_d_in_0;
        inst->tau_d[id] = 1.0 / (a_d + b_d);
        inst->d_inf[id] = a_d / (a_d + b_d);
        return ret_f_rate;
    }


    void check_rate_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        if (inst->global->usetable == 0) {
            return;
        }
        static bool make_table = true;
        static double save_Aalpha_d;
        static double save_Kalpha_d;
        static double save_V0alpha_d;
        static double save_Abeta_d;
        static double save_Kbeta_d;
        static double save_V0beta_d;
        static double save_celsius;
        if (save_Aalpha_d != inst->Aalpha_d[id]) {
            make_table = true;
        }
        if (save_Kalpha_d != inst->Kalpha_d[id]) {
            make_table = true;
        }
        if (save_V0alpha_d != inst->V0alpha_d[id]) {
            make_table = true;
        }
        if (save_Abeta_d != inst->Abeta_d[id]) {
            make_table = true;
        }
        if (save_Kbeta_d != inst->Kbeta_d[id]) {
            make_table = true;
        }
        if (save_V0beta_d != inst->V0beta_d[id]) {
            make_table = true;
        }
        if (save_celsius != *(inst->celsius)) {
            make_table = true;
        }
        if (make_table) {
            make_table = false;
            inst->global->tmin_rate =  -100.0;
            double tmax = 100.0;
            double dx = (tmax-inst->global->tmin_rate) / 200.;
            inst->global->mfac_rate = 1./dx;
            double x = inst->global->tmin_rate;
            for (std::size_t i = 0; i < 201; x += dx, i++) {
                f_rate_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v, x);
                inst->global->t_d_inf[i] = inst->d_inf[id];
                inst->global->t_tau_d[i] = inst->tau_d[id];
            }
            save_Aalpha_d = inst->Aalpha_d[id];
            save_Kalpha_d = inst->Kalpha_d[id];
            save_V0alpha_d = inst->V0alpha_d[id];
            save_Abeta_d = inst->Abeta_d[id];
            save_Kbeta_d = inst->Kbeta_d[id];
            save_V0beta_d = inst->V0beta_d[id];
            save_celsius = *(inst->celsius);
        }
    }


    inline int rate_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v){
        if (inst->global->usetable == 0) {
            f_rate_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v, arg_v);
            return 0;
        }
        double xi = inst->global->mfac_rate * (arg_v - inst->global->tmin_rate);
        if (isnan(xi)) {
            inst->d_inf[id] = xi;
            inst->tau_d[id] = xi;
            return 0;
        }
        if (xi <= 0. || xi >= 200.) {
            int index = (xi <= 0.) ? 0 : 200;
            inst->d_inf[id] = inst->global->t_d_inf[index];
            inst->tau_d[id] = inst->global->t_tau_d[index];
            return 0;
        }
        int i = int(xi);
        double theta = xi - double(i);
        inst->d_inf[id] = inst->global->t_d_inf[i] + theta*(inst->global->t_d_inf[i+1]-inst->global->t_d_inf[i]);
        inst->tau_d[id] = inst->global->t_tau_d[i] + theta*(inst->global->t_tau_d[i+1]-inst->global->t_tau_d[i]);
        return 0;
    }


    inline double alp_d_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_alp_d = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
        ret_alp_d = Q10 * inst->Aalpha_d[id] * exp((arg_v - inst->V0alpha_d[id]) / inst->Kalpha_d[id]);
        return ret_alp_d;
    }


    inline double bet_d_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_bet_d = 0.0;
        double Q10;
        Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
        ret_bet_d = Q10 * inst->Abeta_d[id] * exp((arg_v - inst->V0beta_d[id]) / inst->Kbeta_d[id]);
        return ret_bet_d;
    }


    namespace {
        struct _newton_states_Kir2_3 {
            int operator()(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) const {
                auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);
                double* savstate1 = static_cast<double*>(thread[dith1()].pval);
                auto const& slist1 = inst->global->slist1;
                auto const& dlist1 = inst->global->dlist1;
                double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (1*pnodecount);
                rate_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->Dd[id] = (inst->d_inf[id] - inst->d[id]) / inst->tau_d[id];
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

    int states_Kir2_3(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);
        double* savstate1 = (double*) thread[dith1()].pval;
        auto const& slist1 = inst->global->slist1;
        auto& slist2 = inst->global->slist2;
        double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (1*pnodecount);
        for (int i=0; i<1; i++) {
            savstate1[i*pnodecount+id] = data[slist1[i]*pnodecount+id];
        }
        int reset = nrn_newton_thread(static_cast<NewtonSpace*>(*newtonspace1(thread)), 1, slist2, _newton_states_Kir2_3{}, dlist2, id, pnodecount, data, indexes, thread, nt, ml, v);
        return reset;
    }




    /** initialize channel */
    void nrn_init_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);


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
                inst->d[id] = inst->global->d0;
                rate_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v, v);
                inst->d[id] = inst->d_inf[id];
            }
        }
        deriv_advance_flag = 1;
    }


    inline double nrn_current_Kir2_3(int id, int pnodecount, Kir2_3_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double alp_d_in_1, bet_d_in_1;
        inst->g[id] = inst->gkbar[id] * inst->d[id];
        inst->ik[id] = inst->g[id] * (v - inst->ek[id]);
        {
            double Q10, v_in_2;
            v_in_2 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            alp_d_in_1 = Q10 * inst->Aalpha_d[id] * exp((v_in_2 - inst->V0alpha_d[id]) / inst->Kalpha_d[id]);
        }
        inst->alpha_d[id] = alp_d_in_1;
        {
            double Q10, v_in_3;
            v_in_3 = v;
            Q10 = pow(3.0, ((*(inst->celsius) - 20.0) / 10.0));
            bet_d_in_1 = Q10 * inst->Abeta_d[id] * exp((v_in_3 - inst->V0beta_d[id]) / inst->Kbeta_d[id]);
        }
        inst->beta_d[id] = bet_d_in_1;
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            double g = nrn_current_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Kir2_3(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            states_Kir2_3(id, pnodecount, data, indexes, thread, nt, ml, v);
        }
    }


    static void check_table_thread_Kir2_3 (int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, int tml_id) {
        setup_instance(nt, ml);
        auto* const inst = static_cast<Kir2_3_Instance*>(ml->instance);
        double v = 0;
        check_rate_Kir2_3(id, pnodecount, inst, data, indexes, thread, nt, v);
    }


    /** register channel with the simulator */
    void _Kir23_reg() {

        int mech_type = nrn_get_mechtype("Kir2_3");
        Kir2_3_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kir2_3, nrn_cur_Kir2_3, nullptr, nrn_state_Kir2_3, nrn_init_Kir2_3, nrn_private_constructor_Kir2_3, nrn_private_destructor_Kir2_3, first_pointer_var_index(), 4);
        Kir2_3_global.k_type = nrn_get_mechtype("k_ion");

        thread_mem_init(Kir2_3_global.ext_call_thread);
        _nrn_thread_reg0(mech_type, thread_mem_cleanup);
        _nrn_thread_reg1(mech_type, thread_mem_init);
        _nrn_thread_table_reg(mech_type, check_table_thread_Kir2_3);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
