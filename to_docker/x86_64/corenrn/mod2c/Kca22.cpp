/*********************************************************
Model Name      : Kca2_2
Filename        : Kca22.mod
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
#include <newton/newton.hpp>


namespace coreneuron {
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "Kca2_2",
        "gkbar_Kca2_2",
        0,
        "g_Kca2_2",
        "ik_Kca2_2",
        "tcorr_Kca2_2",
        0,
        "c1_Kca2_2",
        "c2_Kca2_2",
        "c3_Kca2_2",
        "c4_Kca2_2",
        "o1_Kca2_2",
        "o2_Kca2_2",
        0,
        0
    };


    /** all global variables */
    struct Kca2_2_Store {
        int ca_type{};
        int k_type{};
        double c10{};
        double c20{};
        double c30{};
        double c40{};
        double o10{};
        double o20{};
        int reset{};
        int mech_type{};
        double Q10{3};
        double diff{3};
        double invc1{0.08};
        double invc2{0.08};
        double invc3{0.2};
        double invo1{1};
        double invo2{0.1};
        double diro1{0.16};
        double diro2{1.2};
        double dirc2{200};
        double dirc3{160};
        double dirc4{80};
        int slist1[6]{4, 5, 6, 7, 8, 9};
        int dlist1[6]{25, 26, 27, 28, 29, 30};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kca2_2_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kca2_2_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kca2_2_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kca2_2_Store>);
    static_assert(std::is_trivially_destructible_v<Kca2_2_Store>);
    Kca2_2_Store Kca2_2_global;


    /** all mechanism instance variables and global variables */
    struct Kca2_2_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gkbar{};
        double* g{};
        double* ik{};
        double* tcorr{};
        double* c1{};
        double* c2{};
        double* c3{};
        double* c4{};
        double* o1{};
        double* o2{};
        double* cai{};
        double* ek{};
        double* invc1_t{};
        double* invc2_t{};
        double* invc3_t{};
        double* invo1_t{};
        double* invo2_t{};
        double* diro1_t{};
        double* diro2_t{};
        double* dirc2_t{};
        double* dirc3_t{};
        double* dirc4_t{};
        double* dirc2_t_ca{};
        double* dirc3_t_ca{};
        double* dirc4_t_ca{};
        double* Dc1{};
        double* Dc2{};
        double* Dc3{};
        double* Dc4{};
        double* Do1{};
        double* Do2{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_cai{};
        const double* ion_cao{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        Kca2_2_Store* global{&Kca2_2_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"Q10_Kca2_2", &Kca2_2_global.Q10},
        {"diff_Kca2_2", &Kca2_2_global.diff},
        {"invc1_Kca2_2", &Kca2_2_global.invc1},
        {"invc2_Kca2_2", &Kca2_2_global.invc2},
        {"invc3_Kca2_2", &Kca2_2_global.invc3},
        {"invo1_Kca2_2", &Kca2_2_global.invo1},
        {"invo2_Kca2_2", &Kca2_2_global.invo2},
        {"diro1_Kca2_2", &Kca2_2_global.diro1},
        {"diro2_Kca2_2", &Kca2_2_global.diro2},
        {"dirc2_Kca2_2", &Kca2_2_global.dirc2},
        {"dirc3_Kca2_2", &Kca2_2_global.dirc3},
        {"dirc4_Kca2_2", &Kca2_2_global.dirc4},
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
        return 33;
    }


    static inline int int_variables_size() {
        return 5;
    }


    static inline int get_mech_type() {
        return Kca2_2_global.mech_type;
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
    static void nrn_private_constructor_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kca2_2_Instance{};
        assert(inst->global == &Kca2_2_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kca2_2_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kca2_2_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kca2_2_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kca2_2_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kca2_2_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gkbar = ml->data+0*pnodecount;
        inst->g = ml->data+1*pnodecount;
        inst->ik = ml->data+2*pnodecount;
        inst->tcorr = ml->data+3*pnodecount;
        inst->c1 = ml->data+4*pnodecount;
        inst->c2 = ml->data+5*pnodecount;
        inst->c3 = ml->data+6*pnodecount;
        inst->c4 = ml->data+7*pnodecount;
        inst->o1 = ml->data+8*pnodecount;
        inst->o2 = ml->data+9*pnodecount;
        inst->cai = ml->data+10*pnodecount;
        inst->ek = ml->data+11*pnodecount;
        inst->invc1_t = ml->data+12*pnodecount;
        inst->invc2_t = ml->data+13*pnodecount;
        inst->invc3_t = ml->data+14*pnodecount;
        inst->invo1_t = ml->data+15*pnodecount;
        inst->invo2_t = ml->data+16*pnodecount;
        inst->diro1_t = ml->data+17*pnodecount;
        inst->diro2_t = ml->data+18*pnodecount;
        inst->dirc2_t = ml->data+19*pnodecount;
        inst->dirc3_t = ml->data+20*pnodecount;
        inst->dirc4_t = ml->data+21*pnodecount;
        inst->dirc2_t_ca = ml->data+22*pnodecount;
        inst->dirc3_t_ca = ml->data+23*pnodecount;
        inst->dirc4_t_ca = ml->data+24*pnodecount;
        inst->Dc1 = ml->data+25*pnodecount;
        inst->Dc2 = ml->data+26*pnodecount;
        inst->Dc3 = ml->data+27*pnodecount;
        inst->Dc4 = ml->data+28*pnodecount;
        inst->Do1 = ml->data+29*pnodecount;
        inst->Do2 = ml->data+30*pnodecount;
        inst->v_unused = ml->data+31*pnodecount;
        inst->g_unused = ml->data+32*pnodecount;
        inst->ion_cai = nt->_data;
        inst->ion_cao = nt->_data;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
    }



    static void nrn_alloc_Kca2_2(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);

        #endif
    }


    inline double temper_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double Q10, double celsius);
    inline int rates_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double cai);
    inline int rate_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double celsius);


    struct functor_Kca2_2_1 {
        NrnThread* nt;
        Kca2_2_Instance* inst;
        int id, pnodecount;
        double v;
        const Datum* indexes;
        double* data;
        ThreadDatum* thread;
        double old_c1, old_c2, old_c3, old_c4, old_o2;

        void initialize() {
            {
                double cai_in_0;
                cai_in_0 = inst->cai[id] / inst->global->diff;
                inst->dirc2_t_ca[id] = inst->dirc2_t[id] * cai_in_0;
                inst->dirc3_t_ca[id] = inst->dirc3_t[id] * cai_in_0;
                inst->dirc4_t_ca[id] = inst->dirc4_t[id] * cai_in_0;
            }
            old_c1 = inst->c1[id];
            old_c2 = inst->c2[id];
            old_c3 = inst->c3[id];
            old_c4 = inst->c4[id];
            old_o2 = inst->o2[id];
        }

        functor_Kca2_2_1(NrnThread* nt, Kca2_2_Instance* inst, int id, int pnodecount, double v, const Datum* indexes, double* data, ThreadDatum* thread) : nt{nt}, inst{inst}, id{id}, pnodecount{pnodecount}, v{v}, indexes{indexes}, data{data}, thread{thread} {}
        void operator()(const Eigen::Matrix<double, 6, 1>& nmodl_eigen_xm, Eigen::Matrix<double, 6, 1>& nmodl_eigen_fm, Eigen::Matrix<double, 6, 6>& nmodl_eigen_jm) const {
            const double* nmodl_eigen_x = nmodl_eigen_xm.data();
            double* nmodl_eigen_j = nmodl_eigen_jm.data();
            double* nmodl_eigen_f = nmodl_eigen_fm.data();
            nmodl_eigen_f[static_cast<int>(0)] =  -nmodl_eigen_x[static_cast<int>(0)] * inst->dirc2_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(0)] + nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->invc1_t[id] + old_c1;
            nmodl_eigen_j[static_cast<int>(0)] =  -inst->dirc2_t_ca[id] * nt->_dt - 1.0;
            nmodl_eigen_j[static_cast<int>(6)] = nt->_dt * inst->invc1_t[id];
            nmodl_eigen_j[static_cast<int>(12)] = 0.0;
            nmodl_eigen_j[static_cast<int>(18)] = 0.0;
            nmodl_eigen_j[static_cast<int>(24)] = 0.0;
            nmodl_eigen_j[static_cast<int>(30)] = 0.0;
            nmodl_eigen_f[static_cast<int>(1)] = nmodl_eigen_x[static_cast<int>(0)] * inst->dirc2_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * inst->dirc3_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->invc1_t[id] - nmodl_eigen_x[static_cast<int>(1)] + nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->invc2_t[id] + old_c2;
            nmodl_eigen_j[static_cast<int>(1)] = inst->dirc2_t_ca[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(7)] =  -inst->dirc3_t_ca[id] * nt->_dt - nt->_dt * inst->invc1_t[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(13)] = nt->_dt * inst->invc2_t[id];
            nmodl_eigen_j[static_cast<int>(19)] = 0.0;
            nmodl_eigen_j[static_cast<int>(25)] = 0.0;
            nmodl_eigen_j[static_cast<int>(31)] = 0.0;
            nmodl_eigen_f[static_cast<int>(2)] = nmodl_eigen_x[static_cast<int>(1)] * inst->dirc3_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->dirc4_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->diro1_t[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->invc2_t[id] - nmodl_eigen_x[static_cast<int>(2)] + nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->invc3_t[id] + nmodl_eigen_x[static_cast<int>(4)] * nt->_dt * inst->invo1_t[id] + old_c3;
            nmodl_eigen_j[static_cast<int>(2)] = 0.0;
            nmodl_eigen_j[static_cast<int>(8)] = inst->dirc3_t_ca[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(14)] =  -inst->dirc4_t_ca[id] * nt->_dt - inst->diro1_t[id] * nt->_dt - nt->_dt * inst->invc2_t[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(20)] = nt->_dt * inst->invc3_t[id];
            nmodl_eigen_j[static_cast<int>(26)] = nt->_dt * inst->invo1_t[id];
            nmodl_eigen_j[static_cast<int>(32)] = 0.0;
            nmodl_eigen_f[static_cast<int>(3)] = nmodl_eigen_x[static_cast<int>(2)] * inst->dirc4_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * inst->diro2_t[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->invc3_t[id] - nmodl_eigen_x[static_cast<int>(3)] + nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->invo2_t[id] + old_c4;
            nmodl_eigen_j[static_cast<int>(3)] = 0.0;
            nmodl_eigen_j[static_cast<int>(9)] = 0.0;
            nmodl_eigen_j[static_cast<int>(15)] = inst->dirc4_t_ca[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(21)] =  -inst->diro2_t[id] * nt->_dt - nt->_dt * inst->invc3_t[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(27)] = 0.0;
            nmodl_eigen_j[static_cast<int>(33)] = nt->_dt * inst->invo2_t[id];
            nmodl_eigen_f[static_cast<int>(4)] =  -nmodl_eigen_x[static_cast<int>(0)] - nmodl_eigen_x[static_cast<int>(1)] - nmodl_eigen_x[static_cast<int>(2)] - nmodl_eigen_x[static_cast<int>(3)] - nmodl_eigen_x[static_cast<int>(4)] - nmodl_eigen_x[static_cast<int>(5)] + 1.0;
            nmodl_eigen_j[static_cast<int>(4)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(10)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(16)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(22)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(28)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(34)] =  -1.0;
            nmodl_eigen_f[static_cast<int>(5)] = nmodl_eigen_x[static_cast<int>(3)] * inst->diro2_t[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->invo2_t[id] - nmodl_eigen_x[static_cast<int>(5)] + old_o2;
            nmodl_eigen_j[static_cast<int>(5)] = 0.0;
            nmodl_eigen_j[static_cast<int>(11)] = 0.0;
            nmodl_eigen_j[static_cast<int>(17)] = 0.0;
            nmodl_eigen_j[static_cast<int>(23)] = inst->diro2_t[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(29)] = 0.0;
            nmodl_eigen_j[static_cast<int>(35)] =  -nt->_dt * inst->invo2_t[id] - 1.0;
        }

        void finalize() {
        }
    };


    struct functor_Kca2_2_0 {
        NrnThread* nt;
        Kca2_2_Instance* inst;
        int id, pnodecount;
        double v;
        const Datum* indexes;
        double* data;
        ThreadDatum* thread;
        double old_c1, old_c2, old_c3, old_c4, old_o2;

        void initialize() {
            ;
            {
                double cai_in_1;
                cai_in_1 = inst->cai[id] / inst->global->diff;
                inst->dirc2_t_ca[id] = inst->dirc2_t[id] * cai_in_1;
                inst->dirc3_t_ca[id] = inst->dirc3_t[id] * cai_in_1;
                inst->dirc4_t_ca[id] = inst->dirc4_t[id] * cai_in_1;
            }
            old_c1 = inst->c1[id];
            old_c2 = inst->c2[id];
            old_c3 = inst->c3[id];
            old_c4 = inst->c4[id];
            old_o2 = inst->o2[id];
        }

        functor_Kca2_2_0(NrnThread* nt, Kca2_2_Instance* inst, int id, int pnodecount, double v, const Datum* indexes, double* data, ThreadDatum* thread) : nt{nt}, inst{inst}, id{id}, pnodecount{pnodecount}, v{v}, indexes{indexes}, data{data}, thread{thread} {}
        void operator()(const Eigen::Matrix<double, 6, 1>& nmodl_eigen_xm, Eigen::Matrix<double, 6, 1>& nmodl_eigen_fm, Eigen::Matrix<double, 6, 6>& nmodl_eigen_jm) const {
            const double* nmodl_eigen_x = nmodl_eigen_xm.data();
            double* nmodl_eigen_j = nmodl_eigen_jm.data();
            double* nmodl_eigen_f = nmodl_eigen_fm.data();
            nmodl_eigen_f[static_cast<int>(0)] =  -nmodl_eigen_x[static_cast<int>(0)] * inst->dirc2_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(0)] + nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->invc1_t[id] + old_c1;
            nmodl_eigen_j[static_cast<int>(0)] =  -inst->dirc2_t_ca[id] * nt->_dt - 1.0;
            nmodl_eigen_j[static_cast<int>(6)] = nt->_dt * inst->invc1_t[id];
            nmodl_eigen_j[static_cast<int>(12)] = 0.0;
            nmodl_eigen_j[static_cast<int>(18)] = 0.0;
            nmodl_eigen_j[static_cast<int>(24)] = 0.0;
            nmodl_eigen_j[static_cast<int>(30)] = 0.0;
            nmodl_eigen_f[static_cast<int>(1)] = nmodl_eigen_x[static_cast<int>(0)] * inst->dirc2_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * inst->dirc3_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->invc1_t[id] - nmodl_eigen_x[static_cast<int>(1)] + nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->invc2_t[id] + old_c2;
            nmodl_eigen_j[static_cast<int>(1)] = inst->dirc2_t_ca[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(7)] =  -inst->dirc3_t_ca[id] * nt->_dt - nt->_dt * inst->invc1_t[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(13)] = nt->_dt * inst->invc2_t[id];
            nmodl_eigen_j[static_cast<int>(19)] = 0.0;
            nmodl_eigen_j[static_cast<int>(25)] = 0.0;
            nmodl_eigen_j[static_cast<int>(31)] = 0.0;
            nmodl_eigen_f[static_cast<int>(2)] = nmodl_eigen_x[static_cast<int>(1)] * inst->dirc3_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->dirc4_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->diro1_t[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->invc2_t[id] - nmodl_eigen_x[static_cast<int>(2)] + nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->invc3_t[id] + nmodl_eigen_x[static_cast<int>(4)] * nt->_dt * inst->invo1_t[id] + old_c3;
            nmodl_eigen_j[static_cast<int>(2)] = 0.0;
            nmodl_eigen_j[static_cast<int>(8)] = inst->dirc3_t_ca[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(14)] =  -inst->dirc4_t_ca[id] * nt->_dt - inst->diro1_t[id] * nt->_dt - nt->_dt * inst->invc2_t[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(20)] = nt->_dt * inst->invc3_t[id];
            nmodl_eigen_j[static_cast<int>(26)] = nt->_dt * inst->invo1_t[id];
            nmodl_eigen_j[static_cast<int>(32)] = 0.0;
            nmodl_eigen_f[static_cast<int>(3)] = nmodl_eigen_x[static_cast<int>(2)] * inst->dirc4_t_ca[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * inst->diro2_t[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->invc3_t[id] - nmodl_eigen_x[static_cast<int>(3)] + nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->invo2_t[id] + old_c4;
            nmodl_eigen_j[static_cast<int>(3)] = 0.0;
            nmodl_eigen_j[static_cast<int>(9)] = 0.0;
            nmodl_eigen_j[static_cast<int>(15)] = inst->dirc4_t_ca[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(21)] =  -inst->diro2_t[id] * nt->_dt - nt->_dt * inst->invc3_t[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(27)] = 0.0;
            nmodl_eigen_j[static_cast<int>(33)] = nt->_dt * inst->invo2_t[id];
            nmodl_eigen_f[static_cast<int>(4)] =  -nmodl_eigen_x[static_cast<int>(0)] - nmodl_eigen_x[static_cast<int>(1)] - nmodl_eigen_x[static_cast<int>(2)] - nmodl_eigen_x[static_cast<int>(3)] - nmodl_eigen_x[static_cast<int>(4)] - nmodl_eigen_x[static_cast<int>(5)] + 1.0;
            nmodl_eigen_j[static_cast<int>(4)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(10)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(16)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(22)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(28)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(34)] =  -1.0;
            nmodl_eigen_f[static_cast<int>(5)] = nmodl_eigen_x[static_cast<int>(3)] * inst->diro2_t[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->invo2_t[id] - nmodl_eigen_x[static_cast<int>(5)] + old_o2;
            nmodl_eigen_j[static_cast<int>(5)] = 0.0;
            nmodl_eigen_j[static_cast<int>(11)] = 0.0;
            nmodl_eigen_j[static_cast<int>(17)] = 0.0;
            nmodl_eigen_j[static_cast<int>(23)] = inst->diro2_t[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(29)] = 0.0;
            nmodl_eigen_j[static_cast<int>(35)] =  -nt->_dt * inst->invo2_t[id] - 1.0;
        }

        void finalize() {
        }
    };


    inline int rates_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double cai) {
        int ret_rates = 0;
        inst->dirc2_t_ca[id] = inst->dirc2_t[id] * inst->cai[id];
        inst->dirc3_t_ca[id] = inst->dirc3_t[id] * inst->cai[id];
        inst->dirc4_t_ca[id] = inst->dirc4_t[id] * inst->cai[id];
        return ret_rates;
    }


    inline int rate_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double celsius) {
        int ret_rate = 0;
        double temper_in_0;
        {
            double Q10_in_0, celsius_in_0;
            Q10_in_0 = inst->global->Q10;
            celsius_in_0 = *(inst->celsius);
            temper_in_0 = pow(Q10_in_0, ((celsius_in_0 - 23.0) / 10.0));
        }
        inst->tcorr[id] = temper_in_0;
        inst->invc1_t[id] = inst->global->invc1 * inst->tcorr[id];
        inst->invc2_t[id] = inst->global->invc2 * inst->tcorr[id];
        inst->invc3_t[id] = inst->global->invc3 * inst->tcorr[id];
        inst->invo1_t[id] = inst->global->invo1 * inst->tcorr[id];
        inst->invo2_t[id] = inst->global->invo2 * inst->tcorr[id];
        inst->diro1_t[id] = inst->global->diro1 * inst->tcorr[id];
        inst->diro2_t[id] = inst->global->diro2 * inst->tcorr[id];
        inst->dirc2_t[id] = inst->global->dirc2 * inst->tcorr[id];
        inst->dirc3_t[id] = inst->global->dirc3 * inst->tcorr[id];
        inst->dirc4_t[id] = inst->global->dirc4 * inst->tcorr[id];
        return ret_rate;
    }


    inline double temper_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double Q10, double celsius) {
        double ret_temper = 0.0;
        ret_temper = pow(inst->global->Q10, ((*(inst->celsius) - 23.0) / 10.0));
        return ret_temper;
    }


    /** initialize channel */
    void nrn_init_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            double _save_prev_dt = nt->_dt;
            nt->_dt = 1000000000;
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->cai[id] = inst->ion_cai[indexes[0*pnodecount + id]];
                inst->ek[id] = inst->ion_ek[indexes[2*pnodecount + id]];
                inst->c1[id] = inst->global->c10;
                inst->c2[id] = inst->global->c20;
                inst->c3[id] = inst->global->c30;
                inst->c4[id] = inst->global->c40;
                inst->o1[id] = inst->global->o10;
                inst->o2[id] = inst->global->o20;
                {
                    double temper_in_0, celsius_in_1;
                    celsius_in_1 = *(inst->celsius);
                    {
                        double Q10_in_0, celsius_in_0;
                        Q10_in_0 = inst->global->Q10;
                        celsius_in_0 = celsius_in_1;
                        temper_in_0 = pow(Q10_in_0, ((celsius_in_0 - 23.0) / 10.0));
                    }
                    inst->tcorr[id] = temper_in_0;
                    inst->invc1_t[id] = inst->global->invc1 * inst->tcorr[id];
                    inst->invc2_t[id] = inst->global->invc2 * inst->tcorr[id];
                    inst->invc3_t[id] = inst->global->invc3 * inst->tcorr[id];
                    inst->invo1_t[id] = inst->global->invo1 * inst->tcorr[id];
                    inst->invo2_t[id] = inst->global->invo2 * inst->tcorr[id];
                    inst->diro1_t[id] = inst->global->diro1 * inst->tcorr[id];
                    inst->diro2_t[id] = inst->global->diro2 * inst->tcorr[id];
                    inst->dirc2_t[id] = inst->global->dirc2 * inst->tcorr[id];
                    inst->dirc3_t[id] = inst->global->dirc3 * inst->tcorr[id];
                    inst->dirc4_t[id] = inst->global->dirc4 * inst->tcorr[id];
                }
                                
                Eigen::Matrix<double, 6, 1> nmodl_eigen_xm;
                double* nmodl_eigen_x = nmodl_eigen_xm.data();
                nmodl_eigen_x[static_cast<int>(0)] = inst->c1[id];
                nmodl_eigen_x[static_cast<int>(1)] = inst->c2[id];
                nmodl_eigen_x[static_cast<int>(2)] = inst->c3[id];
                nmodl_eigen_x[static_cast<int>(3)] = inst->c4[id];
                nmodl_eigen_x[static_cast<int>(4)] = inst->o1[id];
                nmodl_eigen_x[static_cast<int>(5)] = inst->o2[id];
                // call newton solver
                functor_Kca2_2_0 newton_functor(nt, inst, id, pnodecount, v, indexes, data, thread);
                newton_functor.initialize();
                int newton_iterations = nmodl::newton::newton_solver(nmodl_eigen_xm, newton_functor);
                if (newton_iterations < 0) assert(false && "Newton solver did not converge!");
                inst->c1[id] = nmodl_eigen_x[static_cast<int>(0)];
                inst->c2[id] = nmodl_eigen_x[static_cast<int>(1)];
                inst->c3[id] = nmodl_eigen_x[static_cast<int>(2)];
                inst->c4[id] = nmodl_eigen_x[static_cast<int>(3)];
                inst->o1[id] = nmodl_eigen_x[static_cast<int>(4)];
                inst->o2[id] = nmodl_eigen_x[static_cast<int>(5)];
                newton_functor.finalize();


            }
            nt->_dt = _save_prev_dt;
        }
    }


    inline double nrn_current_Kca2_2(int id, int pnodecount, Kca2_2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->g[id] = inst->gkbar[id] * (inst->o1[id] + inst->o2[id]);
        inst->ik[id] = inst->g[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->cai[id] = inst->ion_cai[indexes[0*pnodecount + id]];
            inst->ek[id] = inst->ion_ek[indexes[2*pnodecount + id]];
            double g = nrn_current_Kca2_2(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kca2_2(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            inst->ion_dikdv[indexes[4*pnodecount + id]] += (dik-inst->ik[id])/0.001;
            inst->ion_ik[indexes[3*pnodecount + id]] += inst->ik[id];
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            vec_rhs[node_id] -= rhs;
            vec_d[node_id] += g;
        }
    }


    /** update state */
    void nrn_state_Kca2_2(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca2_2_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->cai[id] = inst->ion_cai[indexes[0*pnodecount + id]];
            inst->ek[id] = inst->ion_ek[indexes[2*pnodecount + id]];
            
            Eigen::Matrix<double, 6, 1> nmodl_eigen_xm;
            double* nmodl_eigen_x = nmodl_eigen_xm.data();
            nmodl_eigen_x[static_cast<int>(0)] = inst->c1[id];
            nmodl_eigen_x[static_cast<int>(1)] = inst->c2[id];
            nmodl_eigen_x[static_cast<int>(2)] = inst->c3[id];
            nmodl_eigen_x[static_cast<int>(3)] = inst->c4[id];
            nmodl_eigen_x[static_cast<int>(4)] = inst->o1[id];
            nmodl_eigen_x[static_cast<int>(5)] = inst->o2[id];
            // call newton solver
            functor_Kca2_2_1 newton_functor(nt, inst, id, pnodecount, v, indexes, data, thread);
            newton_functor.initialize();
            int newton_iterations = nmodl::newton::newton_solver(nmodl_eigen_xm, newton_functor);
            if (newton_iterations < 0) assert(false && "Newton solver did not converge!");
            inst->c1[id] = nmodl_eigen_x[static_cast<int>(0)];
            inst->c2[id] = nmodl_eigen_x[static_cast<int>(1)];
            inst->c3[id] = nmodl_eigen_x[static_cast<int>(2)];
            inst->c4[id] = nmodl_eigen_x[static_cast<int>(3)];
            inst->o1[id] = nmodl_eigen_x[static_cast<int>(4)];
            inst->o2[id] = nmodl_eigen_x[static_cast<int>(5)];
            newton_functor.finalize();

        }
    }


    /** register channel with the simulator */
    void _Kca22_reg() {

        int mech_type = nrn_get_mechtype("Kca2_2");
        Kca2_2_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kca2_2, nrn_cur_Kca2_2, nullptr, nrn_state_Kca2_2, nrn_init_Kca2_2, nrn_private_constructor_Kca2_2, nrn_private_destructor_Kca2_2, first_pointer_var_index(), 1);
        Kca2_2_global.ca_type = nrn_get_mechtype("ca_ion");
        Kca2_2_global.k_type = nrn_get_mechtype("k_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 1, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_dparam_semantics(mech_type, 3, "k_ion");
        hoc_register_dparam_semantics(mech_type, 4, "k_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
