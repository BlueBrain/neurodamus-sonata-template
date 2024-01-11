/*********************************************************
Model Name      : Kca1_1
Filename        : Kca11.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Thu Jan 11 15:15:21 2024
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


    /** constants used in nmodl from UNITS */
    static const double FARADAY = 0x1.81f0fae775425p+6;
    static const double R = 0x1.0a1013e8990bep+3;
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "Kca1_1",
        "gbar_Kca1_1",
        0,
        "ik_Kca1_1",
        "g_Kca1_1",
        0,
        "C0_Kca1_1",
        "C1_Kca1_1",
        "C2_Kca1_1",
        "C3_Kca1_1",
        "C4_Kca1_1",
        "O0_Kca1_1",
        "O1_Kca1_1",
        "O2_Kca1_1",
        "O3_Kca1_1",
        "O4_Kca1_1",
        0,
        0
    };


    /** all global variables */
    struct Kca1_1_Store {
        int k_type{};
        int ca_type{};
        double C00{};
        double C10{};
        double C20{};
        double C30{};
        double C40{};
        double O00{};
        double O10{};
        double O20{};
        double O30{};
        double O40{};
        int reset{};
        int mech_type{};
        double Qo{0.73};
        double Qc{-0.67};
        double k1{1000};
        double onoffrate{1};
        double L0{1806};
        double Kc{0.011};
        double Ko{0.0011};
        double pf0{0.00239};
        double pf1{0.007};
        double pf2{0.04};
        double pf3{0.295};
        double pf4{0.557};
        double pb0{3.936};
        double pb1{1.152};
        double pb2{0.659};
        double pb3{0.486};
        double pb4{0.092};
        double q10{3};
        int slist1[10]{3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        int dlist1[10]{41, 42, 43, 44, 45, 46, 47, 48, 49, 50};
    };
    static_assert(std::is_trivially_copy_constructible_v<Kca1_1_Store>);
    static_assert(std::is_trivially_move_constructible_v<Kca1_1_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Kca1_1_Store>);
    static_assert(std::is_trivially_move_assignable_v<Kca1_1_Store>);
    static_assert(std::is_trivially_destructible_v<Kca1_1_Store>);
    Kca1_1_Store Kca1_1_global;


    /** all mechanism instance variables and global variables */
    struct Kca1_1_Instance  {
        double* celsius{&coreneuron::celsius};
        const double* gbar{};
        double* ik{};
        double* g{};
        double* C0{};
        double* C1{};
        double* C2{};
        double* C3{};
        double* C4{};
        double* O0{};
        double* O1{};
        double* O2{};
        double* O3{};
        double* O4{};
        double* c01{};
        double* c12{};
        double* c23{};
        double* c34{};
        double* o01{};
        double* o12{};
        double* o23{};
        double* o34{};
        double* f0{};
        double* f1{};
        double* f2{};
        double* f3{};
        double* f4{};
        double* c10{};
        double* c21{};
        double* c32{};
        double* c43{};
        double* o10{};
        double* o21{};
        double* o32{};
        double* o43{};
        double* b0{};
        double* b1{};
        double* b2{};
        double* b3{};
        double* b4{};
        double* cai{};
        double* ek{};
        double* DC0{};
        double* DC1{};
        double* DC2{};
        double* DC3{};
        double* DC4{};
        double* DO0{};
        double* DO1{};
        double* DO2{};
        double* DO3{};
        double* DO4{};
        double* v_unused{};
        double* g_unused{};
        const double* ion_ek{};
        double* ion_ik{};
        double* ion_dikdv{};
        const double* ion_cai{};
        const double* ion_cao{};
        Kca1_1_Store* global{&Kca1_1_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"Qo_Kca1_1", &Kca1_1_global.Qo},
        {"Qc_Kca1_1", &Kca1_1_global.Qc},
        {"k1_Kca1_1", &Kca1_1_global.k1},
        {"onoffrate_Kca1_1", &Kca1_1_global.onoffrate},
        {"L0_Kca1_1", &Kca1_1_global.L0},
        {"Kc_Kca1_1", &Kca1_1_global.Kc},
        {"Ko_Kca1_1", &Kca1_1_global.Ko},
        {"pf0_Kca1_1", &Kca1_1_global.pf0},
        {"pf1_Kca1_1", &Kca1_1_global.pf1},
        {"pf2_Kca1_1", &Kca1_1_global.pf2},
        {"pf3_Kca1_1", &Kca1_1_global.pf3},
        {"pf4_Kca1_1", &Kca1_1_global.pf4},
        {"pb0_Kca1_1", &Kca1_1_global.pb0},
        {"pb1_Kca1_1", &Kca1_1_global.pb1},
        {"pb2_Kca1_1", &Kca1_1_global.pb2},
        {"pb3_Kca1_1", &Kca1_1_global.pb3},
        {"pb4_Kca1_1", &Kca1_1_global.pb4},
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
        return 53;
    }


    static inline int int_variables_size() {
        return 5;
    }


    static inline int get_mech_type() {
        return Kca1_1_global.mech_type;
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
    static void nrn_private_constructor_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Kca1_1_Instance{};
        assert(inst->global == &Kca1_1_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Kca1_1_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kca1_1_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kca1_1_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Kca1_1_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Kca1_1_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->gbar = ml->data+0*pnodecount;
        inst->ik = ml->data+1*pnodecount;
        inst->g = ml->data+2*pnodecount;
        inst->C0 = ml->data+3*pnodecount;
        inst->C1 = ml->data+4*pnodecount;
        inst->C2 = ml->data+5*pnodecount;
        inst->C3 = ml->data+6*pnodecount;
        inst->C4 = ml->data+7*pnodecount;
        inst->O0 = ml->data+8*pnodecount;
        inst->O1 = ml->data+9*pnodecount;
        inst->O2 = ml->data+10*pnodecount;
        inst->O3 = ml->data+11*pnodecount;
        inst->O4 = ml->data+12*pnodecount;
        inst->c01 = ml->data+13*pnodecount;
        inst->c12 = ml->data+14*pnodecount;
        inst->c23 = ml->data+15*pnodecount;
        inst->c34 = ml->data+16*pnodecount;
        inst->o01 = ml->data+17*pnodecount;
        inst->o12 = ml->data+18*pnodecount;
        inst->o23 = ml->data+19*pnodecount;
        inst->o34 = ml->data+20*pnodecount;
        inst->f0 = ml->data+21*pnodecount;
        inst->f1 = ml->data+22*pnodecount;
        inst->f2 = ml->data+23*pnodecount;
        inst->f3 = ml->data+24*pnodecount;
        inst->f4 = ml->data+25*pnodecount;
        inst->c10 = ml->data+26*pnodecount;
        inst->c21 = ml->data+27*pnodecount;
        inst->c32 = ml->data+28*pnodecount;
        inst->c43 = ml->data+29*pnodecount;
        inst->o10 = ml->data+30*pnodecount;
        inst->o21 = ml->data+31*pnodecount;
        inst->o32 = ml->data+32*pnodecount;
        inst->o43 = ml->data+33*pnodecount;
        inst->b0 = ml->data+34*pnodecount;
        inst->b1 = ml->data+35*pnodecount;
        inst->b2 = ml->data+36*pnodecount;
        inst->b3 = ml->data+37*pnodecount;
        inst->b4 = ml->data+38*pnodecount;
        inst->cai = ml->data+39*pnodecount;
        inst->ek = ml->data+40*pnodecount;
        inst->DC0 = ml->data+41*pnodecount;
        inst->DC1 = ml->data+42*pnodecount;
        inst->DC2 = ml->data+43*pnodecount;
        inst->DC3 = ml->data+44*pnodecount;
        inst->DC4 = ml->data+45*pnodecount;
        inst->DO0 = ml->data+46*pnodecount;
        inst->DO1 = ml->data+47*pnodecount;
        inst->DO2 = ml->data+48*pnodecount;
        inst->DO3 = ml->data+49*pnodecount;
        inst->DO4 = ml->data+50*pnodecount;
        inst->v_unused = ml->data+51*pnodecount;
        inst->g_unused = ml->data+52*pnodecount;
        inst->ion_ek = nt->_data;
        inst->ion_ik = nt->_data;
        inst->ion_dikdv = nt->_data;
        inst->ion_cai = nt->_data;
        inst->ion_cao = nt->_data;
    }



    static void nrn_alloc_Kca1_1(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);

        #endif
    }


    inline int rates_Kca1_1(int id, int pnodecount, Kca1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double ca);


    struct functor_Kca1_1_1 {
        NrnThread* nt;
        Kca1_1_Instance* inst;
        int id, pnodecount;
        double v;
        const Datum* indexes;
        double* data;
        ThreadDatum* thread;
        double old_C0, old_C1, old_C2, old_C3, old_C4, old_O0, old_O1, old_O2, old_O3;

        void initialize() {
            {
                double qt, alpha, beta, v_in_0, ca_in_0;
                v_in_0 = v;
                ca_in_0 = inst->cai[id];
                qt = pow(inst->global->q10, ((*(inst->celsius) - 23.0) / 10.0));
                inst->c01[id] = 4.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c12[id] = 3.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c23[id] = 2.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c34[id] = 1.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o01[id] = 4.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o12[id] = 3.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o23[id] = 2.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o34[id] = 1.0 * ca_in_0 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c10[id] = 1.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c21[id] = 2.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c32[id] = 3.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c43[id] = 4.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o10[id] = 1.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o21[id] = 2.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o32[id] = 3.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o43[id] = 4.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                alpha = exp(inst->global->Qo * FARADAY * v_in_0 / R / (273.15 + *(inst->celsius)));
                beta = exp(inst->global->Qc * FARADAY * v_in_0 / R / (273.15 + *(inst->celsius)));
                inst->f0[id] = inst->global->pf0 * alpha * qt;
                inst->f1[id] = inst->global->pf1 * alpha * qt;
                inst->f2[id] = inst->global->pf2 * alpha * qt;
                inst->f3[id] = inst->global->pf3 * alpha * qt;
                inst->f4[id] = inst->global->pf4 * alpha * qt;
                inst->b0[id] = inst->global->pb0 * beta * qt;
                inst->b1[id] = inst->global->pb1 * beta * qt;
                inst->b2[id] = inst->global->pb2 * beta * qt;
                inst->b3[id] = inst->global->pb3 * beta * qt;
                inst->b4[id] = inst->global->pb4 * beta * qt;
            }
            old_C0 = inst->C0[id];
            old_C1 = inst->C1[id];
            old_C2 = inst->C2[id];
            old_C3 = inst->C3[id];
            old_C4 = inst->C4[id];
            old_O0 = inst->O0[id];
            old_O1 = inst->O1[id];
            old_O2 = inst->O2[id];
            old_O3 = inst->O3[id];
        }

        functor_Kca1_1_1(NrnThread* nt, Kca1_1_Instance* inst, int id, int pnodecount, double v, const Datum* indexes, double* data, ThreadDatum* thread) : nt{nt}, inst{inst}, id{id}, pnodecount{pnodecount}, v{v}, indexes{indexes}, data{data}, thread{thread} {}
        void operator()(const Eigen::Matrix<double, 10, 1>& nmodl_eigen_xm, Eigen::Matrix<double, 10, 1>& nmodl_eigen_fm, Eigen::Matrix<double, 10, 10>& nmodl_eigen_jm) const {
            const double* nmodl_eigen_x = nmodl_eigen_xm.data();
            double* nmodl_eigen_j = nmodl_eigen_jm.data();
            double* nmodl_eigen_f = nmodl_eigen_fm.data();
            nmodl_eigen_f[static_cast<int>(0)] =  -nmodl_eigen_x[static_cast<int>(0)] * inst->c01[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(0)] * nt->_dt * inst->f0[id] - nmodl_eigen_x[static_cast<int>(0)] + nmodl_eigen_x[static_cast<int>(1)] * inst->c10[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(5)] * inst->b0[id] * nt->_dt + old_C0;
            nmodl_eigen_j[static_cast<int>(0)] =  -inst->c01[id] * nt->_dt - nt->_dt * inst->f0[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(10)] = inst->c10[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(20)] = 0.0;
            nmodl_eigen_j[static_cast<int>(30)] = 0.0;
            nmodl_eigen_j[static_cast<int>(40)] = 0.0;
            nmodl_eigen_j[static_cast<int>(50)] = inst->b0[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(60)] = 0.0;
            nmodl_eigen_j[static_cast<int>(70)] = 0.0;
            nmodl_eigen_j[static_cast<int>(80)] = 0.0;
            nmodl_eigen_j[static_cast<int>(90)] = 0.0;
            nmodl_eigen_f[static_cast<int>(1)] = nmodl_eigen_x[static_cast<int>(0)] * inst->c01[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * inst->c10[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * inst->c12[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->f1[id] - nmodl_eigen_x[static_cast<int>(1)] + nmodl_eigen_x[static_cast<int>(2)] * inst->c21[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(6)] * inst->b1[id] * nt->_dt + old_C1;
            nmodl_eigen_j[static_cast<int>(1)] = inst->c01[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(11)] =  -inst->c10[id] * nt->_dt - inst->c12[id] * nt->_dt - nt->_dt * inst->f1[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(21)] = inst->c21[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(31)] = 0.0;
            nmodl_eigen_j[static_cast<int>(41)] = 0.0;
            nmodl_eigen_j[static_cast<int>(51)] = 0.0;
            nmodl_eigen_j[static_cast<int>(61)] = inst->b1[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(71)] = 0.0;
            nmodl_eigen_j[static_cast<int>(81)] = 0.0;
            nmodl_eigen_j[static_cast<int>(91)] = 0.0;
            nmodl_eigen_f[static_cast<int>(2)] = nmodl_eigen_x[static_cast<int>(1)] * inst->c12[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->c21[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->c23[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->f2[id] - nmodl_eigen_x[static_cast<int>(2)] + nmodl_eigen_x[static_cast<int>(3)] * inst->c32[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(7)] * inst->b2[id] * nt->_dt + old_C2;
            nmodl_eigen_j[static_cast<int>(2)] = 0.0;
            nmodl_eigen_j[static_cast<int>(12)] = inst->c12[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(22)] =  -inst->c21[id] * nt->_dt - inst->c23[id] * nt->_dt - nt->_dt * inst->f2[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(32)] = inst->c32[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(42)] = 0.0;
            nmodl_eigen_j[static_cast<int>(52)] = 0.0;
            nmodl_eigen_j[static_cast<int>(62)] = 0.0;
            nmodl_eigen_j[static_cast<int>(72)] = inst->b2[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(82)] = 0.0;
            nmodl_eigen_j[static_cast<int>(92)] = 0.0;
            nmodl_eigen_f[static_cast<int>(3)] = nmodl_eigen_x[static_cast<int>(2)] * inst->c23[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * inst->c32[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * inst->c34[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->f3[id] - nmodl_eigen_x[static_cast<int>(3)] + nmodl_eigen_x[static_cast<int>(4)] * inst->c43[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(8)] * inst->b3[id] * nt->_dt + old_C3;
            nmodl_eigen_j[static_cast<int>(3)] = 0.0;
            nmodl_eigen_j[static_cast<int>(13)] = 0.0;
            nmodl_eigen_j[static_cast<int>(23)] = inst->c23[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(33)] =  -inst->c32[id] * nt->_dt - inst->c34[id] * nt->_dt - nt->_dt * inst->f3[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(43)] = inst->c43[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(53)] = 0.0;
            nmodl_eigen_j[static_cast<int>(63)] = 0.0;
            nmodl_eigen_j[static_cast<int>(73)] = 0.0;
            nmodl_eigen_j[static_cast<int>(83)] = inst->b3[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(93)] = 0.0;
            nmodl_eigen_f[static_cast<int>(4)] = nmodl_eigen_x[static_cast<int>(3)] * inst->c34[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(4)] * inst->c43[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(4)] * nt->_dt * inst->f4[id] - nmodl_eigen_x[static_cast<int>(4)] + nmodl_eigen_x[static_cast<int>(9)] * inst->b4[id] * nt->_dt + old_C4;
            nmodl_eigen_j[static_cast<int>(4)] = 0.0;
            nmodl_eigen_j[static_cast<int>(14)] = 0.0;
            nmodl_eigen_j[static_cast<int>(24)] = 0.0;
            nmodl_eigen_j[static_cast<int>(34)] = inst->c34[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(44)] =  -inst->c43[id] * nt->_dt - nt->_dt * inst->f4[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(54)] = 0.0;
            nmodl_eigen_j[static_cast<int>(64)] = 0.0;
            nmodl_eigen_j[static_cast<int>(74)] = 0.0;
            nmodl_eigen_j[static_cast<int>(84)] = 0.0;
            nmodl_eigen_j[static_cast<int>(94)] = inst->b4[id] * nt->_dt;
            nmodl_eigen_f[static_cast<int>(5)] = nmodl_eigen_x[static_cast<int>(0)] * nt->_dt * inst->f0[id] - nmodl_eigen_x[static_cast<int>(5)] * inst->b0[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->o01[id] - nmodl_eigen_x[static_cast<int>(5)] + nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o10[id] + old_O0;
            nmodl_eigen_j[static_cast<int>(5)] = nt->_dt * inst->f0[id];
            nmodl_eigen_j[static_cast<int>(15)] = 0.0;
            nmodl_eigen_j[static_cast<int>(25)] = 0.0;
            nmodl_eigen_j[static_cast<int>(35)] = 0.0;
            nmodl_eigen_j[static_cast<int>(45)] = 0.0;
            nmodl_eigen_j[static_cast<int>(55)] =  -inst->b0[id] * nt->_dt - nt->_dt * inst->o01[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(65)] = nt->_dt * inst->o10[id];
            nmodl_eigen_j[static_cast<int>(75)] = 0.0;
            nmodl_eigen_j[static_cast<int>(85)] = 0.0;
            nmodl_eigen_j[static_cast<int>(95)] = 0.0;
            nmodl_eigen_f[static_cast<int>(6)] = nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->f1[id] + nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->o01[id] - nmodl_eigen_x[static_cast<int>(6)] * inst->b1[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o10[id] - nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o12[id] - nmodl_eigen_x[static_cast<int>(6)] + nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o21[id] + old_O1;
            nmodl_eigen_j[static_cast<int>(6)] = 0.0;
            nmodl_eigen_j[static_cast<int>(16)] = nt->_dt * inst->f1[id];
            nmodl_eigen_j[static_cast<int>(26)] = 0.0;
            nmodl_eigen_j[static_cast<int>(36)] = 0.0;
            nmodl_eigen_j[static_cast<int>(46)] = 0.0;
            nmodl_eigen_j[static_cast<int>(56)] = nt->_dt * inst->o01[id];
            nmodl_eigen_j[static_cast<int>(66)] =  -inst->b1[id] * nt->_dt - nt->_dt * inst->o10[id] - nt->_dt * inst->o12[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(76)] = nt->_dt * inst->o21[id];
            nmodl_eigen_j[static_cast<int>(86)] = 0.0;
            nmodl_eigen_j[static_cast<int>(96)] = 0.0;
            nmodl_eigen_f[static_cast<int>(7)] = nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->f2[id] + nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o12[id] - nmodl_eigen_x[static_cast<int>(7)] * inst->b2[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o21[id] - nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o23[id] - nmodl_eigen_x[static_cast<int>(7)] + nmodl_eigen_x[static_cast<int>(8)] * nt->_dt * inst->o32[id] + old_O2;
            nmodl_eigen_j[static_cast<int>(7)] = 0.0;
            nmodl_eigen_j[static_cast<int>(17)] = 0.0;
            nmodl_eigen_j[static_cast<int>(27)] = nt->_dt * inst->f2[id];
            nmodl_eigen_j[static_cast<int>(37)] = 0.0;
            nmodl_eigen_j[static_cast<int>(47)] = 0.0;
            nmodl_eigen_j[static_cast<int>(57)] = 0.0;
            nmodl_eigen_j[static_cast<int>(67)] = nt->_dt * inst->o12[id];
            nmodl_eigen_j[static_cast<int>(77)] =  -inst->b2[id] * nt->_dt - nt->_dt * inst->o21[id] - nt->_dt * inst->o23[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(87)] = nt->_dt * inst->o32[id];
            nmodl_eigen_j[static_cast<int>(97)] = 0.0;
            nmodl_eigen_f[static_cast<int>(8)] = nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->f3[id] + nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o23[id] - nmodl_eigen_x[static_cast<int>(8)] * inst->b3[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(8)] * nt->_dt * inst->o32[id] - nmodl_eigen_x[static_cast<int>(8)] * nt->_dt * inst->o34[id] - nmodl_eigen_x[static_cast<int>(8)] + nmodl_eigen_x[static_cast<int>(9)] * nt->_dt * inst->o43[id] + old_O3;
            nmodl_eigen_j[static_cast<int>(8)] = 0.0;
            nmodl_eigen_j[static_cast<int>(18)] = 0.0;
            nmodl_eigen_j[static_cast<int>(28)] = 0.0;
            nmodl_eigen_j[static_cast<int>(38)] = nt->_dt * inst->f3[id];
            nmodl_eigen_j[static_cast<int>(48)] = 0.0;
            nmodl_eigen_j[static_cast<int>(58)] = 0.0;
            nmodl_eigen_j[static_cast<int>(68)] = 0.0;
            nmodl_eigen_j[static_cast<int>(78)] = nt->_dt * inst->o23[id];
            nmodl_eigen_j[static_cast<int>(88)] =  -inst->b3[id] * nt->_dt - nt->_dt * inst->o32[id] - nt->_dt * inst->o34[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(98)] = nt->_dt * inst->o43[id];
            nmodl_eigen_f[static_cast<int>(9)] =  -nmodl_eigen_x[static_cast<int>(0)] - nmodl_eigen_x[static_cast<int>(1)] - nmodl_eigen_x[static_cast<int>(2)] - nmodl_eigen_x[static_cast<int>(3)] - nmodl_eigen_x[static_cast<int>(4)] - nmodl_eigen_x[static_cast<int>(5)] - nmodl_eigen_x[static_cast<int>(6)] - nmodl_eigen_x[static_cast<int>(7)] - nmodl_eigen_x[static_cast<int>(8)] - nmodl_eigen_x[static_cast<int>(9)] + 1.0;
            nmodl_eigen_j[static_cast<int>(9)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(19)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(29)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(39)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(49)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(59)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(69)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(79)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(89)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(99)] =  -1.0;
        }

        void finalize() {
        }
    };


    struct functor_Kca1_1_0 {
        NrnThread* nt;
        Kca1_1_Instance* inst;
        int id, pnodecount;
        double v;
        const Datum* indexes;
        double* data;
        ThreadDatum* thread;
        double old_C0, old_C1, old_C2, old_C3, old_C4, old_O0, old_O1, old_O2, old_O3;

        void initialize() {
            ;
            {
                double qt, alpha, beta, v_in_1, ca_in_1;
                v_in_1 = v;
                ca_in_1 = inst->cai[id];
                qt = pow(inst->global->q10, ((*(inst->celsius) - 23.0) / 10.0));
                inst->c01[id] = 4.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c12[id] = 3.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c23[id] = 2.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c34[id] = 1.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o01[id] = 4.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o12[id] = 3.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o23[id] = 2.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o34[id] = 1.0 * ca_in_1 * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c10[id] = 1.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c21[id] = 2.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c32[id] = 3.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->c43[id] = 4.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o10[id] = 1.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o21[id] = 2.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o32[id] = 3.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                inst->o43[id] = 4.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
                alpha = exp(inst->global->Qo * FARADAY * v_in_1 / R / (273.15 + *(inst->celsius)));
                beta = exp(inst->global->Qc * FARADAY * v_in_1 / R / (273.15 + *(inst->celsius)));
                inst->f0[id] = inst->global->pf0 * alpha * qt;
                inst->f1[id] = inst->global->pf1 * alpha * qt;
                inst->f2[id] = inst->global->pf2 * alpha * qt;
                inst->f3[id] = inst->global->pf3 * alpha * qt;
                inst->f4[id] = inst->global->pf4 * alpha * qt;
                inst->b0[id] = inst->global->pb0 * beta * qt;
                inst->b1[id] = inst->global->pb1 * beta * qt;
                inst->b2[id] = inst->global->pb2 * beta * qt;
                inst->b3[id] = inst->global->pb3 * beta * qt;
                inst->b4[id] = inst->global->pb4 * beta * qt;
            }
            old_C0 = inst->C0[id];
            old_C1 = inst->C1[id];
            old_C2 = inst->C2[id];
            old_C3 = inst->C3[id];
            old_C4 = inst->C4[id];
            old_O0 = inst->O0[id];
            old_O1 = inst->O1[id];
            old_O2 = inst->O2[id];
            old_O3 = inst->O3[id];
        }

        functor_Kca1_1_0(NrnThread* nt, Kca1_1_Instance* inst, int id, int pnodecount, double v, const Datum* indexes, double* data, ThreadDatum* thread) : nt{nt}, inst{inst}, id{id}, pnodecount{pnodecount}, v{v}, indexes{indexes}, data{data}, thread{thread} {}
        void operator()(const Eigen::Matrix<double, 10, 1>& nmodl_eigen_xm, Eigen::Matrix<double, 10, 1>& nmodl_eigen_fm, Eigen::Matrix<double, 10, 10>& nmodl_eigen_jm) const {
            const double* nmodl_eigen_x = nmodl_eigen_xm.data();
            double* nmodl_eigen_j = nmodl_eigen_jm.data();
            double* nmodl_eigen_f = nmodl_eigen_fm.data();
            nmodl_eigen_f[static_cast<int>(0)] =  -nmodl_eigen_x[static_cast<int>(0)] * inst->c01[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(0)] * nt->_dt * inst->f0[id] - nmodl_eigen_x[static_cast<int>(0)] + nmodl_eigen_x[static_cast<int>(1)] * inst->c10[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(5)] * inst->b0[id] * nt->_dt + old_C0;
            nmodl_eigen_j[static_cast<int>(0)] =  -inst->c01[id] * nt->_dt - nt->_dt * inst->f0[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(10)] = inst->c10[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(20)] = 0.0;
            nmodl_eigen_j[static_cast<int>(30)] = 0.0;
            nmodl_eigen_j[static_cast<int>(40)] = 0.0;
            nmodl_eigen_j[static_cast<int>(50)] = inst->b0[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(60)] = 0.0;
            nmodl_eigen_j[static_cast<int>(70)] = 0.0;
            nmodl_eigen_j[static_cast<int>(80)] = 0.0;
            nmodl_eigen_j[static_cast<int>(90)] = 0.0;
            nmodl_eigen_f[static_cast<int>(1)] = nmodl_eigen_x[static_cast<int>(0)] * inst->c01[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * inst->c10[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * inst->c12[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->f1[id] - nmodl_eigen_x[static_cast<int>(1)] + nmodl_eigen_x[static_cast<int>(2)] * inst->c21[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(6)] * inst->b1[id] * nt->_dt + old_C1;
            nmodl_eigen_j[static_cast<int>(1)] = inst->c01[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(11)] =  -inst->c10[id] * nt->_dt - inst->c12[id] * nt->_dt - nt->_dt * inst->f1[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(21)] = inst->c21[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(31)] = 0.0;
            nmodl_eigen_j[static_cast<int>(41)] = 0.0;
            nmodl_eigen_j[static_cast<int>(51)] = 0.0;
            nmodl_eigen_j[static_cast<int>(61)] = inst->b1[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(71)] = 0.0;
            nmodl_eigen_j[static_cast<int>(81)] = 0.0;
            nmodl_eigen_j[static_cast<int>(91)] = 0.0;
            nmodl_eigen_f[static_cast<int>(2)] = nmodl_eigen_x[static_cast<int>(1)] * inst->c12[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->c21[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * inst->c23[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->f2[id] - nmodl_eigen_x[static_cast<int>(2)] + nmodl_eigen_x[static_cast<int>(3)] * inst->c32[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(7)] * inst->b2[id] * nt->_dt + old_C2;
            nmodl_eigen_j[static_cast<int>(2)] = 0.0;
            nmodl_eigen_j[static_cast<int>(12)] = inst->c12[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(22)] =  -inst->c21[id] * nt->_dt - inst->c23[id] * nt->_dt - nt->_dt * inst->f2[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(32)] = inst->c32[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(42)] = 0.0;
            nmodl_eigen_j[static_cast<int>(52)] = 0.0;
            nmodl_eigen_j[static_cast<int>(62)] = 0.0;
            nmodl_eigen_j[static_cast<int>(72)] = inst->b2[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(82)] = 0.0;
            nmodl_eigen_j[static_cast<int>(92)] = 0.0;
            nmodl_eigen_f[static_cast<int>(3)] = nmodl_eigen_x[static_cast<int>(2)] * inst->c23[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * inst->c32[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * inst->c34[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->f3[id] - nmodl_eigen_x[static_cast<int>(3)] + nmodl_eigen_x[static_cast<int>(4)] * inst->c43[id] * nt->_dt + nmodl_eigen_x[static_cast<int>(8)] * inst->b3[id] * nt->_dt + old_C3;
            nmodl_eigen_j[static_cast<int>(3)] = 0.0;
            nmodl_eigen_j[static_cast<int>(13)] = 0.0;
            nmodl_eigen_j[static_cast<int>(23)] = inst->c23[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(33)] =  -inst->c32[id] * nt->_dt - inst->c34[id] * nt->_dt - nt->_dt * inst->f3[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(43)] = inst->c43[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(53)] = 0.0;
            nmodl_eigen_j[static_cast<int>(63)] = 0.0;
            nmodl_eigen_j[static_cast<int>(73)] = 0.0;
            nmodl_eigen_j[static_cast<int>(83)] = inst->b3[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(93)] = 0.0;
            nmodl_eigen_f[static_cast<int>(4)] = nmodl_eigen_x[static_cast<int>(3)] * inst->c34[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(4)] * inst->c43[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(4)] * nt->_dt * inst->f4[id] - nmodl_eigen_x[static_cast<int>(4)] + nmodl_eigen_x[static_cast<int>(9)] * inst->b4[id] * nt->_dt + old_C4;
            nmodl_eigen_j[static_cast<int>(4)] = 0.0;
            nmodl_eigen_j[static_cast<int>(14)] = 0.0;
            nmodl_eigen_j[static_cast<int>(24)] = 0.0;
            nmodl_eigen_j[static_cast<int>(34)] = inst->c34[id] * nt->_dt;
            nmodl_eigen_j[static_cast<int>(44)] =  -inst->c43[id] * nt->_dt - nt->_dt * inst->f4[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(54)] = 0.0;
            nmodl_eigen_j[static_cast<int>(64)] = 0.0;
            nmodl_eigen_j[static_cast<int>(74)] = 0.0;
            nmodl_eigen_j[static_cast<int>(84)] = 0.0;
            nmodl_eigen_j[static_cast<int>(94)] = inst->b4[id] * nt->_dt;
            nmodl_eigen_f[static_cast<int>(5)] = nmodl_eigen_x[static_cast<int>(0)] * nt->_dt * inst->f0[id] - nmodl_eigen_x[static_cast<int>(5)] * inst->b0[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->o01[id] - nmodl_eigen_x[static_cast<int>(5)] + nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o10[id] + old_O0;
            nmodl_eigen_j[static_cast<int>(5)] = nt->_dt * inst->f0[id];
            nmodl_eigen_j[static_cast<int>(15)] = 0.0;
            nmodl_eigen_j[static_cast<int>(25)] = 0.0;
            nmodl_eigen_j[static_cast<int>(35)] = 0.0;
            nmodl_eigen_j[static_cast<int>(45)] = 0.0;
            nmodl_eigen_j[static_cast<int>(55)] =  -inst->b0[id] * nt->_dt - nt->_dt * inst->o01[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(65)] = nt->_dt * inst->o10[id];
            nmodl_eigen_j[static_cast<int>(75)] = 0.0;
            nmodl_eigen_j[static_cast<int>(85)] = 0.0;
            nmodl_eigen_j[static_cast<int>(95)] = 0.0;
            nmodl_eigen_f[static_cast<int>(6)] = nmodl_eigen_x[static_cast<int>(1)] * nt->_dt * inst->f1[id] + nmodl_eigen_x[static_cast<int>(5)] * nt->_dt * inst->o01[id] - nmodl_eigen_x[static_cast<int>(6)] * inst->b1[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o10[id] - nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o12[id] - nmodl_eigen_x[static_cast<int>(6)] + nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o21[id] + old_O1;
            nmodl_eigen_j[static_cast<int>(6)] = 0.0;
            nmodl_eigen_j[static_cast<int>(16)] = nt->_dt * inst->f1[id];
            nmodl_eigen_j[static_cast<int>(26)] = 0.0;
            nmodl_eigen_j[static_cast<int>(36)] = 0.0;
            nmodl_eigen_j[static_cast<int>(46)] = 0.0;
            nmodl_eigen_j[static_cast<int>(56)] = nt->_dt * inst->o01[id];
            nmodl_eigen_j[static_cast<int>(66)] =  -inst->b1[id] * nt->_dt - nt->_dt * inst->o10[id] - nt->_dt * inst->o12[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(76)] = nt->_dt * inst->o21[id];
            nmodl_eigen_j[static_cast<int>(86)] = 0.0;
            nmodl_eigen_j[static_cast<int>(96)] = 0.0;
            nmodl_eigen_f[static_cast<int>(7)] = nmodl_eigen_x[static_cast<int>(2)] * nt->_dt * inst->f2[id] + nmodl_eigen_x[static_cast<int>(6)] * nt->_dt * inst->o12[id] - nmodl_eigen_x[static_cast<int>(7)] * inst->b2[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o21[id] - nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o23[id] - nmodl_eigen_x[static_cast<int>(7)] + nmodl_eigen_x[static_cast<int>(8)] * nt->_dt * inst->o32[id] + old_O2;
            nmodl_eigen_j[static_cast<int>(7)] = 0.0;
            nmodl_eigen_j[static_cast<int>(17)] = 0.0;
            nmodl_eigen_j[static_cast<int>(27)] = nt->_dt * inst->f2[id];
            nmodl_eigen_j[static_cast<int>(37)] = 0.0;
            nmodl_eigen_j[static_cast<int>(47)] = 0.0;
            nmodl_eigen_j[static_cast<int>(57)] = 0.0;
            nmodl_eigen_j[static_cast<int>(67)] = nt->_dt * inst->o12[id];
            nmodl_eigen_j[static_cast<int>(77)] =  -inst->b2[id] * nt->_dt - nt->_dt * inst->o21[id] - nt->_dt * inst->o23[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(87)] = nt->_dt * inst->o32[id];
            nmodl_eigen_j[static_cast<int>(97)] = 0.0;
            nmodl_eigen_f[static_cast<int>(8)] = nmodl_eigen_x[static_cast<int>(3)] * nt->_dt * inst->f3[id] + nmodl_eigen_x[static_cast<int>(7)] * nt->_dt * inst->o23[id] - nmodl_eigen_x[static_cast<int>(8)] * inst->b3[id] * nt->_dt - nmodl_eigen_x[static_cast<int>(8)] * nt->_dt * inst->o32[id] - nmodl_eigen_x[static_cast<int>(8)] * nt->_dt * inst->o34[id] - nmodl_eigen_x[static_cast<int>(8)] + nmodl_eigen_x[static_cast<int>(9)] * nt->_dt * inst->o43[id] + old_O3;
            nmodl_eigen_j[static_cast<int>(8)] = 0.0;
            nmodl_eigen_j[static_cast<int>(18)] = 0.0;
            nmodl_eigen_j[static_cast<int>(28)] = 0.0;
            nmodl_eigen_j[static_cast<int>(38)] = nt->_dt * inst->f3[id];
            nmodl_eigen_j[static_cast<int>(48)] = 0.0;
            nmodl_eigen_j[static_cast<int>(58)] = 0.0;
            nmodl_eigen_j[static_cast<int>(68)] = 0.0;
            nmodl_eigen_j[static_cast<int>(78)] = nt->_dt * inst->o23[id];
            nmodl_eigen_j[static_cast<int>(88)] =  -inst->b3[id] * nt->_dt - nt->_dt * inst->o32[id] - nt->_dt * inst->o34[id] - 1.0;
            nmodl_eigen_j[static_cast<int>(98)] = nt->_dt * inst->o43[id];
            nmodl_eigen_f[static_cast<int>(9)] =  -nmodl_eigen_x[static_cast<int>(0)] - nmodl_eigen_x[static_cast<int>(1)] - nmodl_eigen_x[static_cast<int>(2)] - nmodl_eigen_x[static_cast<int>(3)] - nmodl_eigen_x[static_cast<int>(4)] - nmodl_eigen_x[static_cast<int>(5)] - nmodl_eigen_x[static_cast<int>(6)] - nmodl_eigen_x[static_cast<int>(7)] - nmodl_eigen_x[static_cast<int>(8)] - nmodl_eigen_x[static_cast<int>(9)] + 1.0;
            nmodl_eigen_j[static_cast<int>(9)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(19)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(29)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(39)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(49)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(59)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(69)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(79)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(89)] =  -1.0;
            nmodl_eigen_j[static_cast<int>(99)] =  -1.0;
        }

        void finalize() {
        }
    };


    inline int rates_Kca1_1(int id, int pnodecount, Kca1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v, double ca) {
        int ret_rates = 0;
        double qt, alpha, beta;
        qt = pow(inst->global->q10, ((*(inst->celsius) - 23.0) / 10.0));
        inst->c01[id] = 4.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c12[id] = 3.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c23[id] = 2.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c34[id] = 1.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o01[id] = 4.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o12[id] = 3.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o23[id] = 2.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o34[id] = 1.0 * ca * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c10[id] = 1.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c21[id] = 2.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c32[id] = 3.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
        inst->c43[id] = 4.0 * inst->global->Kc * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o10[id] = 1.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o21[id] = 2.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o32[id] = 3.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
        inst->o43[id] = 4.0 * inst->global->Ko * inst->global->k1 * inst->global->onoffrate * qt;
        alpha = exp(inst->global->Qo * FARADAY * arg_v / R / (273.15 + *(inst->celsius)));
        beta = exp(inst->global->Qc * FARADAY * arg_v / R / (273.15 + *(inst->celsius)));
        inst->f0[id] = inst->global->pf0 * alpha * qt;
        inst->f1[id] = inst->global->pf1 * alpha * qt;
        inst->f2[id] = inst->global->pf2 * alpha * qt;
        inst->f3[id] = inst->global->pf3 * alpha * qt;
        inst->f4[id] = inst->global->pf4 * alpha * qt;
        inst->b0[id] = inst->global->pb0 * beta * qt;
        inst->b1[id] = inst->global->pb1 * beta * qt;
        inst->b2[id] = inst->global->pb2 * beta * qt;
        inst->b3[id] = inst->global->pb3 * beta * qt;
        inst->b4[id] = inst->global->pb4 * beta * qt;
        return ret_rates;
    }


    /** initialize channel */
    void nrn_init_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);

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
                inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
                inst->cai[id] = inst->ion_cai[indexes[3*pnodecount + id]];
                inst->C0[id] = inst->global->C00;
                inst->C1[id] = inst->global->C10;
                inst->C2[id] = inst->global->C20;
                inst->C3[id] = inst->global->C30;
                inst->C4[id] = inst->global->C40;
                inst->O0[id] = inst->global->O00;
                inst->O1[id] = inst->global->O10;
                inst->O2[id] = inst->global->O20;
                inst->O3[id] = inst->global->O30;
                inst->O4[id] = inst->global->O40;
                                
                Eigen::Matrix<double, 10, 1> nmodl_eigen_xm;
                double* nmodl_eigen_x = nmodl_eigen_xm.data();
                nmodl_eigen_x[static_cast<int>(0)] = inst->C0[id];
                nmodl_eigen_x[static_cast<int>(1)] = inst->C1[id];
                nmodl_eigen_x[static_cast<int>(2)] = inst->C2[id];
                nmodl_eigen_x[static_cast<int>(3)] = inst->C3[id];
                nmodl_eigen_x[static_cast<int>(4)] = inst->C4[id];
                nmodl_eigen_x[static_cast<int>(5)] = inst->O0[id];
                nmodl_eigen_x[static_cast<int>(6)] = inst->O1[id];
                nmodl_eigen_x[static_cast<int>(7)] = inst->O2[id];
                nmodl_eigen_x[static_cast<int>(8)] = inst->O3[id];
                nmodl_eigen_x[static_cast<int>(9)] = inst->O4[id];
                // call newton solver
                functor_Kca1_1_0 newton_functor(nt, inst, id, pnodecount, v, indexes, data, thread);
                newton_functor.initialize();
                int newton_iterations = nmodl::newton::newton_solver(nmodl_eigen_xm, newton_functor);
                if (newton_iterations < 0) assert(false && "Newton solver did not converge!");
                inst->C0[id] = nmodl_eigen_x[static_cast<int>(0)];
                inst->C1[id] = nmodl_eigen_x[static_cast<int>(1)];
                inst->C2[id] = nmodl_eigen_x[static_cast<int>(2)];
                inst->C3[id] = nmodl_eigen_x[static_cast<int>(3)];
                inst->C4[id] = nmodl_eigen_x[static_cast<int>(4)];
                inst->O0[id] = nmodl_eigen_x[static_cast<int>(5)];
                inst->O1[id] = nmodl_eigen_x[static_cast<int>(6)];
                inst->O2[id] = nmodl_eigen_x[static_cast<int>(7)];
                inst->O3[id] = nmodl_eigen_x[static_cast<int>(8)];
                inst->O4[id] = nmodl_eigen_x[static_cast<int>(9)];
                newton_functor.finalize();


            }
            nt->_dt = _save_prev_dt;
        }
    }


    inline double nrn_current_Kca1_1(int id, int pnodecount, Kca1_1_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        inst->g[id] = inst->gbar[id] * (inst->O0[id] + inst->O1[id] + inst->O2[id] + inst->O3[id] + inst->O4[id]);
        inst->ik[id] = inst->g[id] * (v - inst->ek[id]);
        current += inst->ik[id];
        return current;
    }


    /** update current */
    void nrn_cur_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            inst->cai[id] = inst->ion_cai[indexes[3*pnodecount + id]];
            double g = nrn_current_Kca1_1(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double dik = inst->ik[id];
            double rhs = nrn_current_Kca1_1(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Kca1_1(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Kca1_1_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->ek[id] = inst->ion_ek[indexes[0*pnodecount + id]];
            inst->cai[id] = inst->ion_cai[indexes[3*pnodecount + id]];
            
            Eigen::Matrix<double, 10, 1> nmodl_eigen_xm;
            double* nmodl_eigen_x = nmodl_eigen_xm.data();
            nmodl_eigen_x[static_cast<int>(0)] = inst->C0[id];
            nmodl_eigen_x[static_cast<int>(1)] = inst->C1[id];
            nmodl_eigen_x[static_cast<int>(2)] = inst->C2[id];
            nmodl_eigen_x[static_cast<int>(3)] = inst->C3[id];
            nmodl_eigen_x[static_cast<int>(4)] = inst->C4[id];
            nmodl_eigen_x[static_cast<int>(5)] = inst->O0[id];
            nmodl_eigen_x[static_cast<int>(6)] = inst->O1[id];
            nmodl_eigen_x[static_cast<int>(7)] = inst->O2[id];
            nmodl_eigen_x[static_cast<int>(8)] = inst->O3[id];
            nmodl_eigen_x[static_cast<int>(9)] = inst->O4[id];
            // call newton solver
            functor_Kca1_1_1 newton_functor(nt, inst, id, pnodecount, v, indexes, data, thread);
            newton_functor.initialize();
            int newton_iterations = nmodl::newton::newton_solver(nmodl_eigen_xm, newton_functor);
            if (newton_iterations < 0) assert(false && "Newton solver did not converge!");
            inst->C0[id] = nmodl_eigen_x[static_cast<int>(0)];
            inst->C1[id] = nmodl_eigen_x[static_cast<int>(1)];
            inst->C2[id] = nmodl_eigen_x[static_cast<int>(2)];
            inst->C3[id] = nmodl_eigen_x[static_cast<int>(3)];
            inst->C4[id] = nmodl_eigen_x[static_cast<int>(4)];
            inst->O0[id] = nmodl_eigen_x[static_cast<int>(5)];
            inst->O1[id] = nmodl_eigen_x[static_cast<int>(6)];
            inst->O2[id] = nmodl_eigen_x[static_cast<int>(7)];
            inst->O3[id] = nmodl_eigen_x[static_cast<int>(8)];
            inst->O4[id] = nmodl_eigen_x[static_cast<int>(9)];
            newton_functor.finalize();

        }
    }


    /** register channel with the simulator */
    void _Kca11_reg() {

        int mech_type = nrn_get_mechtype("Kca1_1");
        Kca1_1_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_Kca1_1, nrn_cur_Kca1_1, nullptr, nrn_state_Kca1_1, nrn_init_Kca1_1, nrn_private_constructor_Kca1_1, nrn_private_destructor_Kca1_1, first_pointer_var_index(), 1);
        Kca1_1_global.k_type = nrn_get_mechtype("k_ion");
        Kca1_1_global.ca_type = nrn_get_mechtype("ca_ion");

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "k_ion");
        hoc_register_dparam_semantics(mech_type, 1, "k_ion");
        hoc_register_dparam_semantics(mech_type, 2, "k_ion");
        hoc_register_dparam_semantics(mech_type, 3, "ca_ion");
        hoc_register_dparam_semantics(mech_type, 4, "ca_ion");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
